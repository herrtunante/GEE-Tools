/* Offline support: precaches the app shell and keeps a runtime cache of
 * imagery byte-ranges and basemap tiles, so areas already viewed keep
 * working on a degraded or absent connection in the field.
 */
'use strict';

var SHELL_CACHE = 'nlv-shell-v3';
var IMG_CACHE = 'nlv-imagery-v1';
var TILE_CACHE = 'nlv-tiles-v1';
var IMG_CACHE_MAX = 4000;

var SHELL = [
  './',
  'index.html',
  'css/style.css',
  'js/scenes.js',
  'js/app.js',
  'js/tools.js',
  'vendor/leaflet/leaflet.js',
  'vendor/leaflet/leaflet.css',
  'vendor/leaflet-draw/leaflet.draw.js',
  'vendor/leaflet-draw/leaflet.draw.css',
  'vendor/georaster.browser.bundle.min.js',
  'vendor/0.georaster.browser.bundle.min.worker.js',
  'vendor/georaster-layer-for-leaflet.min.js',
  'data/ndvi_change.json',
  'data/ndvi_change.png'
];

self.addEventListener('install', function (e) {
  e.waitUntil(
    caches.open(SHELL_CACHE).then(function (cache) {
      // best-effort: a missing optional asset must not break install
      return Promise.all(SHELL.map(function (url) {
        return cache.add(url).catch(function () {});
      }));
    }).then(function () { return self.skipWaiting(); })
  );
});

self.addEventListener('activate', function (e) {
  e.waitUntil(
    caches.keys().then(function (keys) {
      return Promise.all(keys.map(function (k) {
        if ([SHELL_CACHE, IMG_CACHE, TILE_CACHE].indexOf(k) === -1) return caches.delete(k);
      }));
    }).then(function () { return self.clients.claim(); })
  );
});

// Range requests all share one URL, so cache under a synthetic URL that
// includes the byte range.
function rangeKey(request) {
  var range = request.headers.get('range') || 'full';
  return new Request(request.url + (request.url.indexOf('?') === -1 ? '?' : '&') +
    '__swr=' + encodeURIComponent(range));
}

function trim(cacheName, max) {
  caches.open(cacheName).then(function (cache) {
    cache.keys().then(function (keys) {
      if (keys.length > max) {
        keys.slice(0, Math.ceil(max / 10)).forEach(function (k) { cache.delete(k); });
      }
    });
  });
}

self.addEventListener('fetch', function (e) {
  var req = e.request;
  if (req.method !== 'GET') return;
  var url = new URL(req.url);

  // Imagery from the S3 mirror: network first, range-keyed cache fallback.
  // The Cache API refuses 206 responses, so partials are stored as 200
  // with the real status kept in a header and reconstructed on the way out.
  if (url.hostname.endsWith('amazonaws.com')) {
    var key = rangeKey(req);
    e.respondWith(
      fetch(req).then(function (resp) {
        if (resp.ok || resp.status === 206) {
          var clone = resp.clone();
          e.waitUntil(clone.arrayBuffer().then(function (buf) {
            var headers = new Headers();
            ['content-type', 'content-range', 'content-length', 'accept-ranges', 'etag']
              .forEach(function (h) {
                var v = resp.headers.get(h);
                if (v) headers.set(h, v);
              });
            headers.set('x-sw-status', String(resp.status));
            return caches.open(IMG_CACHE).then(function (cache) {
              return cache.put(key, new Response(buf, { status: 200, headers: headers }));
            });
          }).then(function () { trim(IMG_CACHE, IMG_CACHE_MAX); }));
        }
        return resp;
      }).catch(function () {
        return caches.match(key).then(function (hit) {
          if (!hit) return Response.error();
          var status = parseInt(hit.headers.get('x-sw-status') || '200', 10);
          return hit.arrayBuffer().then(function (buf) {
            return new Response(buf, { status: status, headers: hit.headers });
          });
        });
      })
    );
    return;
  }

  // basemap tiles: cache first, refresh in background
  if (url.hostname === 'tile.openstreetmap.org') {
    e.respondWith(
      caches.match(req).then(function (hit) {
        var refresh = fetch(req).then(function (resp) {
          if (resp.ok) {
            var clone = resp.clone();
            caches.open(TILE_CACHE).then(function (cache) { cache.put(req, clone); });
          }
          return resp;
        }).catch(function () { return hit || Response.error(); });
        return hit || refresh;
      })
    );
    return;
  }

  // same-origin app shell: cache first, refresh in background
  if (url.origin === location.origin) {
    e.respondWith(
      caches.match(req, { ignoreSearch: true }).then(function (hit) {
        var refresh = fetch(req).then(function (resp) {
          if (resp.ok) {
            var clone = resp.clone();
            caches.open(SHELL_CACHE).then(function (cache) { cache.put(req, clone); });
          }
          return resp;
        }).catch(function () { return hit || Response.error(); });
        return hit || refresh;
      })
    );
  }
});
