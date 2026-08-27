/* Nepal Bhote Koshi–Trishuli outburst flood (26 Aug 2026) before/after viewer.
 *
 * Streams Planet Crisis Response Cloud-Optimized GeoTIFFs straight from
 * Source Cooperative's AWS mirror — the browser fetches only the byte ranges
 * it needs for the current view, so there is no backend to deploy.
 *
 * Imagery (c) Planet Labs PBC, CC-BY-NC-4.0.
 */
'use strict';

/* ---------------------------------------------------------------- config */

var PHASES = {
  pre:         { label: 'Before — PlanetScope, 27 May 2026', color: '#2ecc71', pane: 'beforePane' },
  post_ps:     { label: 'After — PlanetScope, 26 Aug 2026',  color: '#ff5d4d', pane: 'afterPane' },
  post_skysat: { label: 'After — SkySat 0.8 m, 27 Aug 2026', color: '#ffb14d', pane: 'afterHiPane' }
};

// Rescue-relevant settlements along the corridor (approximate coordinates —
// verify precise positions against the imagery itself before tasking teams).
var POIS = [
  { name: 'Rasuwagadhi',     sub: 'Nepal–China border crossing', lat: 28.2795, lng: 85.3768 },
  { name: 'Timure',          sub: 'village + customs yard',      lat: 28.2545, lng: 85.3640 },
  { name: 'Syabrubesi',      sub: 'trailhead town',              lat: 28.1607, lng: 85.3346 },
  { name: 'Dhunche',         sub: 'Rasuwa district HQ',          lat: 28.1103, lng: 85.2966 },
  { name: 'Betrawati',       sub: 'valley-mouth town',           lat: 27.9683, lng: 85.1862 },
  { name: 'Trishuli Bazaar', sub: 'Bidur municipality',          lat: 27.8990, lng: 85.1472 }
];

var CORRIDOR_BOUNDS = L.latLngBounds([27.79, 84.89], [28.66, 85.65]);

/* ------------------------------------------------------------------ map */

var map = L.map('map', {
  zoomControl: true,
  minZoom: 8,
  maxZoom: 20,
  center: [28.16, 85.33],
  zoom: 12
});
L.control.scale({ imperial: false }).addTo(map);

map.createPane('beforePane').style.zIndex = 350;
map.createPane('afterPane').style.zIndex = 360;
map.createPane('afterHiPane').style.zIndex = 365;

var basemap = L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
  maxZoom: 19,
  maxNativeZoom: 19,
  attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors' +
    ' &middot; Imagery &copy; <a href="https://www.planet.com/disasterdata/">Planet Labs PBC</a> CC-BY-NC-4.0' +
    ' via <a href="https://source.coop/planet/disasterdata/nepal-flash-flood-2026-08-26">Source Cooperative</a>'
}).addTo(map);

/* -------------------------------------------------------- swipe control */

var mode = 'swipe';        // 'before' | 'swipe' | 'after'
var swipeFrac = 0.5;       // divider position, 0 = far left
var handle = document.getElementById('swipe-handle');

function applyClips() {
  var nw = map.containerPointToLayerPoint([0, 0]);
  var se = map.containerPointToLayerPoint(map.getSize());
  var clipX = nw.x + (se.x - nw.x) * (mode === 'before' ? 1 : mode === 'after' ? 0 : swipeFrac);

  var beforeClip = 'rect(' + nw.y + 'px ' + clipX + 'px ' + se.y + 'px ' + nw.x + 'px)';
  var afterClip  = 'rect(' + nw.y + 'px ' + se.x  + 'px ' + se.y + 'px ' + clipX + 'px)';

  map.getPane('beforePane').style.clip = beforeClip;
  map.getPane('afterPane').style.clip = afterClip;
  map.getPane('afterHiPane').style.clip = afterClip;

  handle.style.display = (mode === 'swipe') ? '' : 'none';
  handle.style.left = (swipeFrac * 100) + '%';
  document.querySelector('.before-tag').style.display = (mode === 'after') ? 'none' : '';
  document.querySelector('.after-tag').style.display = (mode === 'before') ? 'none' : '';
}

map.on('move zoom viewreset resize', applyClips);

(function makeDraggable() {
  var dragging = false;
  handle.addEventListener('pointerdown', function (e) {
    dragging = true;
    handle.setPointerCapture(e.pointerId);
    e.preventDefault();
  });
  handle.addEventListener('pointermove', function (e) {
    if (!dragging) return;
    var rect = document.getElementById('map-area').getBoundingClientRect();
    swipeFrac = Math.min(0.98, Math.max(0.02, (e.clientX - rect.left) / rect.width));
    applyClips();
  });
  handle.addEventListener('pointerup', function () { dragging = false; });
})();

document.querySelectorAll('.mode-switch button').forEach(function (btn) {
  btn.addEventListener('click', function () {
    mode = btn.dataset.mode;
    document.querySelectorAll('.mode-switch button').forEach(function (b) {
      b.classList.toggle('active', b === btn);
    });
    applyClips();
  });
});

/* ------------------------------------------------------- raster layers */

var statusEl = document.getElementById('status');
var pendingLoads = 0;

function setStatus() {
  statusEl.textContent = pendingLoads > 0
    ? 'Opening ' + pendingLoads + ' scene' + (pendingLoads > 1 ? 's' : '') +
      '… imagery then streams tile-by-tile as you pan and zoom'
    : '';
}

var rasterLayers = {}; // scene id -> GeoRasterLayer (or Promise placeholder)

function sceneLayer(scene) {
  if (rasterLayers[scene.id]) return Promise.resolve(rasterLayers[scene.id]);
  pendingLoads++; setStatus();
  return parseGeoraster(scene.visual)
    .then(function (georaster) {
      var layer = new GeoRasterLayer({
        georaster: georaster,
        pane: PHASES[scene.phase].pane,
        resolution: scene.phase === 'post_skysat' ? 256 : 128,
        pixelValuesToColorFn: function (v) {
          // visual assets are RGBA uint8; alpha 0 marks scene collar
          if (!v || v[3] === 0 || (v[0] === 0 && v[1] === 0 && v[2] === 0)) return null;
          return 'rgb(' + v[0] + ',' + v[1] + ',' + v[2] + ')';
        }
      });
      rasterLayers[scene.id] = layer;
      return layer;
    })
    .catch(function (err) {
      console.error('Failed to open ' + scene.id, err);
      statusEl.textContent = 'Could not open ' + scene.id + ' — check network access to s3.us-west-2.amazonaws.com';
      return null;
    })
    .finally(function () { pendingLoads--; setStatus(); });
}

function setSceneVisible(scene, on) {
  scene._wanted = on;
  if (on) {
    sceneLayer(scene).then(function (layer) {
      if (layer && scene._wanted && !map.hasLayer(layer)) layer.addTo(map);
      applyClips();
    });
  } else if (rasterLayers[scene.id] && map.hasLayer(rasterLayers[scene.id])) {
    map.removeLayer(rasterLayers[scene.id]);
  }
}

/* ---------------------------------------------------------- footprints */

var footprints = L.layerGroup().addTo(map);
SCENES.forEach(function (scene) {
  var ph = PHASES[scene.phase];
  var gj = L.geoJSON({ type: 'Feature', geometry: scene.geometry }, {
    style: { color: ph.color, weight: 1.5, opacity: 0.9, fillOpacity: 0 }
  }).bindPopup(
    '<b>' + scene.id + '</b><br>' + ph.label +
    '<br>Acquired: ' + scene.datetime.replace('T', ' ').slice(0, 16) + ' UTC' +
    '<br>Cloud cover: ' + scene.cloud + '% &middot; GSD ' + scene.gsd + ' m' +
    '<br>Platform: ' + scene.platform + ' (' + scene.constellation + ')'
  );
  gj._sceneId = scene.id;
  footprints.addLayer(gj);
});

/* ---------------------------------------------------------------- POIs */

var poiLayer = L.layerGroup().addTo(map);
POIS.forEach(function (p) {
  L.circleMarker([p.lat, p.lng], {
    radius: 6, color: '#fff', weight: 2, fillColor: '#3fa7ff', fillOpacity: 1
  })
    .bindTooltip(p.name, { permanent: true, direction: 'right', offset: [8, 0], className: 'poi-label' })
    .bindPopup('<b>' + p.name + '</b><br>' + p.sub + '<br><i>Marker position approximate — verify against imagery.</i>')
    .addTo(poiLayer);
});

var poiButtons = document.getElementById('poi-list');
POIS.forEach(function (p) {
  var b = document.createElement('button');
  b.innerHTML = p.name + '<span class="sub">' + p.sub + '</span>';
  b.addEventListener('click', function () { map.setView([p.lat, p.lng], 15); });
  poiButtons.appendChild(b);
});
var bAll = document.createElement('button');
bAll.innerHTML = 'Whole corridor<span class="sub">fit all imagery</span>';
bAll.addEventListener('click', function () { map.fitBounds(CORRIDOR_BOUNDS); });
poiButtons.appendChild(bAll);

/* ------------------------------------------------------- sidebar scenes */

var groupsEl = document.getElementById('scene-groups');
Object.keys(PHASES).forEach(function (phase) {
  var ph = PHASES[phase];
  var wrap = document.createElement('div');
  wrap.className = 'scene-group';
  wrap.innerHTML = '<h3><span class="dot" style="background:' + ph.color + '"></span>' + ph.label + '</h3>';
  SCENES.filter(function (s) { return s.phase === phase; }).forEach(function (scene) {
    var row = document.createElement('div');
    row.className = 'scene';
    var badCloud = scene.cloud >= 60;
    row.innerHTML =
      '<input type="checkbox" checked>' +
      '<img loading="lazy" src="' + scene.thumbnail + '" alt="">' +
      '<div class="meta">' +
        '<div>' + scene.datetime.slice(11, 16) + ' UTC ' +
          '<span class="cloudbadge' + (badCloud ? ' bad' : '') + '">' + scene.cloud + '% cloud</span></div>' +
        '<div class="sid">' + scene.id + '</div>' +
      '</div>' +
      '<button class="zoom-btn" title="Zoom to this scene">⌖</button>';
    row.querySelector('input').addEventListener('change', function (e) {
      setSceneVisible(scene, e.target.checked);
    });
    row.querySelector('.zoom-btn').addEventListener('click', function () {
      map.fitBounds([[scene.bbox[1], scene.bbox[0]], [scene.bbox[3], scene.bbox[2]]]);
    });
    groupsEl.appendChild(wrap);
    wrap.appendChild(row);
  });
});

/* ------------------------------------------------------ overlay toggles */

function bindToggle(id, fn) {
  document.getElementById(id).addEventListener('change', function (e) { fn(e.target.checked); });
}
bindToggle('chk-footprints', function (on) { on ? footprints.addTo(map) : map.removeLayer(footprints); });
bindToggle('chk-pois', function (on) { on ? poiLayer.addTo(map) : map.removeLayer(poiLayer); });
bindToggle('chk-basemap', function (on) { on ? basemap.addTo(map) : map.removeLayer(basemap); });

function bulkToggle(phase, on) {
  SCENES.filter(function (s) { return s.phase === phase; }).forEach(function (s) { setSceneVisible(s, on); });
  // keep the per-scene checkboxes in sync
  document.querySelectorAll('.scene-group').forEach(function (g) {
    if (g.querySelector('h3').textContent.indexOf(PHASES[phase].label) !== -1) {
      g.querySelectorAll('input[type=checkbox]').forEach(function (c) { c.checked = on; });
    }
  });
}
bindToggle('chk-after-ps', function (on) { bulkToggle('post_ps', on); });
bindToggle('chk-after-sky', function (on) { bulkToggle('post_skysat', on); });

/* --------------------------------------------------------------- modal */

var modal = document.getElementById('about');
document.getElementById('btn-about').addEventListener('click', function () { modal.classList.remove('hidden'); });
modal.querySelector('.modal-close').addEventListener('click', function () { modal.classList.add('hidden'); });
modal.addEventListener('click', function (e) { if (e.target === modal) modal.classList.add('hidden'); });

/* ---------------------------------------------------------------- init */

applyClips();
SCENES.forEach(function (scene) { setSceneVisible(scene, true); });
