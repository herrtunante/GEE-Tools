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

// Epoch/sensor groups. Pane z-order stacks sharper or newer imagery on
// top within each epoch; every *Pane listed under before/after is clipped
// by the swipe divider.
var PHASES = {
  pre:          { label: 'Before — PlanetScope, 27 May 2026',        color: '#2ecc71', pane: 'beforePane',   epoch: 'before' },
  pre_vantor:   { label: 'Before — Vantor hi-res, 2021–2024',        color: '#35b8b0', pane: 'beforeHiPane', epoch: 'before' },
  post_ps:      { label: 'After — PlanetScope, 26 Aug 2026',         color: '#ff5d4d', pane: 'afterPane',    epoch: 'after' },
  post_skysat:  { label: 'After — SkySat 0.8 m, 27 Aug 2026',        color: '#ffb14d', pane: 'afterHiPane',  epoch: 'after' },
  post_pelican: { label: 'After — Pelican 0.55 m, 27 Aug 2026',      color: '#c77dff', pane: 'afterHi2Pane', epoch: 'after' },
  post_vantor:  { label: 'After — Vantor Legion ~0.5 m, 27–28 Aug',  color: '#ff6fb5', pane: 'afterHi3Pane', epoch: 'after' }
};
var PANE_Z = {
  beforePane: 350, beforeHiPane: 355,
  afterPane: 360, afterHiPane: 365, afterHi2Pane: 366, afterHi3Pane: 367
};
function phasesOf(epoch) {
  return Object.keys(PHASES).filter(function (k) { return PHASES[k].epoch === epoch; });
}

// Band indices within the 4-band analytic assets. All three products are
// B,G,R,NIR — verified empirically by correlating each pansharpened band
// against the visual RGB asset (the SkySat file's colorinterp tags wrongly
// claim R,G,B).
var BAND_IDX = {
  pre:         { red: 2, green: 1, blue: 0, nir: 3 },
  post_ps:     { red: 2, green: 1, blue: 0, nir: 3 },
  post_skysat: { red: 2, green: 1, blue: 0, nir: 3 }
};

// Rescue-relevant settlements along the corridor (positions checked against
// the pre-event imagery; still verify against imagery before tasking teams).
var POIS = [
  { name: 'Rasuwagadhi',     sub: 'Nepal–China border crossing', lat: 28.2795, lng: 85.3768 },
  { name: 'Timure',          sub: 'village + customs yard',      lat: 28.2545, lng: 85.3640 },
  { name: 'Syabrubesi',      sub: 'trailhead town',              lat: 28.1607, lng: 85.3346 },
  { name: 'Dhunche',         sub: 'Rasuwa district HQ',          lat: 28.1103, lng: 85.2966 },
  { name: 'Betrawati',       sub: 'valley-mouth town',           lat: 27.9683, lng: 85.1862 },
  { name: 'Trishuli Bazaar', sub: 'Bidur municipality',          lat: 27.8990, lng: 85.1472 }
];

var CORRIDOR_BOUNDS = L.latLngBounds([27.79, 84.89], [28.66, 85.65]);
var MAX_ZOOM = 22; // overzoom past native resolution for close inspection

/* ------------------------------------------------------------------ map */

var map = L.map('map', {
  zoomControl: true,
  minZoom: 8,
  maxZoom: MAX_ZOOM,
  center: [28.16, 85.33],
  zoom: 12
});
L.control.scale({ imperial: false }).addTo(map);

Object.keys(PANE_Z).forEach(function (p) { map.createPane(p).style.zIndex = PANE_Z[p]; });

function osmLayer() {
  return L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
    maxZoom: MAX_ZOOM,
    maxNativeZoom: 19,
    attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors' +
      ' &middot; Imagery &copy; <a href="https://www.planet.com/disasterdata/">Planet Labs PBC</a>' +
      ' &amp; &copy; Vantor (Open Data Program), CC-BY-NC-4.0'
  });
}
var basemap = osmLayer().addTo(map);

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

  Object.keys(PANE_Z).forEach(function (p) {
    map.getPane(p).style.clip = p.indexOf('before') === 0 ? beforeClip : afterClip;
  });

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

function setMode(newMode) {
  mode = newMode;
  document.querySelectorAll('.mode-switch button').forEach(function (b) {
    b.classList.toggle('active', b.dataset.mode === newMode);
  });
  applyClips();
  refreshLoupeLayers();
}
document.querySelectorAll('.mode-switch button').forEach(function (btn) {
  btn.addEventListener('click', function () { setMode(btn.dataset.mode); });
});

/* --------------------------------------------------- band combinations */

var bandMode = 'visual';   // 'visual' | 'falsecolor' | 'ndvi'

// NDVI colour ramp: sediment/water blue-greys and browns through to greens.
var NDVI_STOPS = [
  [-1.00, [107, 127, 150]],
  [-0.05, [138, 122,  99]],
  [ 0.10, [201, 185, 138]],
  [ 0.25, [213, 207, 110]],
  [ 0.45, [143, 191,  90]],
  [ 0.70, [ 61, 143,  61]],
  [ 1.00, [ 30, 107,  46]]
];
function ndviColor(v) {
  var i = 1;
  while (i < NDVI_STOPS.length - 1 && NDVI_STOPS[i][0] < v) i++;
  var a = NDVI_STOPS[i - 1], b = NDVI_STOPS[i];
  var t = Math.min(1, Math.max(0, (v - a[0]) / (b[0] - a[0])));
  var c = [0, 1, 2].map(function (k) { return Math.round(a[1][k] + (b[1][k] - a[1][k]) * t); });
  return 'rgb(' + c[0] + ',' + c[1] + ',' + c[2] + ')';
}

// Per-scene contrast stretch for the analytic bands, from a coarse overview
// sample (2nd–98th percentile of non-nodata pixels). Pre-event is surface
// reflectance and post-event is TOA radiance, so a per-scene stretch is the
// only way to get comparable-looking false colour.
var stretchCache = {};

function computeStretch(scene, georaster) {
  if (stretchCache[scene.id]) return Promise.resolve(stretchCache[scene.id]);
  var opts = {
    left: 0, top: 0, right: georaster.width, bottom: georaster.height,
    width: 128, height: 128, resampleMethod: 'nearest'
  };
  return georaster.getValues(opts).then(function (bands) {
    var stretch = bands.map(function (rows) {
      var vals = [];
      for (var r = 0; r < rows.length; r++) {
        var row = rows[r];
        for (var c = 0; c < row.length; c++) {
          if (row[c] > 0) vals.push(row[c]);
        }
      }
      if (!vals.length) return [0, 4000];
      vals.sort(function (a, b) { return a - b; });
      var lo = vals[Math.floor(vals.length * 0.02)];
      var hi = vals[Math.floor(vals.length * 0.98)];
      if (hi <= lo) hi = lo + 1;
      return [lo, hi];
    });
    stretchCache[scene.id] = stretch;
    return stretch;
  }).catch(function (err) {
    console.warn('stretch sampling failed for ' + scene.id + ', using defaults', err);
    var fallback = [[0, 4000], [0, 4000], [0, 4000], [0, 5000]];
    stretchCache[scene.id] = fallback;
    return fallback;
  });
}

function colorFnFor(scene, whichBands, stretch) {
  var idx = BAND_IDX[scene.phase];
  if (whichBands === 'visual') {
    if (scene.ycbcr) {
      // YCbCr-JPEG COGs (Vantor): geotiff.js returns raw Y/Cb/Cr samples
      return function (v) {
        if (!v || (v[0] === 0 && v[1] === 0 && v[2] === 0)) return null;
        var y = v[0], cb = v[1] - 128, cr = v[2] - 128;
        var r = Math.max(0, Math.min(255, Math.round(y + 1.402 * cr)));
        var g = Math.max(0, Math.min(255, Math.round(y - 0.344136 * cb - 0.714136 * cr)));
        var b = Math.max(0, Math.min(255, Math.round(y + 1.772 * cb)));
        if (r === 0 && g === 0 && b === 0) return null;
        return 'rgb(' + r + ',' + g + ',' + b + ')';
      };
    }
    // visual assets are RGBA uint8; alpha 0 marks scene collar
    return function (v) {
      if (!v || v[3] === 0 || (v[0] === 0 && v[1] === 0 && v[2] === 0)) return null;
      return 'rgb(' + v[0] + ',' + v[1] + ',' + v[2] + ')';
    };
  }
  function nodata(v) {
    return !v || (v[0] === 0 && v[1] === 0 && v[2] === 0 && v[3] === 0);
  }
  if (whichBands === 'falsecolor') {
    var sN = stretch[idx.nir], sR = stretch[idx.red], sG = stretch[idx.green];
    // gamma-lifted stretch: on cloudy scenes the upper percentile is set by
    // cloud, so a plain linear stretch crushes terrain into the shadows
    function sc(val, s) {
      var t = Math.max(0, Math.min(1, (val - s[0]) / (s[1] - s[0])));
      return Math.round(Math.pow(t, 0.7) * 255);
    }
    return function (v) {
      if (nodata(v)) return null;
      return 'rgb(' + sc(v[idx.nir], sN) + ',' + sc(v[idx.red], sR) + ',' + sc(v[idx.green], sG) + ')';
    };
  }
  // ndvi
  return function (v) {
    if (nodata(v)) return null;
    var nir = v[idx.nir], red = v[idx.red];
    if (nir + red <= 0) return null;
    return ndviColor((nir - red) / (nir + red));
  };
}

/* ------------------------------------------------------- raster layers */

var statusEl = document.getElementById('status');
var pendingLoads = 0;

function setStatus() {
  statusEl.textContent = pendingLoads > 0
    ? 'Opening ' + pendingLoads + ' scene' + (pendingLoads > 1 ? 's' : '') +
      '… imagery then streams tile-by-tile as you pan and zoom'
    : '';
}

var georasterCache = {};   // url -> Promise<georaster>
function getGeoraster(url) {
  if (!georasterCache[url]) georasterCache[url] = parseGeoraster(url);
  return georasterCache[url];
}

// layer caches: one per (scene, asset) for the main map and the loupe — a
// Leaflet layer instance can only live on one map at a time. False colour
// and NDVI share the analytic layer: switching between them only recolours
// the already-fetched tiles via updateColors(), no re-streaming.
var rasterLayers = {};
var loupeRasterLayers = {};

// Scenes without a 4-band analytic asset (Pelican, Vantor) keep their
// visual rendering in false colour / NDVI mode rather than disappearing.
function assetFor(scene, whichBands) {
  return whichBands === 'visual' || !scene.analytic ? 'visual' : 'analytic';
}

function buildLayer(scene, whichBands, cache, pane) {
  var asset = assetFor(scene, whichBands);
  if (asset === 'visual') whichBands = 'visual';
  cache[scene.id] = cache[scene.id] || {};
  var existing = cache[scene.id][asset];
  if (existing) {
    return existing.then ? existing : Promise.resolve(existing);
  }
  var url = asset === 'visual' ? scene.visual : scene.analytic;
  pendingLoads++; setStatus();
  var promise = getGeoraster(url)
    .then(function (georaster) {
      // stretch stats are needed for false colour; compute them for any
      // analytic layer so later recolours never have to wait on the network
      var needStretch = asset === 'analytic'
        ? computeStretch(scene, georaster)
        : Promise.resolve(null);
      return needStretch.then(function (stretch) {
        // analytic (uint16) reads are heavier than the 8-bit visual COGs, so
        // sample them at a coarser per-tile resolution to keep panning fluid
        var hires = (scene.gsd || 4) < 2;
        var res = asset === 'visual' ? (hires ? 256 : 128) : (hires ? 128 : 64);
        var layer = new GeoRasterLayer({
          georaster: georaster,
          pane: pane,
          resolution: res,
          pixelValuesToColorFn: colorFnFor(scene, whichBands, stretch)
        });
        layer._bands = whichBands;
        cache[scene.id][asset] = layer;
        return layer;
      });
    })
    .catch(function (err) {
      console.error('Failed to open ' + scene.id, err);
      statusEl.textContent = 'Could not open ' + scene.id + ' — check network access to s3.us-west-2.amazonaws.com';
      cache[scene.id][asset] = null;
      return null;
    })
    .finally(function () { pendingLoads--; setStatus(); });
  cache[scene.id][asset] = promise;
  return promise;
}

function retargetLayer(scene, layer, whichBands) {
  // recolour an analytic layer in place when flipping falsecolor <-> ndvi
  if (assetFor(scene, whichBands) === 'visual') whichBands = 'visual';
  if (layer._bands === whichBands) return;
  layer._bands = whichBands;
  layer.updateColors(colorFnFor(scene, whichBands, stretchCache[scene.id]));
}

function showScene(scene) {
  var wantedBands = bandMode;
  // drop a different-asset representation right away so the map never mixes
  // true colour with analytic renders while the new layers stream in
  if (scene._activeLayer && assetFor(scene, scene._activeBands) !== assetFor(scene, wantedBands) &&
      map.hasLayer(scene._activeLayer)) {
    map.removeLayer(scene._activeLayer);
    scene._activeLayer = null;
  }
  buildLayer(scene, wantedBands, rasterLayers, PHASES[scene.phase].pane).then(function (layer) {
    if (!layer) return;
    // still wanted, and the band mode hasn't changed while loading?
    if (scene._wanted && bandMode === wantedBands) {
      if (scene._activeLayer && scene._activeLayer !== layer && map.hasLayer(scene._activeLayer)) {
        map.removeLayer(scene._activeLayer);
      }
      retargetLayer(scene, layer, wantedBands);
      scene._activeLayer = layer;
      scene._activeBands = wantedBands;
      if (!map.hasLayer(layer)) layer.addTo(map);
      applyClips();
    }
  });
}

function setSceneVisible(scene, on) {
  scene._wanted = on;
  if (on) {
    showScene(scene);
  } else if (scene._activeLayer && map.hasLayer(scene._activeLayer)) {
    map.removeLayer(scene._activeLayer);
  }
  refreshLoupeLayers();
}

function setBandMode(newMode) {
  bandMode = newMode;
  document.querySelectorAll('input[name=bands]').forEach(function (r) { r.checked = r.value === newMode; });
  document.getElementById('ndvi-legend').style.display = newMode === 'ndvi' ? '' : 'none';
  SCENES.forEach(function (scene) {
    if (scene._wanted) {
      showScene(scene);
    } else if (scene._activeLayer && map.hasLayer(scene._activeLayer)) {
      map.removeLayer(scene._activeLayer);
    }
  });
  refreshLoupeLayers();
}
document.querySelectorAll('input[name=bands]').forEach(function (r) {
  r.addEventListener('change', function () { if (r.checked) setBandMode(r.value); });
});

/* ------------------------------------------------------- loupe (lens) */

// The loupe follows the cursor and shows the OTHER epoch, magnified: while
// looking at post-event imagery it reveals what stood there before, and
// vice versa.
var loupeOn = false;
var loupeMap = null;
var loupeEl = document.getElementById('loupe');
var loupeLabel = document.getElementById('loupe-label');
var LOUPE_ZOOM_BOOST = 2;

function loupeEpoch() { return mode === 'before' ? 'after' : 'before'; }

function ensureLoupeMap() {
  if (loupeMap) return;
  loupeMap = L.map('loupe-map', {
    zoomControl: false, attributionControl: false,
    dragging: false, scrollWheelZoom: false, doubleClickZoom: false,
    boxZoom: false, keyboard: false, touchZoom: false,
    fadeAnimation: false, zoomAnimation: false, markerZoomAnimation: false, inertia: false,
    minZoom: 8, maxZoom: MAX_ZOOM,
    center: map.getCenter(), zoom: map.getZoom()
  });
  Object.keys(PANE_Z).forEach(function (p) { loupeMap.createPane(p).style.zIndex = PANE_Z[p]; });
  osmLayer().addTo(loupeMap);
}

function refreshLoupeLayers() {
  if (!loupeOn || !loupeMap) return;
  var epoch = loupeEpoch();
  var phases = phasesOf(epoch);
  loupeLabel.textContent = epoch === 'before' ? 'BEFORE THE FLOOD' : 'AFTER · 26–28 AUG';
  loupeLabel.className = epoch === 'before' ? 'lbl-before' : 'lbl-after';
  var wantedBands = bandMode;
  SCENES.forEach(function (scene) {
    var wanted = scene._wanted && phases.indexOf(scene.phase) !== -1;
    if (wanted) {
      buildLayer(scene, wantedBands, loupeRasterLayers, PHASES[scene.phase].pane).then(function (layer) {
        if (!layer) return;
        var stillWanted = loupeOn && bandMode === wantedBands &&
          scene._wanted && phasesOf(loupeEpoch()).indexOf(scene.phase) !== -1;
        if (scene._loupeLayer && scene._loupeLayer !== layer && loupeMap.hasLayer(scene._loupeLayer)) {
          loupeMap.removeLayer(scene._loupeLayer);
        }
        if (stillWanted) {
          retargetLayer(scene, layer, wantedBands);
          scene._loupeLayer = layer;
          if (!loupeMap.hasLayer(layer)) layer.addTo(loupeMap);
        }
      });
    } else if (scene._loupeLayer && loupeMap.hasLayer(scene._loupeLayer)) {
      loupeMap.removeLayer(scene._loupeLayer);
    }
  });
}

var loupeRaf = null;
function moveLoupe(e) {
  if (!loupeOn) return;
  loupeEl.style.display = 'block';
  var pt = e.containerPoint;
  loupeEl.style.left = (pt.x - loupeEl.offsetWidth / 2) + 'px';
  loupeEl.style.top = (pt.y - loupeEl.offsetHeight / 2) + 'px';
  if (loupeRaf) return;
  loupeRaf = requestAnimationFrame(function () {
    loupeRaf = null;
    if (loupeMap) {
      loupeMap.invalidateSize({ animate: false });
      loupeMap.setView(e.latlng, Math.min(MAX_ZOOM, map.getZoom() + LOUPE_ZOOM_BOOST), { animate: false });
    }
  });
}

function setLoupe(on) {
  loupeOn = on;
  document.getElementById('btn-loupe').classList.toggle('active', on);
  if (on) {
    ensureLoupeMap();
    refreshLoupeLayers();
  } else {
    loupeEl.style.display = 'none';
  }
}
map.on('mousemove', moveLoupe);
map.on('mouseout', function () { loupeEl.style.display = 'none'; });
document.getElementById('btn-loupe').addEventListener('click', function () { setLoupe(!loupeOn); });
document.addEventListener('keydown', function (e) {
  if (e.key === 'l' || e.key === 'L') {
    if (/INPUT|TEXTAREA/.test(document.activeElement.tagName)) return;
    setLoupe(!loupeOn);
  }
});

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
  wrap.dataset.phase = phase;
  wrap.innerHTML = '<h3><span class="dot" style="background:' + ph.color + '"></span>' + ph.label + '</h3>';
  SCENES.filter(function (s) { return s.phase === phase; }).forEach(function (scene) {
    var row = document.createElement('div');
    row.className = 'scene';
    var badCloud = scene.cloud >= 60;
    // groups whose scenes share one date show the time; multi-year groups
    // (the Vantor pre-event baseline) show the date instead
    var when = scene.phase === 'pre_vantor'
      ? scene.datetime.slice(0, 10)
      : scene.datetime.slice(11, 16) + ' UTC';
    row.innerHTML =
      '<input type="checkbox" checked>' +
      '<img loading="lazy" src="' + scene.thumbnail + '" alt="">' +
      '<div class="meta">' +
        '<div>' + when + ' ' +
          '<span class="cloudbadge' + (badCloud ? ' bad' : '') + '">' + scene.cloud + '% cloud</span></div>' +
        '<div class="sid">' + scene.id + ' · ' + (scene.gsd || '?') + ' m</div>' +
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
  var g = document.querySelector('.scene-group[data-phase="' + phase + '"]');
  if (g) g.querySelectorAll('input[type=checkbox]').forEach(function (c) { c.checked = on; });
}

// one master toggle per epoch/sensor group
var sourceToggles = document.getElementById('source-toggles');
Object.keys(PHASES).forEach(function (phase) {
  var ph = PHASES[phase];
  var n = SCENES.filter(function (s) { return s.phase === phase; }).length;
  if (!n) return;
  var lab = document.createElement('label');
  lab.className = 'row';
  lab.innerHTML = '<input type="checkbox" checked> <span class="dot" style="background:' +
    ph.color + '"></span>' + ph.label + ' <span class="srcn">(' + n + ')</span>';
  lab.querySelector('input').addEventListener('change', function (e) {
    bulkToggle(phase, e.target.checked);
  });
  sourceToggles.appendChild(lab);
});

/* --------------------------------------------------------------- modal */

var modal = document.getElementById('about');
document.getElementById('btn-about').addEventListener('click', function () { modal.classList.remove('hidden'); });
modal.querySelector('.modal-close').addEventListener('click', function () { modal.classList.add('hidden'); });
modal.addEventListener('click', function (e) { if (e.target === modal) modal.classList.add('hidden'); });

/* ---------------------------------------------------------------- init */

applyClips();
SCENES.forEach(function (scene) { setSceneVisible(scene, true); });
