/* Response tools: shareable view links, annotation layer with GeoJSON
 * exchange, the precomputed vegetation-loss overlay, and offline support.
 * Loads after app.js and uses its globals (map, mode, bandMode, ...).
 */
'use strict';

/* ------------------------------------------------- shareable view links */

// #v=<zoom>/<lat>/<lng>/<mode>/<bands>/<swipe>
function readHash() {
  var m = location.hash.match(/^#v=([\d.]+)\/(-?[\d.]+)\/(-?[\d.]+)\/(\w+)\/(\w+)\/([\d.]+)/);
  if (!m) return null;
  return {
    zoom: parseFloat(m[1]), lat: parseFloat(m[2]), lng: parseFloat(m[3]),
    mode: m[4], bands: m[5], swipe: parseFloat(m[6])
  };
}

var hashTimer = null;
function writeHash() {
  if (hashTimer) return;
  hashTimer = setTimeout(function () {
    hashTimer = null;
    var c = map.getCenter();
    var h = '#v=' + map.getZoom().toFixed(2) + '/' + c.lat.toFixed(5) + '/' +
      c.lng.toFixed(5) + '/' + mode + '/' + bandMode + '/' + swipeFrac.toFixed(2);
    history.replaceState(null, '', h);
  }, 250);
}

(function initHash() {
  var st = readHash();
  if (st) {
    map.setView([st.lat, st.lng], st.zoom);
    if (st.mode !== mode && /^(before|swipe|after)$/.test(st.mode)) setMode(st.mode);
    if (st.bands !== bandMode && /^(visual|falsecolor|ndvi)$/.test(st.bands)) setBandMode(st.bands);
    swipeFrac = Math.min(0.98, Math.max(0.02, st.swipe || 0.5));
    applyClips();
  }
  map.on('moveend zoomend', writeHash);
  document.querySelectorAll('.mode-switch button, input[name=bands]').forEach(function (el) {
    el.addEventListener('click', writeHash);
  });
  handle.addEventListener('pointerup', writeHash);
})();

document.getElementById('btn-share').addEventListener('click', function () {
  writeHash();
  var btn = this;
  setTimeout(function () {
    var url = location.href;
    function done(ok) {
      btn.textContent = ok ? '✓ Link copied' : url;
      setTimeout(function () { btn.textContent = '🔗 Copy view link'; }, 2500);
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(url).then(function () { done(true); }, function () { done(false); });
    } else { window.prompt('Copy this link:', url); done(true); }
  }, 300);
});

/* -------------------------------------------------------- annotations */

var CATS = {
  hazard:  { label: 'Hazard / damage',   color: '#e04545' },
  blocked: { label: 'Blocked route',     color: '#f0952f' },
  access:  { label: 'Access confirmed',  color: '#2ecc71' },
  info:    { label: 'Info / other',      color: '#3fa7ff' }
};
var LS_KEY = 'nlv-annotations';

var annotations = new L.FeatureGroup().addTo(map);

function catOf(layer) {
  return (layer.feature && layer.feature.properties && layer.feature.properties.cat) || 'info';
}

function styleLayer(layer) {
  var c = (CATS[catOf(layer)] || CATS.info).color;
  if (layer.setStyle) layer.setStyle({ color: c, weight: 3, fillOpacity: 0.15 });
}

function popupHtml(layer) {
  var p = layer.feature.properties;
  var cat = CATS[p.cat] || CATS.info;
  return '<span style="color:' + cat.color + '">&#9679;</span> <b>' + cat.label + '</b>' +
    (p.note ? '<br>' + String(p.note).replace(/</g, '&lt;') : '') +
    (p.ts ? '<br><i style="font-size:10px">' + p.ts.slice(0, 16).replace('T', ' ') + ' UTC</i>' : '') +
    '<br><button class="annot-del">Delete</button>';
}

function wireLayer(layer) {
  styleLayer(layer);
  layer.bindPopup(popupHtml(layer));
  layer.on('popupopen', function (e) {
    var btn = e.popup.getElement().querySelector('.annot-del');
    if (btn) btn.onclick = function () {
      annotations.removeLayer(layer);
      saveAnnotations();
    };
  });
}

function addFeature(geojsonFeature) {
  L.geoJSON(geojsonFeature, {
    onEachFeature: function (feature, layer) {
      layer.feature = feature;
      wireLayer(layer);
      annotations.addLayer(layer);
    }
  });
}

function annotationsGeoJSON() {
  var feats = [];
  annotations.eachLayer(function (layer) {
    var gj = layer.toGeoJSON();
    gj.properties = (layer.feature && layer.feature.properties) || {};
    feats.push(gj);
  });
  return { type: 'FeatureCollection', features: feats };
}

function saveAnnotations() {
  try { localStorage.setItem(LS_KEY, JSON.stringify(annotationsGeoJSON())); } catch (e) {}
  var n = annotations.getLayers().length;
  document.getElementById('annot-count').textContent =
    n ? n + ' annotation' + (n > 1 ? 's' : '') + ' (saved in this browser)' : 'No annotations yet';
}

(function loadAnnotations() {
  try {
    var raw = localStorage.getItem(LS_KEY);
    if (raw) JSON.parse(raw).features.forEach(addFeature);
  } catch (e) { console.warn('could not restore annotations', e); }
  saveAnnotations();
})();

var drawControl = new L.Control.Draw({
  position: 'topleft',
  draw: {
    polygon: { shapeOptions: { weight: 3 } },
    polyline: { shapeOptions: { weight: 3 } },
    marker: true,
    circle: false, rectangle: false, circlemarker: false
  },
  edit: { featureGroup: annotations, remove: true }
});
map.addControl(drawControl);

map.on(L.Draw.Event.CREATED, function (e) {
  var cat = document.getElementById('annot-cat').value;
  var note = window.prompt('Note for this ' + (CATS[cat] || CATS.info).label.toLowerCase() +
    ' annotation (optional):', '') || '';
  var layer = e.layer;
  var gj = layer.toGeoJSON();
  gj.properties = { cat: cat, note: note, ts: new Date().toISOString() };
  layer.feature = gj;
  wireLayer(layer);
  annotations.addLayer(layer);
  saveAnnotations();
});
map.on(L.Draw.Event.EDITED, saveAnnotations);
map.on(L.Draw.Event.DELETED, saveAnnotations);

document.getElementById('btn-annot-export').addEventListener('click', function () {
  var blob = new Blob([JSON.stringify(annotationsGeoJSON(), null, 1)],
    { type: 'application/geo+json' });
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = 'nepal-flood-annotations-' + new Date().toISOString().slice(0, 10) + '.geojson';
  a.click();
  URL.revokeObjectURL(a.href);
});

document.getElementById('annot-import').addEventListener('change', function (e) {
  var file = e.target.files[0];
  if (!file) return;
  var reader = new FileReader();
  reader.onload = function () {
    try {
      var gj = JSON.parse(reader.result);
      (gj.features || []).forEach(addFeature);
      saveAnnotations();
    } catch (err) { alert('Could not read that file as GeoJSON: ' + err.message); }
    e.target.value = '';
  };
  reader.readAsText(file);
});

document.getElementById('btn-annot-clear').addEventListener('click', function () {
  if (annotations.getLayers().length &&
      window.confirm('Delete ALL annotations in this browser? Export first if you need them.')) {
    annotations.clearLayers();
    saveAnnotations();
  }
});

/* ------------------------------------------- vegetation-loss overlay */

map.createPane('changePane').style.zIndex = 368;
var changeOverlay = null;

document.getElementById('chk-change').addEventListener('change', function (e) {
  if (e.target.checked) {
    if (changeOverlay) { changeOverlay.addTo(map); return; }
    fetch('data/ndvi_change.json')
      .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then(function (meta) {
        changeOverlay = L.imageOverlay('data/ndvi_change.png', meta.bounds, {
          opacity: 0.8, pane: 'changePane', interactive: false
        }).addTo(map);
      })
      .catch(function (err) {
        console.error('vegetation-loss overlay unavailable', err);
        e.target.checked = false;
        alert('Vegetation-loss overlay data not found (data/ndvi_change.png).');
      });
  } else if (changeOverlay) {
    map.removeLayer(changeOverlay);
  }
});

/* ------------------------------------------------------- offline (SW) */

if ('serviceWorker' in navigator && /^https?:$/.test(location.protocol)) {
  navigator.serviceWorker.register('sw.js').catch(function (err) {
    console.warn('service worker registration failed', err);
  });
}
