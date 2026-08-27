# Nepal Outburst Flood 2026 — Before / After Viewer

Static web viewer comparing satellite imagery from before and after the
26 August 2026 Bhote Koshi–Trishuli outburst flood (Nepal–China border),
built to help direct emergency rescue efforts.

It streams Planet Crisis Response **Cloud-Optimized GeoTIFFs** directly from
[Source Cooperative](https://source.coop/planet/disasterdata/nepal-flash-flood-2026-08-26)'s
AWS Open Data mirror using HTTP range requests — the browser only downloads
the pixels needed for the current view, so **no backend or tile server is
required**. Any static host works.

## Features

- **Swipe comparison** — drag the divider to compare pre-event (27 May 2026)
  and post-event (26–27 Aug 2026) imagery; `Before` / `After` buttons show a
  single epoch full-screen.
- **Two post-event sources** — PlanetScope (~3.8 m, morning of the flood) and
  SkySat (~0.8 m, next morning) over the Rasuwagadhi border crossing and
  Syabrubesi; SkySat draws on top where available and resolves individual
  buildings, bridges and road cuts.
- **Rescue-relevant landmarks** — one-click zoom to Rasuwagadhi, Timure,
  Syabrubesi, Dhunche, Betrawati and Trishuli Bazaar.
- **Scene footprints & metadata** — acquisition time, cloud cover, GSD and
  platform per scene; per-scene visibility toggles with thumbnails.
- **Caveat surfacing** — monsoon cloud cover, non-coincident footprints and
  Planet `test` quality flags are shown prominently, since misreading the
  imagery could misdirect field teams.

## Running

Serve the directory with any static file server, e.g.:

```bash
cd nepal-landslide-viewer
python3 -m http.server 8000
# open http://localhost:8000
```

Opening `index.html` via `file://` will not work (range-request fetches
require an http(s) origin). For GitHub Pages, publish the repository and
point Pages at this folder (or copy it to the Pages root).

Internet access to `s3.us-west-2.amazonaws.com` (imagery) and
`tile.openstreetmap.org` (basemap, optional) is required at view time.

## Data

All imagery comes from the STAC catalog at
`https://source.coop/planet/disasterdata/nepal-flash-flood-2026-08-26`
(mirrored at
`https://s3.us-west-2.amazonaws.com/us-west-2.opendata.source.coop/planet/disasterdata/nepal-flash-flood-2026-08-26/`):

| Epoch | Sensor | Scenes | Acquired | GSD |
|---|---|---|---|---|
| Before | PlanetScope | 5 | 2026-05-27 05:32 UTC | ~3.8 m |
| After | PlanetScope | 9 | 2026-08-26 05:01 / 05:45 UTC | ~3.8 m |
| After | SkySat | 2 | 2026-08-27 02:00 UTC | ~0.8 m |

`js/scenes.js` is generated from the dataset's STAC items (footprints,
cloud cover, asset URLs). Regenerate it against the live catalog if the
dataset is updated.

### Interpretation caveats (from the dataset README)

- Post-event PlanetScope scenes are 62–93% cloud (monsoon); much of the loss
  is thin haze through which terrain remains interpretable.
- The post-event swaths extend beyond the single pre-event strip; some
  post-event area has no baseline. SkySat covers two focal points only.
- 13 of 14 PlanetScope scenes are `quality_category: test` in Planet's API.
- This viewer is an aid, not ground truth — verify road/bridge status
  locally before routing rescue convoys.

## License & attribution

- Imagery © Planet Labs PBC via the
  [Planet Crisis Response Program](https://www.planet.com/disasterdata/),
  licensed [CC-BY-NC-4.0](https://creativecommons.org/licenses/by-nc/4.0/)
  (attribution required, non-commercial use only).
- Basemap © [OpenStreetMap](https://www.openstreetmap.org/copyright) contributors.
- Vendored libraries: [Leaflet](https://leafletjs.com) (BSD-2),
  [georaster](https://github.com/geotiff/georaster) and
  [georaster-layer-for-leaflet](https://github.com/geotiff/georaster-layer-for-leaflet) (MIT),
  which bundle [geotiff.js](https://geotiffjs.github.io/) (MIT).
