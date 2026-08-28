# Nepal Outburst Flood 2026 — Before / After Viewer

Static web viewer comparing satellite imagery from before and after the
26 August 2026 Bhote Koshi–Trishuli outburst flood (Nepal–China border),
built to help direct emergency rescue efforts.

**Live site: <https://herrtunante.github.io/GEE-Tools/>** (served from the
`gh-pages` branch — see [Deploying](#deploying) below).

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
- **Band combinations** — true colour (RGB `visual` assets), false colour
  NIR·R·G and NDVI, both streamed from the 4-band analytic assets
  (`analytic_sr` pre-event, `analytic` post-event, `pansharpened` SkySat)
  with a per-scene 2–98% percentile stretch (gamma-lifted so terrain stays
  readable under the monsoon cloud) sampled from the COG overviews. False
  colour and NDVI share one layer per scene, so flipping between them
  recolours already-fetched tiles instantly. All three products are
  B,G,R,NIR — verified empirically against the visual assets (the SkySat
  file's colorinterp tags wrongly claim R,G,B).
- **Comparison loupe** (button or `L`) — a magnifier that follows the cursor
  and shows the *other* epoch at +2 zoom: hover post-event debris to see
  what stood there before, and vice versa.
- **Deep zoom** — the map overzooms to z22 so the 0.5 m/px SkySat
  pansharpened detail can be inspected at building scale.
- **Rescue-relevant landmarks** — one-click zoom to Rasuwagadhi, Timure,
  Syabrubesi, Dhunche, Betrawati and Trishuli Bazaar.
- **Shareable view links** — the URL hash tracks position, zoom, display
  mode, band combination and swipe position; "Copy view link" hands a
  colleague this exact view.
- **Annotations** — draw markers, lines and areas (hazard / blocked route /
  access confirmed / info) with notes. Stored in the browser
  (localStorage) and exchanged as GeoJSON files via Export/Import, so a
  coordinator can mark up the map and send the file to field teams.
- **Vegetation-loss overlay (ΔNDVI)** — a precomputed screening layer:
  mean clear-sky NDVI drop between 27 May and 26 Aug on a ~24 m grid,
  cloud-masked with UDM2, in two classes (severe / probable). Regenerate
  with `tools/make_ndvi_change.py` (needs `rasterio`, `pillow`) when the
  dataset gains scenes. It is a screening aid — pre-event is surface
  reflectance and post-event TOA radiance, and monsoon agriculture can
  register as change — so cross-check against imagery.
- **Offline support** — a service worker precaches the app shell and keeps
  visited imagery byte-ranges and basemap tiles, so areas already viewed
  keep working on a degraded or absent connection in the field.
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
require an http(s) origin).

## Deploying

The site is published on GitHub Pages from the `gh-pages` branch, which
holds this folder's contents at its root. To redeploy after changing the
viewer, regenerate that branch from the development branch:

```bash
git subtree split --prefix nepal-landslide-viewer -b gh-pages-build
git push -f origin gh-pages-build:gh-pages
git branch -D gh-pages-build
```

Internet access to `s3.us-west-2.amazonaws.com` (imagery) and
`tile.openstreetmap.org` (basemap, optional) is required at view time.

## Data

Imagery comes from two open-data releases, both CC-BY-NC-4.0:

- **Planet Crisis Response** on Source Cooperative
  (`https://source.coop/planet/disasterdata/nepal-flash-flood-2026-08-26`,
  mirrored on AWS Open Data)
- **Vantor (formerly Maxar) Open Data Program**
  (`s3://vantor-opendata/events/Nepal-Flooding-Aug-2026`)

| Epoch | Source / sensor | Scenes | Acquired | GSD |
|---|---|---|---|---|
| Before | PlanetScope | 5 | 2026-05-27 05:32 UTC | ~3.8 m |
| Before | Vantor WorldView | 3 | 2021-10 / 2023-09 / 2024-05 | ~0.5 m |
| After | PlanetScope | 9 | 2026-08-26 05:01 / 05:45 UTC | ~3.8 m |
| After | SkySat | 2 | 2026-08-27 02:00 UTC | ~0.8 m |
| After | Pelican | 3 | 2026-08-27 06:09 UTC | 0.55 m |
| After | Vantor WV/Legion | 9 | 2026-08-27 / 08-28 | ~0.5 m |

`js/scenes.js` is generated by `tools/make_scenes.py`, which walks both
live catalogs (new Planet collections are discovered automatically from
the bucket layout). Rerun it whenever either dataset gains scenes — the
Planet release was restructured once already, which broke hardcoded
paths. Notes: Vantor COGs are YCbCr-JPEG, which geotiff.js returns
unconverted, so the viewer converts Y/Cb/Cr to RGB per pixel; Pelican's
6-band pansharpened layout is unverified, so Pelican and Vantor scenes
stay true-colour in the band-combination modes.

### Interpretation caveats (from the dataset README)

- Pre-event analytic is surface reflectance while post-event is
  top-of-atmosphere radiance: in false colour and NDVI, compare spatial
  patterns across epochs, not absolute values.

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
