#!/usr/bin/env python3
"""Regenerate js/scenes.js from the live open-data catalogs.

Sources:
  * Planet Crisis Response release on Source Cooperative (PlanetScope,
    SkySat, Pelican) — collections discovered from the bucket listing, so
    new collections appear automatically.
  * Vantor (formerly Maxar) Open Data Program event bucket (WorldView /
    Legion visual COGs, pre- and post-event).

Both are CC-BY-NC-4.0. Run whenever either dataset gains scenes:

    python3 make_scenes.py
"""
import json
import os
import re
import ssl
import urllib.parse
import urllib.request
from xml.etree import ElementTree

PLANET_BUCKET = "https://s3.us-west-2.amazonaws.com/us-west-2.opendata.source.coop"
PLANET_PREFIX = "planet/disasterdata/nepal-flash-flood-2026-08-26"
VANTOR_BUCKET = "https://vantor-opendata.s3.amazonaws.com"
VANTOR_PREFIX = "events/Nepal-Flooding-Aug-2026"
EVENT_DATE = "2026-08-26"

_ctx = ssl.create_default_context(cafile=os.environ.get("CURL_CA_BUNDLE") or None)

def fetch(url):
    with urllib.request.urlopen(url, context=_ctx, timeout=60) as r:
        return r.read()

def list_keys(bucket, prefix):
    keys, token = [], None
    while True:
        url = f"{bucket}/?list-type=2&prefix={prefix}/&max-keys=1000"
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token)
        root = ElementTree.fromstring(fetch(url))
        ns = {"s3": root.tag.split("}")[0].strip("{")}
        keys += [k.text for k in root.findall(".//s3:Key", ns)]
        token_el = root.find("s3:NextContinuationToken", ns)
        if token_el is None:
            return keys
        token = token_el.text

def rnd(coords):
    if isinstance(coords[0], (int, float)):
        return [round(coords[0], 5), round(coords[1], 5)]
    return [rnd(c) for c in coords]

# phase per Planet collection: sensor name + pre/post position in the tree
PLANET_PHASE = {
    ("pre-event", "planetscope"): "pre",
    ("post-event", "planetscope"): "post_ps",
    ("post-event", "skysat"): "post_skysat",
    ("post-event", "pelican"): "post_pelican",
}
# which asset key carries the 4-band analytic product usable for false
# colour / NDVI (B,G,R,NIR). Pelican's pansharpened is 6-band with an
# unverified layout, so it stays visual-only for now.
PLANET_ANALYTIC = {"pre": "analytic_sr", "post_ps": "analytic", "post_skysat": "pansharpened"}

def planet_scenes():
    keys = list_keys(PLANET_BUCKET, PLANET_PREFIX)
    item_jsons = sorted(k for k in keys if re.search(r"/items/[^/]+/[^/]+\.json$", k))
    scenes = []
    for key in item_jsons:
        m = re.match(rf"{PLANET_PREFIX}/(pre-event|post-event)/([a-z]+)-[0-9-]+/items/", key)
        if not m:
            continue
        phase = PLANET_PHASE.get((m.group(1), m.group(2)))
        if not phase:
            print("  ! unknown Planet collection for", key, "- skipped, extend PLANET_PHASE")
            continue
        it = json.loads(fetch(f"{PLANET_BUCKET}/{key}"))
        p = it["properties"]
        base = f"{PLANET_BUCKET}/{os.path.dirname(key)}"
        sid = it["id"]
        akey = PLANET_ANALYTIC.get(phase)
        scene = {
            "id": sid,
            "phase": phase,
            "provider": "planet",
            "datetime": p["datetime"],
            "cloud": p.get("eo:cloud_cover"),
            "gsd": p.get("gsd"),
            "platform": p.get("platform"),
            "constellation": p.get("constellation"),
            "quality": p.get("pl:quality_category"),
            "bbox": [round(x, 5) for x in it["bbox"]],
            "geometry": {"type": it["geometry"]["type"], "coordinates": rnd(it["geometry"]["coordinates"])},
            "visual": f"{base}/{sid}_visual.tif",
            "thumbnail": f"{base}/{sid}_thumbnail.png",
        }
        if akey and akey in it["assets"]:
            scene["analytic"] = f"{base}/{sid}_{akey}.tif"
            scene["analyticKind"] = akey
        scenes.append(scene)
    return scenes

def vantor_scenes():
    keys = list_keys(VANTOR_BUCKET, VANTOR_PREFIX)
    scenes = []
    for key in sorted(k for k in keys if k.endswith(".json")):
        it = json.loads(fetch(f"{VANTOR_BUCKET}/{key}"))
        if it.get("type") != "Feature" or "properties" not in it or not it.get("bbox"):
            continue  # event/collection documents
        p = it["properties"]
        if not p.get("datetime"):
            continue
        title = p.get("title", "")
        phase = "pre_vantor" if ("[PRE]" in title or p["datetime"][:10] < EVENT_DATE) else "post_vantor"
        sid = it["id"]
        scenes.append({
            "id": sid,
            "phase": phase,
            "provider": "vantor",
            "datetime": p["datetime"],
            "cloud": p.get("eo:cloud_cover"),
            "gsd": p.get("pan_gsd") or p.get("gsd"),
            "platform": p.get("vehicle_name"),
            "constellation": p.get("constellation"),
            "quality": None,
            "bbox": [round(x, 5) for x in it["bbox"]],
            "geometry": {"type": it["geometry"]["type"], "coordinates": rnd(it["geometry"]["coordinates"])},
            "visual": f"{VANTOR_BUCKET}/{VANTOR_PREFIX}/{sid}.tif",
            "thumbnail": f"{VANTOR_BUCKET}/{VANTOR_PREFIX}/{sid}.jpg",
            # Vantor COGs are YCbCr-JPEG; geotiff.js returns raw Y/Cb/Cr
            # samples, so the viewer converts to RGB per pixel
            "ycbcr": True,
        })
    return scenes

def main():
    scenes = planet_scenes() + vantor_scenes()
    order = {"pre": 0, "pre_vantor": 1, "post_ps": 2, "post_skysat": 3, "post_pelican": 4, "post_vantor": 5}
    scenes.sort(key=lambda s: (order.get(s["phase"], 9), s["datetime"]))
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "js", "scenes.js")
    with open(out, "w") as f:
        f.write("// Auto-generated by tools/make_scenes.py from:\n")
        f.write("//   https://source.coop/planet/disasterdata/nepal-flash-flood-2026-08-26\n")
        f.write("//   s3://vantor-opendata/events/Nepal-Flooding-Aug-2026 (Vantor Open Data Program)\n")
        f.write("// Imagery (c) Planet Labs PBC / (c) Vantor, CC-BY-NC-4.0\n")
        f.write("const SCENES = ")
        f.write(json.dumps(scenes, separators=(",", ":")))
        f.write(";\n")
    counts = {}
    for s in scenes:
        counts[s["phase"]] = counts.get(s["phase"], 0) + 1
    print("wrote", len(scenes), "scenes:", counts)

if __name__ == "__main__":
    main()
