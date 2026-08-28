#!/usr/bin/env python3
"""Precompute the vegetation-loss (NDVI change) overlay for the viewer.

Reads the pre-event (surface reflectance) and post-event PlanetScope (TOA
radiance) analytic COGs plus their UDM2 masks straight from the Source
Cooperative AWS mirror, computes mean clear-sky NDVI per epoch on a common
~24 m EPSG:4326 grid, and writes:

    data/ndvi_change.png   RGBA overlay (red = vegetation lost)
    data/ndvi_change.json  bounds + parameters for the viewer

NDVI from TOA radiance is approximate; combined with the seasonal green-up
between May and August it makes NDVI *loss* a conservative damage signal:
monsoon vegetation should be greener, so a strong drop is meaningful.

Usage:  python3 make_ndvi_change.py   (needs rasterio, pillow, numpy)
"""
import json
import os

import numpy as np
import rasterio
from rasterio.warp import reproject, transform as rtransform, Resampling
from rasterio.transform import from_origin
from PIL import Image

BASE = ("https://s3.us-west-2.amazonaws.com/us-west-2.opendata.source.coop"
        "/planet/disasterdata/nepal-flash-flood-2026-08-26")

PRE = ["20260527_053217_72_254a", "20260527_053219_95_254a",
       "20260527_053221_96_254a", "20260527_053224_18_254a",
       "20260527_053226_41_254a"]
POST = ["20260826_050125_99_255f", "20260826_050128_33_255f",
        "20260826_050130_66_255f", "20260826_050133_00_255f",
        "20260826_050135_34_255f", "20260826_054456_67_251f",
        "20260826_054458_74_251f", "20260826_054500_80_251f",
        "20260826_054502_86_251f"]

# target grid: whole catalog bbox at ~24 m
WEST, SOUTH, EAST, NORTH = 84.894, 27.795, 85.648, 28.659
RES = 0.00022  # degrees, ~24 m
WIDTH = int((EAST - WEST) / RES)
HEIGHT = int((NORTH - SOUTH) / RES)
DST_TRANSFORM = from_origin(WEST, NORTH, RES, RES)
DST_CRS = "EPSG:4326"

# thresholds: pixel counted as vegetation loss when it was vegetated
# pre-event and its NDVI dropped despite the seasonal green-up
PRE_VEG = 0.30
MODERATE = -0.25
SEVERE = -0.45

def scene_ndvi(phase, sid):
    """Return (ndvi, clear) reprojected onto the target grid, or None."""
    aname = "analytic_sr" if phase == "pre-event/planetscope-2026-05-27" else "analytic"
    aurl = f"{BASE}/{phase}/items/{sid}/{sid}_{aname}.tif"
    murl = f"{BASE}/{phase}/items/{sid}/{sid}_udm2.tif"
    with rasterio.open(aurl) as ds:
        # read decimated (uses COG overviews): ~4x coarser than native
        d = 8
        oh, ow = ds.height // d, ds.width // d
        red = ds.read(3, out_shape=(oh, ow)).astype("float32")
        nir = ds.read(4, out_shape=(oh, ow)).astype("float32")
        src_transform = ds.transform * ds.transform.scale(ds.width / ow, ds.height / oh)
        src_crs = ds.crs
    with rasterio.open(murl) as dm:
        clear = dm.read(1, out_shape=(oh, ow))  # UDM2 band 1: 1 = clear

    denom = nir + red
    valid = (denom > 0) & (clear == 1)
    ndvi = np.where(valid, (nir - red) / np.where(denom == 0, 1, denom), np.nan)

    dst_ndvi = np.full((HEIGHT, WIDTH), np.nan, dtype="float32")
    reproject(ndvi, dst_ndvi, src_transform=src_transform, src_crs=src_crs,
              dst_transform=DST_TRANSFORM, dst_crs=DST_CRS,
              src_nodata=np.nan, dst_nodata=np.nan,
              resampling=Resampling.average)
    return dst_ndvi

def epoch_mean(phase, ids):
    total = np.zeros((HEIGHT, WIDTH), dtype="float32")
    count = np.zeros((HEIGHT, WIDTH), dtype="uint8")
    for sid in ids:
        print("  ", sid)
        nd = scene_ndvi(phase, sid)
        good = ~np.isnan(nd)
        total[good] += nd[good]
        count[good] += 1
    mean = np.full((HEIGHT, WIDTH), np.nan, dtype="float32")
    has = count > 0
    mean[has] = total[has] / count[has]
    return mean

def main():
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
    os.makedirs(out_dir, exist_ok=True)

    print("pre-event NDVI")
    pre = epoch_mean("pre-event/planetscope-2026-05-27", PRE)
    print("post-event NDVI")
    post = epoch_mean("post-event/planetscope-2026-08-26", POST)

    both = ~np.isnan(pre) & ~np.isnan(post)
    delta = post - pre
    was_veg = both & (pre >= PRE_VEG)
    moderate = was_veg & (delta <= MODERATE) & (delta > SEVERE)
    severe = was_veg & (delta <= SEVERE)

    rgba = np.zeros((HEIGHT, WIDTH, 4), dtype="uint8")
    rgba[moderate] = [255, 120, 40, 175]   # orange: probable loss
    rgba[severe] = [220, 30, 30, 205]      # red: severe loss
    Image.fromarray(rgba).save(os.path.join(out_dir, "ndvi_change.png"))

    meta = {
        "bounds": [[SOUTH, WEST], [NORTH, EAST]],
        "res_deg": RES,
        "pre_scenes": PRE, "post_scenes": POST,
        "pre_veg_threshold": PRE_VEG,
        "moderate_delta": MODERATE, "severe_delta": SEVERE,
        "compared_px": int(both.sum()),
        "moderate_px": int(moderate.sum()), "severe_px": int(severe.sum()),
        "note": ("Mean clear-sky NDVI (UDM2 band 1) per epoch on a ~24 m grid. "
                 "Pre-event is surface reflectance, post-event TOA radiance: "
                 "treat as a screening layer, not a measurement."),
    }
    with open(os.path.join(out_dir, "ndvi_change.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(json.dumps({k: meta[k] for k in ("compared_px", "moderate_px", "severe_px")}))

if __name__ == "__main__":
    main()
