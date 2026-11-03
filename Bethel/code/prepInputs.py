"""Prepare AvaFrame inputs for Bethel: source DEM -> 1 m .asc, release polygon -> shapefile.

Usage (from /home/ron/AvaFrame):
    uv run python Bethel/code/prepInputs.py              # DEM and release polygon
    uv run python Bethel/code/prepInputs.py --dem-only   # DEM only; keeps Inputs/REL (e.g. the adjustRelease.py polygon)

The source DEM (bethel.SOURCE_DEM) has 4 cm cells. It is averaged to a 1 m grid aligned to whole metres; AvaFrame
remeshes that to its computational mesh (5 m) anyway. Cells outside the surveyed area stay nodata, so the
simulation domain is the surveyed area.
"""
import argparse
import math

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.warp import reproject

import bethel as b
base = b.BASE
relThickness = 1.0  # m, PLACEHOLDER release thickness - set to your design value
cell = 1.0          # m, cell size of Inputs/Bethel_DEM.asc
# optional crop box (xmin, ymin, xmax, ymax, EPSG:26913); None = whole source DEM
cropBounds = None

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dem-only", action="store_true", help="do not rewrite Inputs/REL")
a = ap.parse_args()

if not a.dem_only:
    # release area: geojson -> shapefile in the DEM CRS, 2D geometry, with the name/thickness attributes AvaFrame expects
    rel = gpd.read_file(base / "data/release_area/Bethel_W_Start_Zone.geojson")
    rel = rel.to_crs(b.CRS)
    rel["geometry"] = rel.geometry.force_2d()
    rel = rel[["geometry"]].assign(name="Bethel_W_Start_Zone", thickness=relThickness)
    (base / "Inputs/REL").mkdir(parents=True, exist_ok=True)
    rel.to_file(base / "Inputs/REL/Bethel_W_Start_Zone.shp")

with rasterio.open(b.SOURCE_DEM) as src:
    if src.crs.to_epsg() != b.CRS:
        raise SystemExit(f"{b.SOURCE_DEM.name} is {src.crs}, expected EPSG:{b.CRS}")
    l, bt, r, t = src.bounds
    if cropBounds:
        l, bt, r, t = max(l, cropBounds[0]), max(bt, cropBounds[1]), min(r, cropBounds[2]), min(t, cropBounds[3])
    # target grid: whole-metre origin inside the source bounds
    x0, y1 = math.ceil(l / cell) * cell, math.floor(t / cell) * cell
    w, h = int((math.floor(r / cell) * cell - x0) / cell), int((y1 - math.ceil(bt / cell) * cell) / cell)
    transform = from_origin(x0, y1, cell, cell)
    arr = np.full((h, w), -9999, "float32")
    # block average of the source cells (GDAL skips source nodata, so edge cells average the surveyed part)
    reproject(rasterio.band(src, 1), arr, src_transform=src.transform, src_crs=src.crs, src_nodata=src.nodata,
              dst_transform=transform, dst_crs=src.crs, dst_nodata=-9999, resampling=Resampling.average)
# write an ESRI ASCII grid, which is what AvaFrame reads as DEM
with rasterio.open(base / "Inputs/Bethel_DEM.asc", "w", driver="AAIGrid", height=h, width=w, count=1, dtype="float32",
                   crs=f"EPSG:{b.CRS}", nodata=-9999, transform=transform) as dst:
    dst.write(arr, 1)
print("done", (w, h), "cellsize", cell, "nodata cells:", int((arr == -9999).sum()), f"({(arr == -9999).mean():.0%})")
