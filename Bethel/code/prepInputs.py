"""Prepare AvaFrame inputs for Bethel: crop DEM (1 m, no resampling) to .asc, reproject release polygon to shapefile."""
import pathlib
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.windows import from_bounds

import bethel as b
base = b.BASE
relThickness = 1.0  # m, PLACEHOLDER release thickness - set to your design value

rel = gpd.read_file(base / "data/release_area/Bethel_W_Start_Zone.geojson")
rel = rel.to_crs(b.CRS)
rel["geometry"] = rel.geometry.force_2d()
rel = rel[["geometry"]].assign(name="Bethel_W_Start_Zone", thickness=relThickness)
(base / "Inputs/REL").mkdir(parents=True, exist_ok=True)
rel.to_file(base / "Inputs/REL/Bethel_W_Start_Zone.shp")

# Crop only (no resampling): keep the source 1 m grid, snap crop box to whole source cells
cropBounds = (424300, 4394200, 426400, 4396400)  # xmin, ymin, xmax, ymax in EPSG:26913 (flow path + ~500 m margin)
with rasterio.open(next((base / "data/DEM").glob("*.tif"))) as src:
    win = from_bounds(*cropBounds, src.transform).round_offsets().round_lengths()
    arr = src.read(1, window=win)
    transform = src.window_transform(win)
    nodata = src.nodata
arr = np.where(arr == nodata, -9999, arr).astype("float32")
with rasterio.open(base / "Inputs/Bethel_DEM.asc", "w", driver="AAIGrid", height=arr.shape[0], width=arr.shape[1],
                   count=1, dtype="float32", crs=f"EPSG:{b.CRS}", nodata=-9999, transform=transform) as dst:
    dst.write(arr, 1)
print("done", arr.shape[::-1], "cellsize", transform.a, "nodata cells:", int((arr == -9999).sum()))
