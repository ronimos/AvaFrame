"""Prepare AvaFrame inputs for a site: source DEM -> 1 m .asc, release polygon -> shapefile, baseline config.

Usage (from the repo root):
    uv run python siteTools/prepInputs.py --site McClure              # DEM and release polygon
    uv run python siteTools/prepInputs.py --site McClure --dem-only   # DEM only; keeps Inputs/REL (e.g. an adjusted polygon)

The source DEM (site.toml `sourceDem`) is block-averaged to a `demCell` grid (default 1 m) aligned to whole metres; AvaFrame
remeshes that to its computational mesh (5 m) anyway. Cells outside the surveyed area stay nodata, so the
simulation domain is the surveyed area. Optional `cropBounds = [xmin, ymin, xmax, ymax]` (DEM CRS) limits the grid.
If the site has no Inputs/CFGs/com1DFACfg.ini yet, one is written with `rho` and `friction` from site.toml.
With a [vegetationFilter] table in site.toml, tree/shrub crowns are removed from the 1 m grid (groundFilter.py), with a
record in Inputs/ground_filter.json. With [shrubs], the shrub area (whole surveyed domain minus the road strip) is written to
Inputs/RES for the resistance runs.
"""
import json
import argparse
import math

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.warp import reproject

import avasite as b

cell = float(b.CFG.get("demCell", 1.0))
cropBounds = b.CFG.get("cropBounds")

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dem-only", action="store_true", help="do not rewrite Inputs/REL")
a = ap.parse_args()

# baseline com1DFA settings for this site (only if missing, so hand edits are kept)
cfgFile = b.INPUTS / "CFGs/com1DFACfg.ini"
if not cfgFile.exists():
    cfgFile.parent.mkdir(parents=True, exist_ok=True)
    cfgFile.write_text(f"[GENERAL]\nrho = {b.RHO:g}\nfrictModel = {b.FRICTION}\n")
    print("wrote", cfgFile)

if not a.dem_only:
    # release area -> shapefile in the DEM CRS, 2D geometry, with the name/thickness attributes AvaFrame expects
    # (thickness here is a placeholder; every run sets relTh itself)
    rel = b.readVector(b.BASE / b.CFG["release"]).to_crs(b.CRS)
    rel["geometry"] = rel.geometry.force_2d()
    rel = rel[["geometry"]].assign(name=b.RELEASE_NAME, thickness=1.0)
    b.RELEASE.parent.mkdir(parents=True, exist_ok=True)
    rel.to_file(b.RELEASE)
    print("release:", b.RELEASE, f"{rel.area.sum():,.0f} m2 planimetric")

for sid, sc in b.SCENARIOS.items():
    # named release scenarios (e.g. a D2 slab) -> Inputs/scenarios/<id>/REL/<id>.shp
    r = b.readVector(b.BASE / sc["release"]).to_crs(b.CRS)
    r["geometry"] = r.geometry.force_2d()
    r = r[["geometry"]].assign(name=sid, thickness=float(sc["thickness"]))
    b.scenarioRel(sid).mkdir(parents=True, exist_ok=True)
    r.to_file(b.scenarioRel(sid) / f"{sid}.shp")
    print(f"scenario {sid}:", b.scenarioRel(sid), f"{r.area.sum():,.0f} m2 planimetric")

with rasterio.open(b.SOURCE_DEM) as src:
    # compare the horizontal CRS (a compound CRS such as UTM 13N + NAVD88 height has no single EPSG code)
    from pyproj import CRS
    c = CRS(src.crs.to_wkt())
    if (c.sub_crs_list[0] if c.is_compound else c).to_epsg() != b.CRS:
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
vf = b.CFG.get("vegetationFilter")
if vf is not None:
    from groundFilter import groundFilter
    z = np.where(arr == -9999, np.nan, arr)
    opts = {k: vf[k] for k in ("windows", "dh0", "slope") if k in vf}
    g, veg = groundFilter(z, cell, **opts)
    removed = (z - g)[veg]
    rec = dict(source=b.SOURCE_DEM.name, cellM=cell, vegetationFraction=round(float(veg.sum() / (~np.isnan(z)).sum()), 3),
               removedMedianM=round(float(np.median(removed)), 2), removedP95M=round(float(np.percentile(removed, 95)), 2),
               removedMaxM=round(float(removed.max()), 1), **{k: v for k, v in opts.items()})
    (b.INPUTS / "ground_filter.json").write_text(json.dumps(rec, indent=1))
    print("vegetation filter:", rec)
    arr = np.where(np.isnan(g), -9999, g).astype("float32")
# write an ESRI ASCII grid, which is what AvaFrame reads as DEM
with rasterio.open(b.DEM, "w", driver="AAIGrid", height=h, width=w, count=1, dtype="float32",
                   crs=f"EPSG:{b.CRS}", nodata=-9999, transform=transform) as dst:
    dst.write(arr, 1)
print("DEM:", b.DEM, (w, h), "cellsize", cell, "nodata cells:", int((arr == -9999).sum()), f"({(arr == -9999).mean():.0%})")

if b.SHRUBS:
    # shrub resistance area as a raster (com1DFA: values > 0 = resistance): every surveyed cell except the road strip.
    # It is written on AvaFrame's computational grid (the DEM remeshed to meshCellSize by AvaFrame's own function), so
    # com1DFA uses it as is. A raster avoids AvaFrame's polygon rules (no holes, every vertex on a valid cell), and
    # matching the grid avoids its own remeshing of the raster, which trips over the .prj file in remeshedRasters.
    import configparser
    from avaframe.in3Utils import geoTrans
    from rasterio.features import rasterize
    from shapely.geometry import Polygon
    c = configparser.ConfigParser(interpolation=None)
    c.read([b.REPO / "avaframe/com1DFA/com1DFACfg.ini", cfgFile])
    mesh = c.getfloat("GENERAL", "meshCellSize")
    hd = geoTrans.remeshDataRio(b.DEM, mesh, larger=False)["header"]  # same call as AvaFrame's remeshRaster
    rT = from_origin(hd["xllcenter"] - mesh / 2, hd["yllcenter"] - mesh / 2 + hd["nrows"] * mesh, mesh, mesh)
    # shrubs where the 1 m DEM is valid (average of the valid flag > 0.5) ...
    frac = np.zeros((hd["nrows"], hd["ncols"]), "float32")
    reproject((arr != -9999).astype("float32"), frac, src_transform=transform, src_crs=f"EPSG:{b.CRS}", dst_transform=rT,
              dst_crs=f"EPSG:{b.CRS}", resampling=Resampling.average)
    mask = frac > .5
    # ... except on the road
    if b.hasRoad():
        hillL, farL = b.roadLines()
        strip = Polygon(list(hillL.coords) + list(farL.coords)[::-1]).buffer(0)
        mask &= ~rasterize([strip], out_shape=mask.shape, transform=rT).astype(bool)
    # optionally not in the release area ([shrubs] excludeRelease): otherwise com1DFA's detrainment (flow slower than
    # forestVMin) drains thin releases before they start moving
    if b.SHRUBS.get("excludeRelease", False):
        mask &= ~rasterize(list(b.readVector(b.RELEASE).geometry), out_shape=mask.shape, transform=rT, all_touched=True).astype(bool)
    for old in b.RES.parent.glob("*"):     # only one resistance file may be in RES/
        old.unlink()
    b.RES.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(b.RES, "w", driver="AAIGrid", height=mask.shape[0], width=mask.shape[1], count=1, dtype="float32",
                       crs=f"EPSG:{b.CRS}", transform=rT) as dst:
        dst.write(mask.astype("float32"), 1)
    for prj in b.RES.parent.glob("*.prj"):  # AvaFrame globs RES/ for rasters; keep only the .asc
        prj.unlink()
    print("shrub resistance area:", b.RES, f"{mask.sum() * mesh**2 / 1e4:.1f} ha on the {mesh:g} m grid")
