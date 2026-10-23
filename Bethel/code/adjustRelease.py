"""Adjust the Bethel release polygon downslope so it ends where slope drops below 28 deg.

Steepest-descent (D8) paths are traced from every cell of the original polygon on a 5 m block-mean DEM.
A path is followed while slope >= minSlope and stops after nBelow consecutive cells below it.
New zone = traced steep cells; original-polygon cells that are not steep and not reached from above are dropped.
"""
import json
import pathlib
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize, shapes
from scipy import ndimage as ndi
from shapely.geometry import shape

import bethel as b
base = b.BASE
minSlope, nBelow, f = 28.0, 3, 5

with rasterio.open(base / "Inputs/Bethel_DEM.asc") as r:
    z = r.read(1).astype(float); T = r.transform
n0, n1 = z.shape[0] // f * f, z.shape[1] // f * f
z5 = z[:n0, :n1].reshape(n0 // f, f, n1 // f, f).mean(axis=(1, 3))
T5 = T * T.scale(f, f)
gy, gx = np.gradient(z5, f)
slope = ndi.uniform_filter(np.degrees(np.arctan(np.hypot(gx, gy))), 3)
steep = slope >= minSlope

orig = gpd.read_file(base / "Inputs/REL/Bethel_W_Start_Zone.shp") if not (base / "data/Bethel_W_orig.shp").exists() \
    else gpd.read_file(base / "data/Bethel_W_orig.shp")
origMask = rasterize(orig.geometry, out_shape=z5.shape, transform=T5).astype(bool)

nb = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
dist = np.array([np.hypot(*d) for d in nb])
zone = np.zeros_like(steep)
for r0, c0 in zip(*np.nonzero(origMask & steep)):
    r_, c_, below, path = r0, c0, 0, [(r0, c0)]
    for _ in range(400):
        best, bd = None, 0.0
        for (dr, dc), L in zip(nb, dist):
            rr, cc = r_ + dr, c_ + dc
            if 0 <= rr < z5.shape[0] and 0 <= cc < z5.shape[1]:
                d = (z5[r_, c_] - z5[rr, cc]) / L
                if d > bd: best, bd = (rr, cc), d
        if best is None: break
        r_, c_ = best
        below = 0 if steep[r_, c_] else below + 1
        if below >= nBelow: path = path[:len(path) - (below - 1)]; break
        path.append((r_, c_))
    for p in path: zone[p] = True

zone = ndi.binary_closing(zone, iterations=2)
zone = ndi.binary_fill_holes(ndi.binary_opening(zone, iterations=1))
lab, n = ndi.label(zone)
if n > 1:  # keep components that overlap the original polygon
    keep = np.unique(lab[origMask & zone]); keep = keep[keep > 0]
    zone = np.isin(lab, keep)
polys = [shape(g) for g, v in shapes(zone.astype("uint8"), mask=zone, transform=T5)]
geom = gpd.GeoSeries(polys, crs=b.CRS).union_all().simplify(5).buffer(0)
new = gpd.GeoDataFrame({"name": ["Bethel_W_Start_Zone"], "thickness": [1.0]}, geometry=[geom], crs=b.CRS)

# keep original, write adjusted
(base / "data").mkdir(exist_ok=True)
if not (base / "data/Bethel_W_orig.shp").exists():
    orig.to_file(base / "data/Bethel_W_orig.shp")
new.to_file(base / "Inputs/REL/Bethel_W_Start_Zone.shp")
a0, a1 = orig.geometry.area.sum(), new.geometry.area.sum()
print(f"area original {a0/1e4:.2f} ha -> adjusted {a1/1e4:.2f} ha  ({(a1-a0)/a0*100:+.0f}%)")
inter = orig.geometry.iloc[0].intersection(geom).area
print(f"kept from original: {inter/1e4:.2f} ha; added: {(a1-inter)/1e4:.2f} ha; removed: {(a0-inter)/1e4:.2f} ha")
print("geometry type:", geom.geom_type, " bounds:", [round(b) for b in geom.bounds])

# parameters and areas, picked up by runThicknessSeries.py for the methodology text
(base / "Inputs/release_adjustment.json").write_text(json.dumps(dict(
    minSlopeDeg=minSlope, nBelowCells=nBelow, blockCellM=f, slopeSmoothingCells=3, simplifyM=5,
    areaOriginalHa=a0 / 1e4, areaAdjustedHa=a1 / 1e4, addedHa=(a1 - inter) / 1e4, removedHa=(a0 - inter) / 1e4), indent=1))
