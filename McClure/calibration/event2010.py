"""Friction test against the 22 Feb 2010 avalanche.

Usage (from the repo root):
    uv run python McClure/calibration/event2010.py --site McClure            # run what is missing, table, figure
    uv run python McClure/calibration/event2010.py --site McClure --no-run   # table and figure only

Release: data/release_area/LoneTree_2010_release.kml (drawn by overlaying the incident report's Figure 3 in Google
Earth; the site's release in site.toml) at 1.0, 1.4 and 1.8 m. Observed extent: Feb2010_extent in
LoneTree_2010_fig3_resection.geojson (photo resection of Figure 3); its track and runout part (outside the resection's
release part) is used for the overlap scores. Friction: samosATSmall with mu varied, rho 200. Runs go to
runs/event2010/, the figure to Outputs/event2010_calibration.png and the table to Outputs/event2010_calibration.json.
Observed (report): crown 9,400 ft (2,865 m), runout 7,700 ft (2,347 m), 3-6 ft of debris over 300 ft (91 m) of road.
"""
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "siteTools"))
import avasite as b  # noqa: E402
import matplotlib  # noqa: E402
import numpy as np  # noqa: E402
import rasterio.transform as rt  # noqa: E402
from rasterio.features import rasterize  # noqa: E402
from shapely.geometry import Polygon  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

THICK = [1.0, 1.4, 1.8]
MU = [0.18, 0.22, 0.26, 0.30, 0.34, 0.38, 0.42]
RHO = 200
THR = 0.1                                   # m, footprint threshold
D = b.RUNS / "event2010"
OBS = b.BASE / "data/release_area/LoneTree_2010_fig3_resection.geojson"

g = b.readVector(OBS).to_crs(b.CRS)
g["geometry"] = g.geometry.force_2d()
ext = g[g["name"] == "Feb2010_extent"].geometry.iloc[0]
track = ext.difference(g[g["name"] == "Feb2010_release"].geometry.iloc[0])     # observed track + runout
rel = b.readVector(b.BASE / b.CFG["release"]).to_crs(b.CRS)
rel["geometry"] = rel.geometry.force_2d()
runDir = lambda th, mu: D / f"th{th:g}_mu{mu:g}"

if "--no-run" not in sys.argv:
    (D / "REL").mkdir(parents=True, exist_ok=True)
    rel[["geometry"]].assign(name="Feb2010", thickness=1.0).to_file(D / "REL/Feb2010.shp")
    for th in THICK:
        for mu in MU:
            if not (runDir(th, mu) / "Outputs/com1DFA/peakFiles").is_dir():
                b.runCom1DFA(runDir(th, mu), dict(rho=RHO, frictModel="samosATSmall", musamosatsmall=mu,
                                                  relThFromFile=False, relTh=th), rel=D / "REL")
                print("ran", th, mu, flush=True)

hill, far = b.roadLines()
strip = Polygon(list(hill.coords) + list(far.coords)[::-1]).buffer(0)
roadWidth = strip.area / hill.length
dem, TD, _ = b.readRaster(b.DEM)
rows = []
for th in THICK:
    for mu in MU:
        pft, T, _ = b.peakRaster(runDir(th, mu), "pft")
        fp = np.nan_to_num(pft) > THR
        rr, cc = np.nonzero(fp)
        r2, c2 = rt.rowcol(TD, *(T * (cc + .5, rr + .5)))
        z = dem[np.clip(r2, 0, dem.shape[0] - 1), np.clip(c2, 0, dem.shape[1] - 1)]
        obs = rasterize([track], out_shape=fp.shape, transform=T).astype(bool)
        road = rasterize([strip], out_shape=fp.shape, transform=T).astype(bool)
        mass = b.releaseMass(runDir(th, mu))
        rows.append(dict(th=th, mu=mu, massT=round(mass / 1000), d=b.dSize(mass), runoutZ=round(float(np.nanmin(z))),
                         roadHitM=round(float((fp & road).sum() * abs(T.a) ** 2 / roadWidth)),
                         obsCovered=round(float((fp & obs).sum() / obs.sum()), 2),
                         simInObs=round(float((fp & obs).sum() / fp.sum()), 2)))
        print(rows[-1])
b.OUTPUTS.mkdir(exist_ok=True)
(b.OUTPUTS / "event2010_calibration.json").write_text(json.dumps(rows, indent=1))

# figure: three runs on the observed outline
t = b.loadTerrain()
sz = b.readVector(b.BASE / "data/release_area/LoneTree_start_zone.kml").to_crs(b.CRS)
show = [(1.0, 0.26), (1.0, 0.30), (1.8, 0.34)]
x0, y0, x1, y1 = ext.bounds
fig, ax = plt.subplots(1, 3, figsize=(15, 3.7), layout="constrained")
for a, (th, mu) in zip(ax, show):
    a.imshow(t["hill"], extent=t["ext"], cmap="gray")
    pft, T, B = b.peakRaster(runDir(th, mu), "pft")
    im = a.imshow(np.where(pft > THR, pft, np.nan), extent=b.extent(B), cmap="viridis", vmax=4)
    g[g["name"] == "Feb2010_extent"].boundary.plot(ax=a, color="red", lw=1.5)
    rel.boundary.plot(ax=a, color="orange", lw=1.5)
    sz.boundary.plot(ax=a, color="white", lw=1, ls="--")
    for line in (hill, far):
        a.plot(*line.xy, "k", lw=1)
    r = next(r for r in rows if r["th"] == th and r["mu"] == mu)
    a.set_title(f"{th:g} m, μ {mu:g} ({r['massT'] / 1000:.1f} kt): runout {r['runoutZ']:,} m")
    a.set_xlim(x0 - 150, x1 + 150)
    a.set_ylim(y0 - 120, y1 + 60)
    a.set_xticks([])
    a.set_yticks([])
fig.colorbar(im, ax=ax, shrink=.8, label="peak flow thickness (m)")
fig.suptitle("22 Feb 2010: release (orange, drawn from report Fig. 3 in Google Earth), observed extent (red, photo resection "
             "of Fig. 3; runout 2,347 m), whole start zone (white dashed), road edges (black)", fontsize=10)
fig.savefig(b.OUTPUTS / "event2010_calibration.png", dpi=110)
print("saved", b.OUTPUTS / "event2010_calibration.png")
