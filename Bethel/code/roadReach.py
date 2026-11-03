"""Does the avalanche reach the road? Scenario matrix: release thickness x friction class.

Road = strip between the hill-side and far-side edge lines in Bethel/data/Road. For every run:
  reaches road  - any cell with peak flow thickness (pft) > thr touches the hill-side edge line
  crosses road  - same for the far-side edge line
plus max pft / pfv / ppr on the road strip, road area covered, and the gap to the road when it is not reached.
Usage (from /home/ron/AvaFrame):  uv run python Bethel/code/roadReach.py [--render-only] [--thr 0.1]
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rasterio.features import rasterize
from shapely import distance, points
from shapely.geometry import Polygon

import bethel as b

# Scenario matrix: release thicknesses [m] x friction classes. Fixed classes (not samosATAuto) so the matrix
# shows how sensitive the road result is to friction.
THICK = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
FRICT = ["samosATMedium", "samosAT"]
# snow density [kg/m3], also used to turn the simulated release mass into a volume
RHO = 300


def runDir(fr, t):
    return b.RUNS / f"road_{fr}" / f"th_{t:g}"


def main():
    b.OUTPUTS.mkdir(exist_ok=True)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--render-only", action="store_true")
    ap.add_argument("--thr", type=float, default=0.1, help="peak flow thickness [m] counted as 'flow'")
    a = ap.parse_args()

    # 1) simulations: one avalanche dir per (friction, thickness) under runs/road_<friction>/
    if not a.render_only:
        for fr in FRICT:
            for t in THICK:
                b.runCom1DFA(runDir(fr, t), dict(rho=RHO, frictModel=fr, relThFromFile=False, relTh=f"{t:g}"))
                print(f"{fr}, relTh {t:g} m done")

    # 2) road geometry: the two edge lines (hill side / far side) and the strip between them, in the DEM's CRS.
        hillL, farL = b.roadLines()
    strip = Polygon(list(hillL.coords) + list(farL.coords)[::-1]).buffer(0)

    # 3) rasterise the lines/strip onto the simulation grid (all_touched so thin lines still hit cells)
    first, T, B = b.peakRaster(runDir(FRICT[0], THICK[0]), "pft")
    shape = first.shape
    cell = abs(T.a)
    rast = lambda g, touch: rasterize([g], out_shape=shape, transform=T, all_touched=touch).astype(bool)
    mHill, mFar, mStrip = rast(hillL, True), rast(farL, True), rast(strip, True)
    rows, res = [], {}
    # cell-centre coordinates, used for the distance to the road when the flow stops short
    rr, cc = np.mgrid[0:shape[0], 0:shape[1]]
    X, Y = T * (cc + .5, rr + .5)
    for fr in FRICT:
        for t in THICK:
            d = runDir(fr, t)
            pft, pfv, ppr = (b.peakRaster(d, f)[0] for f in ("pft", "pfv", "ppr"))
            # kept for the map matrix below
            res[fr, t] = pft
            # 'flow' = cells whose peak flow thickness exceeds --thr
            flow = pft > a.thr
            reaches, crosses = bool((flow & mHill).any()), bool((flow & mFar).any())
            onRoad = flow & mStrip
            # gap to the road: 0 if the flow touches the hill-side edge, otherwise the distance [m] from the nearest flow cell
            if reaches:
                gap = 0.0
            else:  # distance from the nearest flow cell to the hill-side edge
                gap = float(distance(points(X[flow], Y[flow]), hillL).min())
            rows.append(dict(friction=fr, relTh_m=t, reachesRoad=reaches, crossesRoad=crosses, gapToRoad_m=round(gap),
                             roadAreaHit_m2=int(onRoad.sum() * cell**2), maxPFT_onRoad_m=round(float(np.nanmax(np.where(mStrip, pft, np.nan))), 2),
                             maxPFV_onRoad_ms=round(float(np.nanmax(np.where(mStrip, pfv, np.nan))), 1),
                             maxPPR_onRoad_kPa=round(float(np.nanmax(np.where(mStrip, ppr, np.nan)))),
                             volume_m3=round(b.releaseMass(d) / RHO)))
    # 4) results table (csv + printed)
    df = pd.DataFrame(rows)
    df.to_csv(b.OUTPUTS / "road_reach.csv", index=False)
    print(f"flow = peak thickness > {a.thr} m; road = strip between the two edge lines\n")
    print(df.to_string(index=False))

    # map matrix zoomed on the road
    # 5) map matrix, zoomed on the road; xl/yl are Bethel-specific map limits (EPSG:26913 easting/northing)
    terrain = b.loadTerrain()
    rel, hill, dExt, ext = terrain["rel"], terrain["hill"], terrain["ext"], b.extent(B)
    xl, yl = (425250, 426050), (4394500, 4395250)
    vmax = max(np.nanmax(v) for v in res.values())
    fig, axs = plt.subplots(len(FRICT), len(THICK), figsize=(4.6 * len(THICK) + 1, 4.6 * len(FRICT) + .8), constrained_layout=True)
    for i, fr in enumerate(FRICT):
        for j, t in enumerate(THICK):
            ax = axs[i, j]; ax.imshow(hill, cmap="gray", extent=dExt, alpha=.85)
            im = ax.imshow(np.where(res[fr, t] > a.thr, res[fr, t], np.nan), cmap="Blues", vmin=0, vmax=vmax, extent=ext)
            for g, c in ((hillL, "#d95f02"), (farL, "#7b3294")):
                ax.plot(*np.array(g.coords)[:, :2].T, color=c, lw=2)
            row = df[(df.friction == fr) & (df.relTh_m == t)].iloc[0]
            status = "CROSSES road" if row.crossesRoad else ("reaches road" if row.reachesRoad else f"stops {row.gapToRoad_m:.0f} m short")
            ax.set_title(f"{fr}, {t:g} m: {status}", fontsize=9, loc="left")
            ax.set_xlim(xl); ax.set_ylim(yl); ax.set_aspect("equal"); ax.tick_params(labelsize=6)
            ax.ticklabel_format(style="plain", useOffset=False)
    fig.colorbar(im, ax=axs, shrink=.6, label="peak flow thickness (m)")
    fig.suptitle(f"Peak flow thickness near the road (orange = hill-side edge, purple = far-side edge; flow = pft > {a.thr:g} m)", fontsize=11)
    fig.savefig(b.OUTPUTS / "road_reach.png", dpi=90)
    print("\nsaved", b.OUTPUTS / "road_reach.png", "and road_reach.csv")


if __name__ == "__main__":
    main()
