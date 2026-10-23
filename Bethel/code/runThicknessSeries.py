"""Run com1DFA for Bethel at several release thicknesses and plot peak flow velocity / thickness.

Usage (from /home/ron/AvaFrame):
    uv run python Bethel/code/runThicknessSeries.py                 # run all default thicknesses, then plot
    uv run python Bethel/code/runThicknessSeries.py --plot-only     # reuse existing runs
    uv run python Bethel/code/runThicknessSeries.py --clean         # delete the series' run directories first
    uv run python Bethel/code/runThicknessSeries.py --thickness 1 2 --rho 250 --friction samosATMedium

Writes a static PNG and an interactive HTML (plotly, hover shows velocity/thickness/pressure per cell) to Outputs/.
Each thickness is its own AvaFrame avalanche dir under runs/ (settings = Inputs/CFGs/com1DFACfg.ini plus the
overrides of that run), run in its own process.
The release shapefile (Inputs/REL) is used as is, run adjustRelease.py first if it needs adjusting.
"""
import argparse
import configparser
import pathlib
import shutil
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio

import bethel as b

footprintThr = 0.1  # m, pft threshold used to define the footprint / zoom extent


def seriesDir(rho, friction):
    return b.RUNS / f"thickness_{friction}_rho{rho:g}"


def runSeries(thicknesses, rho, friction):
    """One AvaFrame run (own process, own avalanche dir) per release thickness; returns {thickness: avaDir}."""
    runs = {}
    for t in thicknesses:
        runs[t] = seriesDir(rho, friction) / f"th_{t:g}"
        b.runCom1DFA(runs[t], dict(rho=f"{rho:g}", frictModel=friction, relThFromFile=False, relTh=f"{t:g}"))
        print(f"relTh {t:g} m done")
    return runs


def findRuns(thicknesses, rho, friction):
    """Existing run dirs {thickness: avaDir} of a series."""
    runs = {t: seriesDir(rho, friction) / f"th_{t:g}" for t in thicknesses}
    missing = [t for t, d in runs.items() if not (d / "Outputs/com1DFA/peakFiles").is_dir()]
    if missing:
        sys.exit(f"no results found for relTh {sorted(missing)}")
    return runs


def loadSeries(runs):
    """Peak rasters of all runs plus hillshade, release polygon and a zoom window around the flow paths."""
    th = sorted(runs)
    pft, pfv, ppr = ({t: b.peakRaster(runs[t], f)[0] for t in th} for f in ("pft", "pfv", "ppr"))
    _, T, B = b.peakRaster(runs[th[0]], "pft")
    terrain = b.loadTerrain()
    rel = terrain["rel"]
    xl, yl = b.zoomWindow([pft[t] > footprintThr for t in th], T, rel, pad=120)
    return dict(th=th, pft=pft, pfv=pfv, ppr=ppr, T=T, meshExt=b.extent(B), rel=rel, hill=terrain["hill"],
                demExt=terrain["ext"], xl=xl, yl=yl)


# HTML skeleton of the report; {title}, {div} (plotly figure), {anim} (video block) and {method} (methodology) are filled in plotInteractive()
PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
 body {{ font-family: system-ui, sans-serif; margin: 0 auto; max-width: 1300px; padding: 0 16px 48px; color: #222; line-height: 1.5; }}
 h1 {{ font-size: 1.4rem; margin: 18px 0 2px; }} h2 {{ font-size: 1.25rem; margin-top: 28px; border-bottom: 1px solid #ccc; }}
 h3 {{ font-size: 1rem; margin: 18px 0 2px; }}
 .sub {{ color: #555; margin: 0 0 8px; }} .method {{ max-width: 900px; }}
 table {{ border-collapse: collapse; margin: 8px 0; }} td, th {{ border: 1px solid #ccc; padding: 3px 10px; text-align: left; }}
 code {{ background: #f2f2f2; padding: 0 3px; }}
</style></head><body>
<h1>{title}</h1>
<p class="sub">By release thickness. Hover for values; zoom and pan are linked across panels. Orange outline = release area.</p>
{div}
{anim}
<div class="method">{method}</div>
</body></html>"""


def runInfo(runs, rho):
    """Per thickness: friction class AvaFrame used, release mass [kg] and volume [m3] (initial mass / density)."""
    info = {}
    for t, tag in runs.items():
        mass = b.releaseMass(tag)
        info[t] = dict(cls=b.runConfig(tag)["GENERAL"]["frictModel"], mass=mass, vol=mass / rho)
    return info


def defaultCfg():
    """AvaFrame's default friction parameters and volume classes."""
    c = configparser.ConfigParser(interpolation=None)
    c.read(b.REPO / "avaframe/com1DFA/com1DFACfg.ini")
    g = c["GENERAL"]
    return dict(volSmall=float(g["volClassSmall"]), volMedium=float(g["volClassMedium"]),
                mu={"samosATSmall": float(g["musamosatsmall"]), "samosATMedium": float(g["musamosatmedium"]),
                    "samosAT": float(g["musamosat"])})


def animationHtml(friction):
    """Embedded MP4 of the side-by-side animation, if it exists."""
    import base64
    d = b.OUTPUTS / f"animation_{friction}"
    srcs = "".join(f'<source src="data:{mime};base64,{base64.b64encode(f.read_bytes()).decode()}" type="{mime}">'
                   for f, mime in ((d / "flow.webm", "video/webm"), (d / "flow.mp4", "video/mp4")) if f.exists())
    if not srcs:
        return ""
    return f"""<h2>Animation</h2>
<p class="sub">Flow thickness (top) and velocity (bottom) for each release thickness on a shared clock; a run that has stopped holds its last frame.
 Saved every 2&nbsp;s of simulated time. Colour scales are fixed across columns and frames.
 Also saved as <code>flow.gif</code> and individual frames in <code>Outputs/animation_{friction}/</code>.</p>
<video controls loop muted autoplay playsinline style="width:100%;max-width:1300px">{srcs}Your browser cannot play this video; see flow.gif in the output folder.</video>"""


def frictionParagraph(info, dc, th, auto, friction):
    """Explain the friction choice from the actual classes/volumes of this run set."""
    if not auto:
        return (f"<p>All runs use <code>{friction}</code> (mu = {dc['mu'].get(friction, float('nan')):g}). It is the same for every depth, so "
                "differences between rows come from release thickness (and so volume) only.</p>")
    vs, vm = dc["volSmall"], dc["volMedium"]
    mu = dc["mu"]
    txt = f"""<p><b>How.</b> Friction was set with <code>samosATAuto</code>. AvaFrame computes each simulation's release volume and picks one of
 three calibrations of the samosAT friction model by volume: <code>samosATSmall</code> (&mu; = {mu['samosATSmall']:g}) below {vs:,.0f}&nbsp;m&sup3;,
 <code>samosATMedium</code> (&mu; = {mu['samosATMedium']:g}) from {vs:,.0f} to {vm:,.0f}&nbsp;m&sup3;, and standard <code>samosAT</code>
 (&mu; = {mu['samosAT']:g}) above {vm:,.0f}&nbsp;m&sup3;. The class is chosen per depth, not by hand; the table shows what each depth got
 (release volume = initial simulated mass / density).</p>
<p><b>Why.</b> The calibrations assume that larger avalanches are less resisted (lower &mu;). Release thickness changes the volume
 from {min(i['vol'] for i in info.values()):,.0f} to {max(i['vol'] for i in info.values()):,.0f}&nbsp;m&sup3;, so one fixed class would be mis-calibrated for some of
 the depths. In particular <code>samosATSmall</code> is meant for volumes below {vs:,.0f}&nbsp;m&sup3;
 and every scenario here exceeds that, so it was not used. An earlier version of this analysis used <code>samosATSmall</code> for all depths and was
 replaced for that reason.</p>
<p><b>Caveat.</b> The class boundaries are step changes. Where a depth sits just on one side of a boundary, the friction jumps between neighbouring
 depths, so part of the difference between those rows is a friction change, not only a thickness change."""
    near = [t for t in th if min(abs(info[t]["vol"] - vs), abs(info[t]["vol"] - vm)) < 5000]
    if near:
        txt += " Closest to a boundary: " + ", ".join(
            f"{t:g}&nbsp;m ({info[t]['vol']:,.0f}&nbsp;m&sup3; vs the {vm if abs(info[t]['vol'] - vm) < abs(info[t]['vol'] - vs) else vs:,.0f}&nbsp;m&sup3; limit)" for t in near) + "."
    txt += "</p>"
    chk = b.OUTPUTS / "road_reach.csv"
    if chk.exists():
        import pandas as pd
        d = pd.read_csv(chk)
        pm = d[d.friction == "samosATMedium"].set_index("relTh_m"); ps = d[d.friction == "samosAT"].set_index("relTh_m")
        same = all((pm.loc[x, "reachesRoad"] == ps.loc[x, "reachesRoad"]) and (pm.loc[x, "crossesRoad"] == ps.loc[x, "crossesRoad"]) for x in pm.index)
        both = [x for x in pm.index if pm.loc[x, "maxPPR_onRoad_kPa"] > 0 and ps.loc[x, "maxPPR_onRoad_kPa"] > 0]
        rat = max(ps.loc[x, "maxPPR_onRoad_kPa"] / pm.loc[x, "maxPPR_onRoad_kPa"] for x in both) if both else float("nan")
        txt += (f"""<p><b>Sensitivity check.</b> A separate test (<code>roadReach.py</code>, thicknesses {pm.index.min():g}&ndash;{pm.index.max():g}&nbsp;m)
 ran every depth with both <code>samosATMedium</code> and <code>samosAT</code>. """
                + ("Whether the flow reached or crossed the road was the same under both classes at every depth tested; " if same else
                   "Whether the flow reached or crossed the road differed between the classes for at least one depth; ")
                + f"the peak pressure on the road differed by up to a factor of {rat:.2f}."
                + ("" if set(th) <= set(pm.index) else " The test did not include every depth shown above.") + "</p>")
    return txt


def methodologyHtml(runs, rho, friction):
    """Methodology section; numbers come from the run configs and Inputs/release_adjustment.json."""
    import datetime
    import importlib.metadata as md
    th = sorted(runs)
    mesh = float(b.runConfig(runs[th[0]])["GENERAL"]["meshCellSize"])
    try:
        ver = md.version("avaframe").rsplit(".d2", 1)[0]  # drop uv's dirty-tree date suffix
    except md.PackageNotFoundError:
        ver = "unknown"
    demName = next((b.DATA / "DEM").glob("*.tif")).name
    with rasterio.open(b.DEM) as r:
        demInfo = f"{r.width} x {r.height} cells at {abs(r.transform.a):g} m"
    adj = b.INPUTS / "release_adjustment.json"
    if adj.exists():
        import json
        a = json.loads(adj.read_text())
        release = f"""<p>The release area started from a manually drawn polygon ({a['areaOriginalHa']:.2f} ha). Its downslope
 extent was then set by slope: the zone follows terrain that is <b>&ge; {a['minSlopeDeg']:g}&deg; in the downhill direction</b>.</p>
<ol>
<li>Slope was computed from the 1&nbsp;m DEM after averaging to {a['blockCellM']}&nbsp;m blocks, using a finite-difference gradient and a
 {a['slopeSmoothingCells']}&times;{a['slopeSmoothingCells']} mean filter.</li>
<li>From every cell of the original polygon that is &ge; {a['minSlopeDeg']:g}&deg;, a steepest-descent path (D8 flow direction on the
 {a['blockCellM']}&nbsp;m DEM) was followed downhill. A path continues while slope stays &ge; {a['minSlopeDeg']:g}&deg; and stops after
 {a['nBelowCells']} consecutive cells ({a['nBelowCells'] * a['blockCellM']} m) below it; the trailing sub-threshold cells are dropped.</li>
<li>The new zone is the union of those paths. Polygon cells that are below {a['minSlopeDeg']:g}&deg; and not reached from a steeper cell above
 are removed (shrink); cells reached beyond the old boundary are added (expand). Expansion is therefore limited to the
 downslope flow corridors of the original zone, it does not spread sideways over other steep terrain.</li>
<li>The mask was cleaned (binary closing/opening, holes filled, only parts overlapping the original polygon kept), polygonised and
 simplified to {a['simplifyM']} m.</li>
</ol>
<p>Result: {a['areaAdjustedHa']:.2f} ha (planimetric), {a['addedHa']:.2f} ha added and {a['removedHa']:.2f} ha removed relative to the
 original ({(a['areaAdjustedHa'] / a['areaOriginalHa'] - 1) * 100:+.0f}%). Release thickness is uniform over the zone and is the varied parameter.</p>"""
    else:
        release = "<p>The release area is the shapefile in <code>Inputs/REL</code>; no adjustment record was found.</p>"
    info, dc = runInfo(runs, rho), defaultCfg()
    rows = "".join(f"<tr><td>{t:g} m</td><td>{info[t]['vol']:,.0f} m&sup3;</td><td><code>{info[t]['cls']}</code></td>"
                   f"<td>{dc['mu'][info[t]['cls']]:g}</td><td><code>{runs[t].relative_to(b.BASE)}</code></td></tr>" for t in th)
    auto = friction == "samosATAuto"
    frictionText = frictionParagraph(info, dc, th, auto, friction)
    return f"""
<h2>Methodology</h2>
<h3>Terrain</h3>
<p>Elevation: USGS 1&nbsp;m lidar DEM (<code>{demName}</code>, EPSG:26913 / UTM 13N), cropped to the flow path plus a margin without resampling
 ({demInfo}). AvaFrame then remeshes it internally to a {mesh:g}&nbsp;m computational mesh, so the 1&nbsp;m detail is averaged before the solver runs.
 The hillshade in the plots is from the 1&nbsp;m DEM.</p>
<h3>Release area</h3>
{release}
<h3>Simulations</h3>
<p>Dense-flow avalanche simulations with AvaFrame <code>com1DFA</code> (version {ver}). Constant for all runs: snow density
 {rho:g}&nbsp;kg/m&sup3;, {mesh:g}&nbsp;m computational mesh, release area as above, all other parameters at AvaFrame defaults.
 The only input that differs between rows is the release thickness ({", ".join(f"{t:g}" for t in th)}&nbsp;m); each value is a separate simulation.
 Friction is set per run by AvaFrame from the release volume (next section).</p>
<h3>Friction model: how and why each depth got its class</h3>
{frictionText}
<table><tr><th>Release thickness</th><th>Release volume</th><th>Friction class used</th><th>&mu;</th><th>Run directory</th></tr>{rows}</table>
<h3>Outputs shown</h3>
<p>Peak flow velocity (m/s) and peak flow thickness (m): the maximum reached at each cell during the simulation, not a snapshot. Cells are
 shown where peak thickness exceeds {footprintThr:g}&nbsp;m. Colour scales are shared per column. Hovering shows easting/northing (EPSG:26913) and
 the cell's peak velocity, thickness and pressure.</p>
<h3>Limitations</h3>
<ul>
<li>Results are from a single deterministic model setup and have not been calibrated or validated against observed events at this site.</li>
<li>Friction classes come from volume bands and change in steps; a class change between depths is a modelling choice, not a measured difference.</li>
<li>Release thickness is a scenario input, not a measurement; slope and release-area extent depend on the 5&nbsp;m averaging and the {a['minSlopeDeg'] if adj.exists() else 28:g}&deg; threshold.</li>
</ul>
<p class="sub">Generated {datetime.date.today().isoformat()} by <code>Bethel/runThicknessSeries.py</code>.</p>"""


def plotInteractive(runs, rho, friction, fname):
    """Plotly version: hover shows easting/northing and thickness, velocity, pressure of that cell for the run."""
    import base64
    import io
    from PIL import Image
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    s = loadSeries(runs)
    th, T, xl, yl = s["th"], s["T"], s["xl"], s["yl"]
    nrow, ncol = s["pft"][th[0]].shape
    xc = T.c + T.a * (np.arange(ncol) + .5)   # cell-centre easting
    yc = T.f + T.e * (np.arange(nrow) + .5)   # cell-centre northing (decreasing)
    cx = (xc >= xl[0]) & (xc <= xl[1])
    cy = (yc >= yl[0]) & (yc <= yl[1])
    # heatmaps only need the footprint's bounding box (everything else is empty), which keeps the HTML small
    foot = np.any([s["pft"][t] > footprintThr for t in th], axis=0)
    fr, fc = np.nonzero(foot)
    cx &= (np.arange(ncol) >= fc.min() - 1) & (np.arange(ncol) <= fc.max() + 1)
    cy &= (np.arange(nrow) >= fr.min() - 1) & (np.arange(nrow) <= fr.max() + 1)
    xc, yc = xc[cx], yc[cy]
    crop = lambda a: a[np.ix_(cy, cx)].astype("float32")

    # hillshade as one PNG (reused as background in every panel), clipped to the zoom window and thinned to 2 m
    dx0, dx1, dy0, dy1 = s["demExt"]
    hs = s["hill"]
    cs = (dx1 - dx0) / hs.shape[1]
    c0, c1 = int((xl[0] - dx0) / cs), int((xl[1] - dx0) / cs)
    r0, r1 = int((dy1 - yl[1]) / cs), int((dy1 - yl[0]) / cs)
    hs = hs[r0:r1:2, c0:c1:2]
    buf = io.BytesIO()
    Image.fromarray((np.clip(np.nan_to_num(hs), 0, 1) * 255).astype("uint8")).save(buf, format="PNG")
    uri = "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()
    imgX0, imgY1 = dx0 + c0 * cs, dy1 - r0 * cs
    imgW, imgH = (c1 - c0) * cs, (r1 - r0) * cs

    # two columns per thickness: peak velocity (left) and peak thickness (right); coloraxis/coloraxis2 give each column its own scale
    cols = [("Peak flow velocity", "pfv", "m/s", "coloraxis"), ("Peak flow thickness", "pft", "m", "coloraxis2")]
    # friction class used and release volume per thickness, shown in the panel titles
    info = runInfo(runs, rho)
    titles = [f"{n} - release {t:g} m ({info[t]['cls']})" for t in th for n, *_ in cols]
    fig = make_subplots(rows=len(th), cols=2, subplot_titles=titles, horizontal_spacing=0.1, vertical_spacing=0.04)
    relXY = [np.array(g.exterior.coords) for g in s["rel"].geometry]
    cmax = {}
    for _, key, *_ in cols:
        cmax[key] = max(np.nanmax(np.where(s["pft"][t] > footprintThr, s[key][t], np.nan)) for t in th)

    # one row per thickness; cells outside the footprint are masked (nan) so the hillshade shows through
    for i, t in enumerate(th, start=1):
        mask = crop(s["pft"][t]) > footprintThr
        dec = {"pfv": 1, "pft": 2, "ppr": 0}  # rounding keeps the HTML small
        z = {k: np.round(np.where(mask, crop(s[k][t]), np.nan), dec[k]).astype("float32") for k in dec}
        for j, (name, key, unit, cax) in enumerate(cols, start=1):
            # customdata: velocity / thickness / pressure of the same cell, so either panel shows all three
            cd = np.dstack([z["pfv"], z["pft"], z["ppr"]])
            fig.add_trace(go.Heatmap(
                x=xc, y=yc, z=z[key], customdata=cd, coloraxis=cax, opacity=0.9, name="",
                hovertemplate=("E %{x:.0f}  N %{y:.0f}<br>"
                               f"<b>{name}: %{{z:.2f}} {unit}</b><br>"
                               "velocity %{customdata[0]:.1f} m/s<br>thickness %{customdata[1]:.2f} m<br>"
                               "pressure %{customdata[2]:.0f} kPa<extra>release " + f"{t:g} m</extra>"),
            ), row=i, col=j)
            for xy in relXY:
                fig.add_trace(go.Scatter(x=xy[:, 0], y=xy[:, 1], mode="lines", line=dict(color="#d95f02", width=2),
                                         hoverinfo="skip", showlegend=False), row=i, col=j)
            ax = (i - 1) * 2 + j
            suffix = "" if ax == 1 else str(ax)
            # hillshade as a background image in every panel (one PNG data URI, cropped to the zoom window and thinned to 2 m to keep the HTML small)
            fig.add_layout_image(dict(source=uri, xref=f"x{suffix}", yref=f"y{suffix}", x=imgX0, y=imgY1,
                                      sizex=imgW, sizey=imgH, xanchor="left", yanchor="top", layer="below", sizing="stretch"))
            fig.update_xaxes(range=xl, matches="x", tickformat="d", tickfont_size=9, showgrid=False, row=i, col=j)
            fig.update_yaxes(range=yl, matches="y", tickformat="d", tickfont_size=9, showgrid=False, scaleanchor=f"x{suffix}",
                             constrain="domain", row=i, col=j)

    # shared colour scales per column: the same value gets the same colour in every row
    bar = dict(len=0.9, thickness=14, y=0.5)
    fig.update_layout(
        coloraxis=dict(colorscale="Blues", cmin=0, cmax=cmax["pfv"], colorbar=dict(title="m/s", x=0.455, **bar)),
        coloraxis2=dict(colorscale="Blues", cmin=0, cmax=cmax["pft"], colorbar=dict(title="m", x=1.0, **bar)),
        height=520 * len(th), width=1250, margin=dict(l=70, r=40, t=40, b=40), plot_bgcolor="rgba(0,0,0,0)",
        hoverlabel=dict(font_size=12),
    )
    fig.update_annotations(font_size=12)
    # plotly.js is embedded, so the HTML works offline as a single file
    div = fig.to_html(full_html=False, include_plotlyjs=True)
    pathlib.Path(fname).write_text(PAGE.format(
        title="Bethel avalanche: peak flow velocity and thickness", div=div, anim=animationHtml(friction),
        method=methodologyHtml(runs, rho, friction)))
    print("saved", fname)


def plot(runs, rho, friction, fname):
    s = loadSeries(runs)
    th, pft, pfv, hill, demExt, meshExt, rel, xl, yl = (s[k] for k in
                                                       ("th", "pft", "pfv", "hill", "demExt", "meshExt", "rel", "xl", "yl"))
    cols = [("Peak flow velocity", pfv, "m/s"), ("Peak flow thickness", pft, "m")]
    # shared colour maximum per column, taken over footprint cells of all thicknesses
    vmax = [max(np.nanmax(np.where(pft[t] > footprintThr, d[t], np.nan)) for t in th) for _, d, _ in cols]
    # static version of the same figure (PNG)
    fig, axs = plt.subplots(len(th), 2, figsize=(11, 4.3 * len(th)), constrained_layout=True, squeeze=False)
    for i, t in enumerate(th):
        for j, (name, d, unit) in enumerate(cols):
            ax = axs[i, j]
            ax.imshow(hill, cmap="gray", extent=demExt, alpha=.8)
            im = ax.imshow(np.where(pft[t] > footprintThr, d[t], np.nan), cmap="Blues", vmin=0, vmax=vmax[j], extent=meshExt)
            rel.boundary.plot(ax=ax, color="#d95f02", lw=1.5)
            ax.set_xlim(xl); ax.set_ylim(yl); ax.set_aspect("equal")
            ax.set_title(f"{name} - release thickness {t:g} m", fontsize=10, loc="left")
            ax.tick_params(labelsize=7); ax.ticklabel_format(style="plain", useOffset=False)
            fig.colorbar(im, ax=ax, shrink=.75, label=unit)
    fig.suptitle(f"Bethel com1DFA: rho {rho:g}, {friction} (orange = release area; colour scales shared per column)", fontsize=11)
    fig.savefig(fname, dpi=100)
    print("saved", fname)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--thickness", type=float, nargs="+", default=[0.5, 0.75, 1.0, 1.5, 2.0], help="release thicknesses [m]")
    ap.add_argument("--rho", type=float, default=300, help="snow density [kg/m3]")
    ap.add_argument("--friction", default="samosATAuto", help="com1DFA frictModel (samosATAuto picks the class per release volume)")
    ap.add_argument("--clean", action="store_true", help="delete the series' run directories before running")
    ap.add_argument("--plot-only", action="store_true", help="skip simulations, only (re)plot existing results")
    ap.add_argument("--out", default=str(b.OUTPUTS / "thickness_series.png"), help="static figure file")
    ap.add_argument("--html", default=str(b.OUTPUTS / "thickness_series.html"), help="interactive plotly figure file")
    a = ap.parse_args()

    b.OUTPUTS.mkdir(exist_ok=True)
    # --plot-only reuses existing runs; otherwise (re)run the series first. The HTML is written last because it embeds
    # the animation and road-reach results if they already exist in Outputs/.
    if a.plot_only:
        runs = findRuns(a.thickness, a.rho, a.friction)
    else:
        if a.clean:
            shutil.rmtree(seriesDir(a.rho, a.friction), ignore_errors=True)
        runs = runSeries(a.thickness, a.rho, a.friction)
    plot(runs, a.rho, a.friction, a.out)
    plotInteractive(runs, a.rho, a.friction, a.html)


if __name__ == "__main__":
    main()
