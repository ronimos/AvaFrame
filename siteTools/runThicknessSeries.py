"""Run com1DFA for a site at several release thicknesses and plot peak flow velocity / thickness.

Usage (from the repo root; defaults for thickness, rho and friction come from site.toml):
    uv run python siteTools/runThicknessSeries.py --site McClure               # run all thicknesses, then plot
    uv run python siteTools/runThicknessSeries.py --site McClure --plot-only   # reuse existing runs
    uv run python siteTools/runThicknessSeries.py --site McClure --clean       # delete the series' run directories first
    uv run python siteTools/runThicknessSeries.py --site McClure --thickness 1 2 --rho 250 --friction samosATMedium
    uv run python siteTools/runThicknessSeries.py --site McClure --shrubs      # with the [shrubs] resistance area
                                                                               # -> Outputs/thickness_series_shrubs.*

Writes a static PNG and a self-contained interactive HTML report (linked canvas maps, see interactiveMaps.py) to Outputs/.
Each thickness is its own AvaFrame avalanche dir under runs/ (settings = Inputs/CFGs/com1DFACfg.ini plus the
overrides of that run), run in its own process.
The release shapefile (Inputs/REL) is used as is.
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
import geopandas as gpd
import rasterio

import avasite as b

footprintThr = 0.1  # m, pft threshold used to define the footprint / zoom extent
SHRUBS = False      # set by --shrubs: series with the shrub resistance area


def skey(k):
    """Sort key for series entries: named scenarios (str, e.g. "D2") first, then release thicknesses (float)."""
    return (0, k, 0) if isinstance(k, str) else (1, "", k)


def lab(k, html=True):
    """Label of a series entry: '0.5 m' for a thickness, the scenario id otherwise."""
    return k if isinstance(k, str) else f"{k:g}{'&nbsp;' if html else ' '}m"


def runName(k):
    return k if isinstance(k, str) else f"th_{k:g}"


def seriesDir(rho, friction, shrubs=None):
    shrubs = SHRUBS if shrubs is None else shrubs
    return b.RUNS / f"thickness_{friction}_rho{rho:g}{'_shrubs' if shrubs else ''}"


def runSeries(keys, rho, friction):
    """One AvaFrame run (own process, own avalanche dir) per release thickness or named scenario; returns {key: avaDir}."""
    runs = {}
    for k in keys:
        runs[k] = seriesDir(rho, friction) / runName(k)
        th = b.SCENARIOS[k]["thickness"] if isinstance(k, str) else k
        o = dict(rho=f"{rho:g}", frictModel=friction, relThFromFile=False, relTh=f"{th:g}")
        b.runCom1DFA(runs[k], {**o, **(b.shrubOverrides() if SHRUBS else {})}, res=SHRUBS,
                     rel=b.scenarioRel(k) if isinstance(k, str) else None)
        print(f"{lab(k, False)} done")
    return runs


def findRuns(keys, rho, friction, shrubs=None):
    """Existing run dirs {key: avaDir} of a series."""
    runs = {k: seriesDir(rho, friction, shrubs) / runName(k) for k in keys}
    missing = [t for t, d in runs.items() if not (d / "Outputs/com1DFA/peakFiles").is_dir()]
    if missing:
        sys.exit(f"no results found for relTh {sorted(missing, key=skey)}")
    return runs


def loadSeries(runs):
    """Peak rasters of all runs plus hillshade, release polygon and a zoom window around the flow paths."""
    th = sorted(runs, key=skey)
    pft, pfv, ppr = ({t: b.peakRaster(runs[t], f)[0] for t in th} for f in ("pft", "pfv", "ppr"))
    _, T, B = b.peakRaster(runs[th[0]], "pft")
    terrain = b.loadTerrain()
    rel = terrain["rel"]
    xl, yl = b.zoomWindow([pft[t] > footprintThr for t in th], T, rel, pad=120)
    return dict(th=th, pft=pft, pfv=pfv, ppr=ppr, T=T, meshExt=b.extent(B), rel=rel, hill=terrain["hill"],
                demExt=terrain["ext"], xl=xl, yl=yl)


# HTML skeleton of the report; {css}, {div} (linked maps), {hist} (closest-approach histogram), {anim} (video block) and
# {method} (methodology) are filled in plotInteractive()
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
{css}
</style></head><body>
<h1>{title}</h1>
<p class="sub">One map per release thickness. Hover for values; zoom, pan and the crosshair are linked across panels.</p>
{div}
{hist}
{anim}
<div class="method">{method}</div>
</body></html>"""


def runInfo(runs, rho):
    """Per thickness: friction class AvaFrame used, release mass [kg] and volume [m3] (initial mass / density)."""
    info = {}
    for t, tag in runs.items():
        mass = b.releaseMass(tag)
        info[t] = dict(cls=b.runConfig(tag)["GENERAL"]["frictModel"], mass=mass, vol=mass / rho, d=b.dSize(mass))
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
    """Embedded MP4 of the side-by-side animation, if it exists (the animation runs have no shrubs)."""
    import base64
    if SHRUBS:
        return ""
    d = b.OUTPUTS / f"animation_{friction}"
    srcs = "".join(f'<source src="data:{mime};base64,{base64.b64encode(f.read_bytes()).decode()}" type="{mime}">'
                   for f, mime in ((d / "flow.webm", "video/webm"), (d / "flow.mp4", "video/mp4")) if f.exists())
    if not srcs:
        return ""
    return f"""<h2>Animation</h2>
<p class="sub">Flow thickness for each release thickness on a shared clock; a run that has stopped holds its last frame.
 Saved every 2&nbsp;s of simulated time. The colour scale is the same in every panel and frame.
 Also saved as <code>flow.gif</code> and individual frames in <code>Outputs/animation_{friction}/</code>.</p>
<video controls loop muted autoplay playsinline style="width:100%;max-width:1300px">{srcs}Your browser cannot play this video; see flow.gif in the output folder.</video>"""


def frictionParagraph(info, dc, th, auto, friction):
    """Explain the friction choice from the actual classes/volumes of this run set."""
    if not auto:
        return (f"<p>All runs use <code>{friction}</code> (mu = {dc['mu'].get(friction, float('nan')):g}). It is the same for every depth, so "
                "differences between rows come from release thickness (and so volume) only.</p>")
    vs, vm = dc["volSmall"], dc["volMedium"]
    mu = dc["mu"]
    small = [t for t in th if info[t]["vol"] < vs]   # scenarios inside the samosATSmall volume class
    txt = f"""<p><b>How.</b> Friction was set with <code>samosATAuto</code>. AvaFrame computes each simulation's release volume and picks one of
 three calibrations of the samosAT friction model by volume: <code>samosATSmall</code> (&mu; = {mu['samosATSmall']:g}) below {vs:,.0f}&nbsp;m&sup3;,
 <code>samosATMedium</code> (&mu; = {mu['samosATMedium']:g}) from {vs:,.0f} to {vm:,.0f}&nbsp;m&sup3;, and standard <code>samosAT</code>
 (&mu; = {mu['samosAT']:g}) above {vm:,.0f}&nbsp;m&sup3;. The class is chosen per depth, not by hand; the table shows what each depth got
 (release volume = initial simulated mass / density).</p>
<p><b>Why.</b> The calibrations assume that larger avalanches are less resisted (lower &mu;). Release thickness changes the volume
 from {min(i['vol'] for i in info.values()):,.0f} to {max(i['vol'] for i in info.values()):,.0f}&nbsp;m&sup3;, so one fixed class would be mis-calibrated for some of
 the depths. In particular <code>samosATSmall</code> is meant for volumes below {vs:,.0f}&nbsp;m&sup3;; """ + (
        "no scenario here is that small, so it was not used. " if not small else
        f"only the {', '.join(lab(t) for t in small)} scenario{'s are' if len(small) > 1 else ' is'} below that and use{'s' if len(small) == 1 else ''} it. ") + """An earlier
 version of this analysis used <code>samosATSmall</code> for all depths, including ones far above that volume, and was replaced for that reason.</p>
<p><b>Caveat.</b> The class boundaries are step changes. Where a depth sits just on one side of a boundary, the friction jumps between neighbouring
 depths, so part of the difference between those rows is a friction change, not only a thickness change."""
    near = [t for t in th if min(abs(info[t]["vol"] - vs), abs(info[t]["vol"] - vm)) < 5000]
    if near:
        txt += " Closest to a boundary: " + ", ".join(
            f"{lab(t)} ({info[t]['vol']:,.0f}&nbsp;m&sup3; vs the {vm if abs(info[t]['vol'] - vm) < abs(info[t]['vol'] - vs) else vs:,.0f}&nbsp;m&sup3; limit)" for t in near) + "."
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
    th = sorted(runs, key=skey)
    road = b.hasRoad()
    mesh = float(b.runConfig(runs[th[0]])["GENERAL"]["meshCellSize"])
    try:
        ver = md.version("avaframe").rsplit(".d2", 1)[0]  # drop uv's dirty-tree date suffix
    except md.PackageNotFoundError:
        ver = "unknown"
    with rasterio.open(b.SOURCE_DEM) as r:
        srcRes = abs(r.transform.a)
    with rasterio.open(b.DEM) as r:
        demInfo = f"{r.width} x {r.height} cells at {abs(r.transform.a):g} m"
    adj = b.INPUTS / "release_adjustment.json"
    if adj.exists():
        import json
        a = json.loads(adj.read_text())
        release = f"""<p>The release area started from a manually drawn polygon ({a['areaOriginalHa']:.2f} ha). Its downslope
 extent was then set by slope: the zone follows terrain that is <b>&ge; {a['minSlopeDeg']:g}&deg; in the downhill direction</b>.</p>
<ol>
<li>Slope was computed from the USGS 1&nbsp;m DEM after averaging to {a['blockCellM']}&nbsp;m blocks, using a finite-difference gradient and a
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
        release = (f"<p>The release area is the polygon <code>{pathlib.Path(b.CFG['release']).name}</code> "
                   f"({gpd.read_file(b.RELEASE).area.sum() / 1e4:.2f}&nbsp;ha planimetric), used as drawn. "
                   "Release thickness is uniform over the zone and is the varied parameter.</p>")
    release += b.CFG.get("report", {}).get("releaseNote", "")
    info, dc = runInfo(runs, rho), defaultCfg()
    rows = "".join(f"<tr><td>{lab(t)}</td><td>{info[t]['vol']:,.0f} m&sup3;</td><td>{info[t]['mass'] / 1000:,.0f} t</td><td>D{info[t]['d']:g}</td><td><code>{info[t]['cls']}</code></td>"
                   f"<td>{dc['mu'][info[t]['cls']]:g}</td><td><code>{runs[t].relative_to(b.BASE)}</code></td></tr>" for t in th)
    auto = friction == "samosATAuto"
    frictionText = frictionParagraph(info, dc, th, auto, friction)
    return f"""
<h2>Methodology</h2>
<h3>Terrain</h3>
<p>Elevation: <code>{b.SOURCE_DEM.name}</code> (EPSG:{b.CRS}, {srcRes * 100:.2f}&nbsp;cm cells),
 block-averaged to a 1&nbsp;m grid ({demInfo}). Outside the surveyed area the grid is nodata, so the surveyed area is the simulation
 domain: flow that leaves it is removed. AvaFrame then remeshes the DEM to a {mesh:g}&nbsp;m computational mesh, so detail finer than
 that is averaged before the solver runs. The hillshade in the maps is from the 1&nbsp;m grid.</p>
{b.CFG.get("report", {}).get("terrainNote", "")}
{filterHtml()}
<h3>Release area</h3>
{release}
<h3>Simulations</h3>
<p>Dense-flow avalanche simulations with AvaFrame <code>com1DFA</code> (version {ver}). Constant for all runs: snow density
 {rho:g}&nbsp;kg/m&sup3;, {mesh:g}&nbsp;m computational mesh, release area as above, all other parameters at AvaFrame defaults{" except the shrub settings below" if SHRUBS else ""}.
 The only input that differs between rows is the release thickness ({", ".join(f"{t:g}" for t in th if not isinstance(t, str))}&nbsp;m); each value is a separate simulation.{scenarioHtml(th)}
 Friction is set per run by AvaFrame from the release volume (next section).</p>
{shrubHtml()}
<h3>Friction model: how and why each depth got its class</h3>
{frictionText}
<table><tr><th>Release</th><th>Release volume</th><th>Mass</th><th>Size (by mass)</th><th>Friction class used</th><th>&mu;</th><th>Run directory</th></tr>{rows}</table>
<h3>Outputs shown</h3>
<p>Peak flow thickness (m), velocity (m/s) and pressure (kPa): the maximum reached at each cell during the simulation, not a
 snapshot. Cells are shown where peak thickness exceeds {footprintThr:g}&nbsp;m. Colour scales are shared across panels. The readout
 gives each run's values at the cursor, its easting/northing (EPSG:{b.CRS}), and how many runs reach that cell.</p>
""" + ("" if not road else f"""<p><b>Road reference.</b> "Reaches" and "crosses" mean the flow touches the hill-side and far-side road edge lines (same test as
 <code>roadReach.py</code>). Distances in the maps and the histogram are measured to one point: where the avalanche path meets the
 hill-side edge. The path is traced by steepest descent (D8, on the pit-filled {mesh:g}&nbsp;m DEM) from the highest cell reached by at
 least half of the runs. The trace stops 40&nbsp;m from the road, because closer in it would turn along the road-cut ditch, and its
 direction over the last 100&nbsp;m is projected onto the road. The "m short of the road edge" in the panel titles is instead the gap to
 the nearest point of the edge line, so the two distances differ.</p>""") + f"""
<h3>Limitations</h3>
<ul>
<li>Results are from a single deterministic model setup and have not been calibrated or validated against observed events at this site.</li>
<li>Friction classes come from volume bands and change in steps; a class change between depths is a modelling choice, not a measured difference.</li>
<li>Release thickness is a scenario input, not a measurement""" + (f"; slope and release-area extent depend on the 5&nbsp;m averaging and the {a['minSlopeDeg']:g}&deg; threshold" if adj.exists() else "; it is uniform over the whole release polygon") + """.</li>
""" + "".join(f"<li>{x}</li>" for x in b.CFG.get("report", {}).get("limitations", [])) + f"""
</ul>
<p class="sub">Generated {datetime.date.today().isoformat()} by <code>siteTools/runThicknessSeries.py --site {b.SITE}</code>.</p>"""


def scenarioHtml(keys):
    """Sentences on the named scenarios (e.g. D2) and the size estimate, for the Simulations section."""
    txt = """ The <b>size</b> column is the destructive size implied by the release mass alone, using the typical masses of the
 scale (D1 10&nbsp;t, D2 100&nbsp;t, D3 1,000&nbsp;t, D4 10,000&nbsp;t), rounded to the nearest half size. Real size classes also
 depend on destructive potential, so this is a rough guide."""
    for k in keys:
        if isinstance(k, str):
            sc = b.SCENARIOS[k]
            txt += (f" The <b>{k}</b> row releases a smaller polygon, <code>{pathlib.Path(sc['release']).name}</code>, with "
                    f"{float(sc['thickness']):g}&nbsp;m thickness: {sc.get('description', '')}")
    return txt


def filterHtml():
    """Vegetation filter paragraph from Inputs/ground_filter.json (prepInputs.py), if the filter was used."""
    f = b.INPUTS / "ground_filter.json"
    if not f.exists():
        return ""
    import json
    r = json.loads(f.read_text())
    return f"""<p><b>Vegetation removed.</b> Tree and shrub crowns were filtered out of the 1&nbsp;m grid before the runs
 (<code>siteTools/groundFilter.py</code>): grey-scale openings with windows of {", ".join(f"{w:g}" for w in r["windows"])}&nbsp;m flag a cell as
 vegetation when it stands more than {r["dh0"]:g}&nbsp;m + {r["slope"]:g}&nbsp;&times;&nbsp;window above the opened surface (the window term keeps
 convex ground such as ridges and gully banks); flagged cells are filled by interpolation from the surrounding ground. {r["vegetationFraction"]:.0%} of the
 cells were flagged; the height removed there has a median of {r["removedMedianM"]:g}&nbsp;m and a 95th percentile of {r["removedP95M"]:g}&nbsp;m
 (maximum {r["removedMaxM"]:g}&nbsp;m). Low shrubs (about 1&ndash;2&nbsp;m) are partly left in the surface. The filter is a simple
 morphological one, not a lidar ground classification; the open gully floor, which has little vegetation, changed by a few centimetres.</p>"""


def shrubHtml():
    """Shrub resistance paragraph for the --shrubs series."""
    if not SHRUBS:
        return ""
    import configparser
    c = configparser.ConfigParser(interpolation=None)
    c.read(b.REPO / "avaframe/com1DFA/com1DFACfg.ini")
    g, h = c["GENERAL"], float(b.SHRUBS["height"])
    excl = b.SHRUBS.get("excludeRelease", False)
    return f"""<h3>Shrubs (resistance area)</h3>
<p>Shrubs about 1&ndash;{h:g}&nbsp;m high cover the whole path. In this series they are a com1DFA resistance area covering the whole
 surveyed domain except the road{" and the release polygon" if excl else ""}. Where flow in that area is faster than {float(g['forestVMin']):g}&nbsp;m/s and deeper than
 {float(g['forestThMin']):g}&nbsp;m, an extra velocity-squared drag is applied (<code>cResH</code> = {float(g['cResH']):g}&nbsp;1/m); where it is slower
 or shallower, snow is taken out of the flow (detrainment, <code>detK</code> = {float(g['detK']):g}). Flow deeper than the shrub height
 ({h:g}&nbsp;m, <code>forestThMax</code>) or faster than {float(g['forestVMax']):g}&nbsp;m/s is taken to override the shrubs, which then stop acting
 in that cell for the rest of the run. Apart from <code>forestThMax</code>, these are AvaFrame's default values, which were set up for forest,
 not shrubs, and are marked in AvaFrame as still needing calibration. The shrubs are treated as standing above the snow; a snowpack
 that buries them would remove most of their effect. Compare with the series without shrubs.</p>""" + ("""
<p>The release polygon is left out because snow starting from rest is slower than the detrainment threshold: with shrubs there,
 com1DFA takes thin releases out of the flow before they move. That would express shrubs anchoring the slab, which is a question of
 whether it releases at all, not of how far it runs; the release thickness scenarios assume it has released.</p>""" if excl else "")


def roadStatus(pft, hillL, farL, T, thr):
    """Same test as roadReach.py: flow (pft > thr) touching the hill-side / far-side edge line, else gap to the hill-side edge."""
    from rasterio.features import rasterize
    from shapely import distance, points
    flow = pft > thr
    hit = lambda g: bool((flow & rasterize([g], out_shape=flow.shape, transform=T, all_touched=True).astype(bool)).any())
    if hit(farL):
        return "<b>crosses the road</b>"
    if hit(hillL):
        return "<b>reaches the road</b>"
    rr, cc = np.nonzero(flow)
    x, y = T * (cc + .5, rr + .5)
    return f"stops {distance(points(x, y), hillL).min():.0f} m short of the road edge"


def roadGroups(runs, friction, rho):
    """Histogram groups: this series plus the roadReach.py runs (fixed friction classes) if they exist."""
    import roadReach as rr
    label = lambda shrubs: f"{friction}{' + shrubs' if shrubs else ''}"
    groups = [dict(name=label(SHRUBS), runs=[dict(id=lab(t, False), pft=b.peakRaster(runs[t], "pft")[0]) for t in sorted(runs, key=skey)])]
    if b.SHRUBS:  # the other variant of the same series, if it has been run
        other = {t: seriesDir(rho, friction, not SHRUBS) / runName(t) for t in sorted(runs, key=skey)}
        if all((d / "Outputs/com1DFA/peakFiles").is_dir() for d in other.values()):
            groups.append(dict(name=label(not SHRUBS), runs=[dict(id=lab(t, False), pft=b.peakRaster(d, "pft")[0]) for t, d in other.items()]))
    for fr in rr.FRICT:
        dirs = {t: rr.runDir(fr, t) for t in rr.THICK}
        if all((d / "Outputs/com1DFA/peakFiles").is_dir() for d in dirs.values()):
            groups.append(dict(name=fr, runs=[dict(id=f"{t:g} m", pft=b.peakRaster(d, "pft")[0]) for t, d in dirs.items()]))
    return groups


def plotInteractive(runs, rho, friction, fname, png):
    """Self-contained HTML: linked canvas maps per thickness (interactiveMaps.py), animation, methodology."""
    import base64
    import interactiveMaps as im

    s = loadSeries(runs)
    th, T = s["th"], s["T"]
    info = runInfo(runs, rho)
    hillL, farL = b.roadLines()
    panels = [dict(id=str(t) if isinstance(t, str) else f"{t:g}", title=f"Release {lab(t)}",
                   sub=f"(D{info[t]['d']:g} by mass, {info[t]['cls']}, {info[t]['vol']:,.0f} m&sup3;)",
                   status=roadStatus(s["pft"][t], hillL, farL, T, footprintThr),
                   fields={k: s[k][t] for k in ("pft", "pfv", "ppr")}) for t in th]
    fields = [dict(key="pft", name="Peak flow thickness", unit="m", scale=100, dec=2),
              dict(key="pfv", name="Peak flow velocity", unit="m/s", scale=100, dec=1),
              dict(key="ppr", name="Peak flow pressure", unit="kPa", scale=10, dec=0)]
    dem, demT, _ = b.readRaster(b.DEM)
    printPng = "data:image/png;base64," + base64.b64encode(pathlib.Path(png).read_bytes()).decode()
    css, maps, hist = im.section(printPng, panels=panels, T=T, fields=fields, thr=footprintThr, roads=[hillL, farL],
                                 roadRef=hillL, release=list(s["rel"].geometry), dem1=dem, demT=demT,
                                 groups=roadGroups(runs, friction, rho))
    pathlib.Path(fname).write_text(PAGE.format(
        title=f"{b.NAME} avalanche: peak flow by release thickness{' (with shrub resistance)' if SHRUBS else ''}", css=css, div=maps, hist=hist, anim=animationHtml(friction),
        method=methodologyHtml(runs, rho, friction)))
    print("saved", fname, f"({pathlib.Path(fname).stat().st_size / 1e6:.1f} MB)")


def plot(runs, rho, friction, fname):
    s = loadSeries(runs)
    info = runInfo(runs, rho)
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
            ax.set_title(f"{name} - release {lab(t, False)} (D{info[t]['d']:g})", fontsize=10, loc="left")
            ax.tick_params(labelsize=7); ax.ticklabel_format(style="plain", useOffset=False)
            fig.colorbar(im, ax=ax, shrink=.75, label=unit)
    fig.suptitle(f"{b.NAME} com1DFA: rho {rho:g}, {friction}{' + shrub resistance' if SHRUBS else ''} (orange = release area; colour scales shared per column)", fontsize=11)
    fig.savefig(fname, dpi=100)
    print("saved", fname)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--thickness", type=float, nargs="+", default=b.THICK, help="release thicknesses [m]")
    ap.add_argument("--rho", type=float, default=b.RHO, help="snow density [kg/m3]")
    ap.add_argument("--friction", default=b.FRICTION, help="com1DFA frictModel (samosATAuto picks the class per release volume)")
    ap.add_argument("--clean", action="store_true", help="delete the series' run directories before running")
    ap.add_argument("--plot-only", action="store_true", help="skip simulations, only (re)plot existing results")
    ap.add_argument("--shrubs", action="store_true", help="add the site.toml [shrubs] resistance area (separate series and report)")
    ap.add_argument("--out", help="static figure file (default Outputs/thickness_series[_shrubs].png)")
    ap.add_argument("--html", help="interactive HTML report (default Outputs/thickness_series[_shrubs].html)")
    a = ap.parse_args()
    global SHRUBS
    SHRUBS = a.shrubs
    if SHRUBS and not b.SHRUBS:
        sys.exit("--shrubs needs a [shrubs] table in site.toml (and prepInputs.py run after adding it)")
    suffix = "_shrubs" if SHRUBS else ""
    a.out = a.out or str(b.OUTPUTS / f"thickness_series{suffix}.png")
    a.html = a.html or str(b.OUTPUTS / f"thickness_series{suffix}.html")

    b.OUTPUTS.mkdir(exist_ok=True)
    # --plot-only reuses existing runs; otherwise (re)run the series first. The HTML is written last because it embeds
    # the animation and road-reach results if they already exist in Outputs/.
    keys = list(b.SCENARIOS) + list(a.thickness)        # named scenarios ([[scenario]] in site.toml) + thicknesses
    if a.plot_only:
        runs = findRuns(keys, a.rho, a.friction)
    else:
        if a.clean:
            shutil.rmtree(seriesDir(a.rho, a.friction), ignore_errors=True)
        runs = runSeries(keys, a.rho, a.friction)
    plot(runs, a.rho, a.friction, a.out)
    plotInteractive(runs, a.rho, a.friction, a.html, a.out)


if __name__ == "__main__":
    main()
