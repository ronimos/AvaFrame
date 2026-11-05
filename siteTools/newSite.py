"""Set up a new avalanche path from its files, check them, and build a draft report.

Usage (from the repo root):
    # files given directly (they are copied into <Name>/data/):
    uv run python siteTools/newSite.py Loveland --title "Loveland Pass (No Name)" \\
        --dem survey.tif --release start_zone.kml --road edge_a.kml edge_b.kmz [--shrubs 1.5] [--d2 d2.kml]
    # or drop the files into <Name>/data/DEM, <Name>/data/release_area and <Name>/data/Road first, then:
    uv run python siteTools/newSite.py Loveland --title "Loveland Pass (No Name)"

What it does:
  1. copies the files into <Name>/data/ and writes <Name>/site.toml and <Name>/.gitignore
  2. checks the inputs (see CHECKS below); a FAIL stops here, WARNs are printed and kept in <Name>/checks.json
  3. if the surface model still holds vegetation, turns on the vegetation filter
  4. adds a D2 scenario: --d2 <polygon> if given, else an assumed ~100 t slab at the top of the start zone (--no-d2: none)
  5. runs the whole pipeline (runSite.py) -> draft reports in <Name>/Outputs/ (--no-run: stop after the checks)
Nothing is published: look at the draft, then run  uv run python siteTools/publish.py --site <Name>

CHECKS: DEM has a projected CRS in metres; release and road files read as one polygon / one line each; the release
polygon and both road edges lie on valid DEM cells; the DEM extends past the road; vegetation bumps in the surface;
which road edge is on the hill side (the one nearer the release); travel angle from the top of the start zone to the
road and the size class of each release thickness.
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
from pyproj import CRS
from rasterio.features import rasterize
from scipy import ndimage
from shapely.geometry import Point

REPO = Path(__file__).resolve().parent.parent
HERE = REPO / "siteTools"
VECTOR = (".kml", ".kmz", ".geojson", ".json", ".shp", ".gpkg")
RASTER = (".tif", ".tiff")
THICK = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
RHO = 300.0
VEG_P95_M = 2.5          # start-zone bumps (95th percentile above a 15 m median surface) above this = vegetation in the DEM

results = []             # (level, check, message)


def note(level, check, msg):
    results.append(dict(level=level, check=check, message=msg))
    print(f"  [{level}] {check}: {msg}")


def fail_if_any():
    if any(r["level"] == "FAIL" for r in results):
        sys.exit("\nstopped: fix the FAIL items above and run again")


def readVector(path):
    path = Path(path)
    if path.suffix.lower() == ".kmz":
        import zipfile
        with zipfile.ZipFile(path) as z:
            inner = next(n for n in z.namelist() if n.lower().endswith(".kml"))
        return gpd.read_file(f"/vsizip/{path}/{inner}")
    return gpd.read_file(path)


def place(src, dstDir):
    """Copy (hard-link when possible, for big DEMs) src into dstDir; a file already there is used as is."""
    src = Path(src).resolve()
    dstDir.mkdir(parents=True, exist_ok=True)
    dst = dstDir / src.name
    if src.parent == dstDir.resolve():
        return dst
    if dst.exists():
        sys.exit(f"{dst} already exists; remove it or use the file that is there")
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)
    if src.suffix.lower() == ".shp":  # shapefile sidecars
        for side in src.parent.glob(src.stem + ".*"):
            if side.suffix.lower() != ".shp" and not (dstDir / side.name).exists():
                shutil.copy2(side, dstDir / side.name)
    return dst


def found(folder, exts):
    return sorted(p for p in folder.glob("*") if p.suffix.lower() in exts) if folder.is_dir() else []


def tomlStr(s):
    return json.dumps(s, ensure_ascii=False)  # a JSON string is a valid TOML basic string


def writeToml(site, cfg):
    """site.toml from the collected settings (same layout as McClure/site.toml)."""
    L = [f"# Written by siteTools/newSite.py. Run: uv run python siteTools/runSite.py --site {site.name}. Paths are relative to this folder.",
         f"title = {tomlStr(cfg['title'])}",
         f"crs = {cfg['crs']}                       # EPSG of the DEM (horizontal part); vector files are reprojected to it",
         f"sourceDem = {tomlStr(cfg['sourceDem'])}",
         "demCell = 1.0",
         f"release = {tomlStr(cfg['release'])}",
         f"releaseName = {tomlStr(cfg['releaseName'])}",
         f"thicknesses = {json.dumps(cfg['thicknesses'])}   # m",
         f"rho = {cfg['rho']:g}                         # kg/m3",
         'friction = "samosATAuto"',
         'roadFrictions = ["samosATMedium", "samosAT"]',
         "",
         "[road]",
         f"hillSide = {tomlStr(cfg['hillSide'])}     # edge nearer the release (chosen by newSite.py)",
         f"farSide = {tomlStr(cfg['farSide'])}"]
    if cfg.get("vegetationFilter"):
        L += ["", "# the surface model still holds tree/shrub crowns (newSite.py check); see siteTools/groundFilter.py",
              "[vegetationFilter]", "windows = [3, 5, 9, 15, 21]", "dh0 = 0.5", "slope = 0.15"]
    if cfg.get("shrubs"):
        L += ["", "[shrubs]", f"height = {cfg['shrubs']:g}                      # m; flow deeper than this overrides the shrubs",
              "excludeRelease = true"]
    if cfg.get("terrainNote"):
        L += ["", "[report]", f'terrainNote = """{cfg["terrainNote"]}"""']
    for sc in cfg.get("scenarios", []):
        L += ["", "[[scenario]]", f"id = {tomlStr(sc['id'])}", f"release = {tomlStr(sc['release'])}",
              f"thickness = {sc['thickness']:g}", f'description = """{sc["description"]}"""']
    (site / "site.toml").write_text("\n".join(L) + "\n")


def prep(site):
    r = subprocess.run([sys.executable, str(HERE / "prepInputs.py"), "--site", site.name], cwd=REPO, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"prepInputs.py failed:\n{r.stdout}\n{r.stderr}")


def readDem(site):
    with rasterio.open(site / f"Inputs/{site.name}_DEM.asc") as r:
        z = r.read(1).astype(float)
        z[z == r.nodata] = np.nan
        return z, r.transform


def bumps(z):
    """Height above a 15 m median surface: tree/shrub crowns stand out, smooth ground does not."""
    f = np.where(np.isnan(z), np.nanmean(z), z)
    return f - ndimage.median_filter(f, 15)


def assumedD2(z, T, relPoly, crs, th=0.5, massT=100.0):
    """~massT tonne slab of thickness th at the top of the start zone, centred on steep (>= 30 deg) ground."""
    m = rasterize([relPoly], out_shape=z.shape, transform=T).astype(bool)
    zs = ndimage.uniform_filter(np.where(np.isnan(z), np.nanmean(z), z), 5)
    gy, gx = np.gradient(zs, abs(T.e), T.a)
    slope = ndimage.uniform_filter(np.degrees(np.arctan(np.hypot(gx, gy))), 25)
    for minSlope in (35, 30, 0):
        rad = math.sqrt(massT * 1000 / RHO / th * math.cos(math.radians(max(minSlope, 30))) / math.pi)
        cand = (ndimage.distance_transform_edt(m) * T.a >= rad) & (slope >= minSlope)
        if cand.any():
            break
    else:
        return None
    r, c = np.unravel_index(np.argmax(np.where(cand, z, -np.inf)), z.shape)
    x, y = T * (c + .5, r + .5)
    rad = math.sqrt(massT * 1000 / RHO / th * math.cos(math.radians(slope[r, c])) / math.pi)
    geom = Point(x, y).buffer(rad, 32).intersection(relPoly)
    return dict(geom=geom, x=x, y=y, z=float(z[r, c]), slope=float(slope[r, c]), radius=rad)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("name", help="site folder name in the repo root (letters, digits, _)")
    ap.add_argument("--title", help="name used in report titles (default: the folder name)")
    ap.add_argument("--dem", help="surface model GeoTIFF")
    ap.add_argument("--release", help="start zone polygon (kml, kmz, geojson, shp)")
    ap.add_argument("--road", nargs=2, metavar=("EDGE1", "EDGE2"), help="the two road edge lines, in any order")
    ap.add_argument("--shrubs", type=float, metavar="HEIGHT_M", help="add the shrub resistance series (shrub height, m)")
    ap.add_argument("--d2", help="outline of a typical D2 release (polygon); default: an assumed slab at the top of the start zone")
    ap.add_argument("--no-d2", action="store_true", help="no D2 scenario")
    ap.add_argument("--thickness", type=float, nargs="+", default=THICK, help="release thicknesses [m]")
    ap.add_argument("--no-run", action="store_true", help="stop after the checks and input preparation")
    a = ap.parse_args()

    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", a.name):
        sys.exit("name: letters, digits and _ only, starting with a letter")
    site = REPO / a.name
    if (site / "site.toml").exists():
        sys.exit(f"{site / 'site.toml'} exists already; edit it and use runSite.py, or remove it to start over")
    data = site / "data"
    print(f"== {a.name}: files")

    # 1) files: given ones are copied in; otherwise found in the data folders
    dem = place(a.dem, data / "DEM") if a.dem else None
    rel = place(a.release, data / "release_area") if a.release else None
    roads = [place(p, data / "Road") for p in a.road] if a.road else None
    d2 = place(a.d2, data / "release_area") if a.d2 else None
    if dem is None:
        c = found(data / "DEM", RASTER)
        if len(c) != 1:
            sys.exit(f"put exactly one DEM (.tif) in {data / 'DEM'} or use --dem (found {len(c)})")
        dem = c[0]
    if rel is None:
        c = [p for p in found(data / "release_area", VECTOR) if p != d2 and "d2" not in p.stem.lower()]
        if len(c) != 1:
            sys.exit(f"put exactly one start zone polygon in {data / 'release_area'} or use --release (found {[p.name for p in c]})")
        rel = c[0]
    if roads is None:
        roads = found(data / "Road", VECTOR)
        if len(roads) != 2:
            sys.exit(f"put the two road edge lines in {data / 'Road'} or use --road (found {[p.name for p in roads]})")
    print(f"  DEM {dem.name}\n  start zone {rel.name}\n  road edges {roads[0].name}, {roads[1].name}")

    # 2) checks that need no simulation grid
    print(f"\n== {a.name}: checks")
    with rasterio.open(dem) as r:
        c = CRS(r.crs.to_wkt()) if r.crs else None
        bounds = r.bounds
        res = abs(r.res[0])
    if c is None:
        note("FAIL", "DEM CRS", "the DEM has no coordinate system")
        fail_if_any()
    h = c.sub_crs_list[0] if c.is_compound else c
    epsg = h.to_epsg()
    if not h.is_projected or h.axis_info[0].unit_name not in ("metre", "meter"):
        note("FAIL", "DEM CRS", f"{h.name} is not a projected CRS in metres; reproject the DEM (e.g. to UTM)")
    elif epsg is None:
        note("FAIL", "DEM CRS", f"{h.name} has no EPSG code")
    else:
        note("OK", "DEM CRS", f"EPSG:{epsg} ({h.name}), {res * 100:.1f} cm cells")
    fail_if_any()
    try:
        relG = readVector(rel).to_crs(epsg)
        relPoly = relG.geometry.union_all()
        relPoly = relPoly if relPoly.geom_type == "Polygon" else max(relPoly.geoms, key=lambda g: g.area)
        note("OK" if len(relG) == 1 else "WARN", "start zone", f"{relPoly.area / 1e4:.2f} ha planimetric"
             + ("" if len(relG) == 1 else f"; {len(relG)} features in the file, the largest is used"))
    except Exception as e:  # noqa: BLE001
        note("FAIL", "start zone", f"cannot read a polygon from {rel.name}: {e}")
    lines = []
    for p in roads:
        try:
            g = readVector(p).to_crs(epsg).geometry.union_all()
            if g.geom_type == "MultiLineString":
                from shapely.ops import linemerge
                g = linemerge(g)
            if g.geom_type != "LineString":
                note("FAIL", "road edge", f"{p.name} is a {g.geom_type}, expected one line")
            lines.append(g)
        except Exception as e:  # noqa: BLE001
            note("FAIL", "road edge", f"cannot read a line from {p.name}: {e}")
    fail_if_any()
    from shapely.geometry import box
    if not box(*bounds).contains(relPoly.union(lines[0]).union(lines[1])):
        note("FAIL", "coverage", "the start zone or a road edge lies outside the DEM")
    fail_if_any()
    # hill side = the edge nearer the start zone, compared across the road: for points along edge 0, the start-zone
    # distance of that point vs of the nearest point on edge 1 (so different edge lengths do not matter)
    votes = []
    for i in np.linspace(0, lines[0].length, 25):
        p0 = lines[0].interpolate(i)
        p1 = lines[1].interpolate(lines[1].project(p0))
        if p0.distance(p1) < 50:                         # only where the two edges face each other
            votes.append(relPoly.distance(p0) - relPoly.distance(p1))
    if not votes:
        note("FAIL", "road sides", "the two road edges are nowhere within 50 m of each other")
        fail_if_any()
    hi = 0 if np.median(votes) < 0 else 1
    md = [relPoly.distance(ln) for ln in lines]
    note("OK", "road sides", f"hill side = {roads[hi].name}, far side = {roads[1 - hi].name} "
         f"(the hill side is {abs(np.median(votes)):.1f} m nearer the start zone across the road, {len(votes)} points)")
    if abs(np.median(votes)) < 2:
        note("WARN", "road sides", "the two edges are about equally far from the start zone; check hillSide/farSide in site.toml")

    # 3) site.toml, then the 1 m DEM for the grid checks
    rp = lambda p: str(Path(p).resolve().relative_to(site.resolve()))
    cfg = dict(title=a.title or a.name, crs=epsg, sourceDem=rp(dem), release=rp(rel),
               releaseName=re.sub(r"\W", "_", Path(rel).stem), thicknesses=a.thickness, rho=RHO,
               hillSide=rp(roads[hi]), farSide=rp(roads[1 - hi]), shrubs=a.shrubs)
    writeToml(site, cfg)
    gi = (REPO / "McClure/.gitignore").read_text().replace("McClure", a.name)
    (site / ".gitignore").write_text(gi)
    print("  preparing the 1 m DEM ...")
    prep(site)
    z, T = readDem(site)
    ok = ~np.isnan(z)
    mRel = rasterize([relPoly], out_shape=z.shape, transform=T).astype(bool)
    cov = ok[mRel].mean()
    note("OK" if cov > .99 else "FAIL", "start zone on DEM", f"{cov:.0%} of the start zone has DEM data")
    for ln, p in zip(lines, roads):
        mL = rasterize([ln], out_shape=z.shape, transform=T, all_touched=True).astype(bool)
        cv = ok[mL].mean()
        note("OK" if cv > .95 else ("WARN" if cv > .5 else "FAIL"), "road on DEM", f"{cv:.0%} of {p.name} has DEM data"
             + ("" if cv > .95 else "; road hits are only detected where there is DEM"))
    # cells within 150 m of the far edge that are farther from the start zone than that edge is (= past the road)
    near = rasterize([lines[1 - hi].buffer(150)], out_shape=z.shape, transform=T).astype(bool)
    dRel = ndimage.distance_transform_edt(~mRel) * T.a
    beyond = near & (dRel > relPoly.distance(lines[1 - hi]) + 10)
    bv = ok[beyond].mean() if beyond.any() else 0
    note("OK" if bv > .6 else "WARN", "runout room", f"{bv:.0%} of the 150 m around the far road edge has DEM data"
         + ("" if bv > .6 else "; flow that leaves the DEM is lost, so runout past the road may be cut short"))
    fail_if_any()

    # 4) vegetation in the surface model
    bz = bumps(z)
    p95, p99 = np.percentile(bz[mRel & ok], [95, 99])
    veg = p95 > VEG_P95_M
    note("WARN" if veg else "OK", "vegetation", f"start zone bumps {p95:.1f} m (95th pct) / {p99:.1f} m (99th) above a 15 m median surface"
         + ("; looks like tree/shrub crowns, vegetation filter turned on" if veg else ""))
    if veg:
        cfg["vegetationFilter"] = True
        cfg["terrainNote"] = (f"<p><b>Vegetation in the surface.</b> The source surface shows tree or shrub crowns: relative to a 15&nbsp;m "
                              f"median-filtered surface the start zone has bumps of +{p95:.1f}&nbsp;m (95th percentile) to +{p99:.1f}&nbsp;m (99th). "
                              "They were removed before the runs (next paragraph).</p>")
        writeToml(site, cfg)
        print("  re-preparing the DEM with the vegetation filter ...")
        prep(site)
        z, T = readDem(site)
        b2 = np.percentile(bumps(z)[mRel & ok], 95)
        note("OK" if b2 < VEG_P95_M else "WARN", "vegetation filter", f"start zone bumps now {b2:.1f} m (95th pct)")

    # 5) D2 scenario
    if not a.no_d2:
        if d2:
            cfg["scenarios"] = [dict(id="D2", release=rp(d2), thickness=0.5,
                                     description=f"a typical D2 outline (<code>{d2.name}</code>), 0.5&nbsp;m thick.")]
            note("OK", "D2", f"from {d2.name}")
        else:
            s = assumedD2(z, T, relPoly, epsg)
            if s is None:
                note("WARN", "D2", "no room for an assumed D2 slab in the start zone; none added")
            else:
                f = data / "release_area" / f"{a.name}_D2_assumed.geojson"
                gpd.GeoDataFrame({"name": ["D2_assumed"]}, geometry=[s["geom"]], crs=f"EPSG:{epsg}").to_file(f, driver="GeoJSON")
                cfg["scenarios"] = [dict(id="D2", release=rp(f), thickness=0.5, description=(
                    f"an assumed D2-sized slab (about 100&nbsp;t): a {s['radius']:.1f}&nbsp;m-radius patch at the top of the start zone "
                    f"({s['z']:,.0f}&nbsp;m, about {s['slope']:.0f}&deg;). Its location is an assumption, not an observed D2 outline."))]
                note("OK", "D2", f"assumed slab at the top of the start zone ({s['z']:.0f} m, {s['slope']:.0f} deg), "
                     f"radius {s['radius']:.1f} m -> {f.name}; replace it with a real D2 outline if you have one")
        writeToml(site, cfg)

    # 6) context: travel angle to the road and size per thickness
    r0, c0 = np.unravel_index(np.nanargmax(np.where(mRel, z, np.nan)), z.shape)
    x0, y0 = T * (c0 + .5, r0 + .5)
    hill = lines[hi]
    pr = hill.interpolate(hill.project(Point(x0, y0)))
    cc, rr = ~T * (pr.x, pr.y)            # inverse transform gives (col, row)
    zr = z[int(rr), int(cc)]
    alpha = math.degrees(math.atan((z[r0, c0] - zr) / math.hypot(pr.x - x0, pr.y - y0)))
    note("INFO", "travel angle", f"top of start zone {z[r0, c0]:.0f} m to the nearest road point {zr:.0f} m: alpha {alpha:.1f} deg "
         "(straight line; the path is longer, so the true angle is lower)")
    gy, gx = np.gradient(np.where(ok, z, np.nanmean(z)), abs(T.e), T.a)
    surf = (np.sqrt(1 + gx ** 2 + gy ** 2)[mRel & ok]).sum() * T.a * abs(T.e)
    sizes = ", ".join(f"{t:g} m D{min(5, max(1, round(2 * math.log10(surf * t * RHO / 1000)) / 2)):g}" for t in a.thickness)
    note("INFO", "size by mass", f"whole start zone ({surf / 1e4:.2f} ha on the slope): {sizes}")

    prep(site)  # final inputs (scenarios, shrubs)
    (site / "checks.json").write_text(json.dumps(results, indent=1))
    print(f"\n  checks saved to {site / 'checks.json'}; settings in {site / 'site.toml'}")
    if a.no_run:
        print(f"\nnext: uv run python siteTools/runSite.py --site {a.name} --skip-prep")
        return
    print(f"\n== {a.name}: running the pipeline (about 10-20 minutes) ...")
    subprocess.run([sys.executable, str(HERE / "runSite.py"), "--site", a.name, "--skip-prep"], cwd=REPO, check=True)
    print(f"\ndraft reports: {site / 'Outputs'}/thickness_series*.html\n"
          f"look them over, then publish:  uv run python siteTools/publish.py --site {a.name}")


if __name__ == "__main__":
    main()
