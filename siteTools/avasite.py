"""Shared helpers for the site scripts: which site, its paths and settings, raster/terrain loading, running com1DFA.

Every script in siteTools/ works on one site folder in the repo root (Bethel/, McClure/, ...). Pick it with
`--site <Name>` on the command line (removed from sys.argv here, so the scripts' own argparse never sees it) or the
AVA_SITE environment variable. The site's settings are in <Site>/site.toml (see McClure/site.toml for all keys).

Site folder layout:
    site.toml  settings: CRS, source files, thicknesses, density, report text
    data/      raw source data (DEM tif, release polygon, road edge lines; geojson, kml, kmz or shp)
    Inputs/    AvaFrame inputs (DEM, REL, CFGs/com1DFACfg.ini = baseline settings for this site)
    runs/      one AvaFrame avalanche dir per simulation (generated)
    Outputs/   figures, tables, animations, report (generated)
"""
import configparser
import os
import pathlib
import shutil
import subprocess
import sys
import tomllib
import zipfile

import geopandas as gpd
import numpy as np
import rasterio
from matplotlib.colors import LightSource


def _siteName():
    if "--site" in sys.argv:
        i = sys.argv.index("--site")
        name = sys.argv[i + 1]
        del sys.argv[i:i + 2]
        os.environ["AVA_SITE"] = name  # so child scripts started from here get the same site
        return name
    if os.environ.get("AVA_SITE"):
        return os.environ["AVA_SITE"]
    sys.exit("which site? add --site <Name> (a folder in the repo root with a site.toml)")


REPO = pathlib.Path(__file__).resolve().parent.parent
SITE = _siteName()
BASE = REPO / SITE
if not (BASE / "site.toml").exists():
    sys.exit(f"{BASE / 'site.toml'} not found")
CFG = tomllib.loads((BASE / "site.toml").read_text())

NAME = CFG.get("title", SITE)                     # name used in titles
DATA = BASE / "data"
INPUTS = BASE / "Inputs"
RUNS = BASE / "runs"
OUTPUTS = BASE / "Outputs"
CRS = int(CFG["crs"])                             # EPSG code of the DEM; vector data is reprojected to it
SOURCE_DEM = BASE / CFG["sourceDem"]              # what prepInputs.py turns into DEM
RELEASE_NAME = CFG.get("releaseName", f"{SITE}_release")
RELEASE = INPUTS / f"REL/{RELEASE_NAME}.shp"
DEM = INPUTS / f"{SITE}_DEM.asc"
THICK = [float(t) for t in CFG.get("thicknesses", [0.25, 0.5, 0.75, 1.0, 1.5, 2.0])]
RHO = float(CFG.get("rho", 300))
FRICTION = CFG.get("friction", "samosATAuto")     # main series and animation
ROAD_FRICT = CFG.get("roadFrictions", ["samosATMedium", "samosAT"])  # fixed classes for the road sensitivity matrix
# optional shrub cover as a com1DFA resistance area: [shrubs] height (m) and the polygon file prepInputs.py writes
SHRUBS = CFG.get("shrubs")
RES = INPUTS / f"RES/{SITE}_shrubs.asc"
# extra named release scenarios ([[scenario]] in site.toml: id, release file, thickness, description), e.g. a D2 slab;
# prepInputs.py writes each to Inputs/scenarios/<id>/REL/<id>.shp
SCENARIOS = {s["id"]: s for s in CFG.get("scenario", [])}


def scenarioRel(sid):
    return INPUTS / "scenarios" / sid / "REL"


def dSize(massKg):
    """Destructive size from mass: typical masses D1 10 t, D2 100 t, D3 1,000 t, D4 10,000 t, D5 100,000 t (log scale),
    rounded to the nearest half size."""
    return min(5.0, max(1.0, round(2 * np.log10(massKg / 1000)) / 2))


def shrubOverrides():
    """com1DFA settings of the shrub variant: resistance simulation, shrubs overridden by flow deeper than their height."""
    return dict(simTypeList="res", forestThMax=f"{float(SHRUBS['height']):g}")


def readVector(path):
    """GeoDataFrame from geojson / shp / kml / kmz (a kmz is a zip holding a .kml)."""
    path = pathlib.Path(path)
    if path.suffix.lower() == ".kmz":
        with zipfile.ZipFile(path) as z:
            inner = next(n for n in z.namelist() if n.lower().endswith(".kml"))
        return gpd.read_file(f"/vsizip/{path}/{inner}")
    return gpd.read_file(path)


def readRaster(path):
    """Raster band 1 as float array (nodata -> nan), its transform and bounds."""
    with rasterio.open(path) as r:
        a = r.read(1).astype(float)
        if r.nodata is not None:
            a[a == r.nodata] = np.nan
        return a, r.transform, r.bounds


def extent(bounds):
    """imshow extent from rasterio bounds."""
    return bounds.left, bounds.right, bounds.bottom, bounds.top


def loadTerrain():
    """Hillshade of the input DEM (+ its extent) and the release polygon."""
    dem, T, B = readRaster(DEM)
    cs = abs(T.a)
    hill = LightSource(315, 45).hillshade(dem, vert_exag=1, dx=cs, dy=cs)
    return dict(dem=dem, T=T, hill=hill, ext=extent(B), rel=gpd.read_file(RELEASE))


def roadLines():
    """Road edge lines (hill side, far side) as 2D shapely lines in CRS, oriented the same way."""
    from shapely import force_2d
    from shapely.geometry import Point
    get = lambda key: force_2d(readVector(BASE / CFG["road"][key]).to_crs(CRS).geometry.union_all())
    hill, far = get("hillSide"), get("farSide")
    if hill.geom_type != "LineString" or far.geom_type != "LineString":
        sys.exit("road edges must each be one line (LineString)")
    # the road strip is hill + reversed far, so both must run in the same direction
    if Point(hill.coords[0]).distance(Point(far.coords[0])) > Point(hill.coords[0]).distance(Point(far.coords[-1])):
        far = far.reverse()
    return hill, far


def hasRoad():
    return "road" in CFG


def zoomWindow(masks, T, rel, pad):
    """x/y limits around all flow footprints (list of boolean arrays on transform T) and the release area."""
    rr, cc = np.nonzero(np.any(masks, axis=0))
    xs, ys = T * (cc + .5, rr + .5)
    x0, y0, x1, y1 = rel.total_bounds
    return (min(np.min(xs), x0) - pad, max(np.max(xs), x1) + pad), (min(np.min(ys), y0) - pad, max(np.max(ys), y1) + pad)


def peakRaster(avaDir, field):
    """Peak file (pft, pfv, ppr, ...) of the simulation in avaDir -> (array, transform, bounds)."""
    peakDir = pathlib.Path(avaDir) / "Outputs/com1DFA/peakFiles"
    return readRaster(next(peakDir.glob(f"*_dfa_{field}.asc")))


def runConfig(avaDir):
    """Config AvaFrame saved for the (single) simulation in avaDir, e.g. to see the friction class samosATAuto chose."""
    ini = next((pathlib.Path(avaDir) / "Outputs/com1DFA/configurationFiles").glob("*_dfa.ini"))
    c = configparser.ConfigParser()
    c.read(ini)
    return c


def releaseMass(avaDir):
    """Initial simulated mass [kg] of the simulation in avaDir."""
    f = next((pathlib.Path(avaDir) / "Outputs/com1DFA").glob("mass_*.txt"))
    return float(f.read_text().splitlines()[1].split(",")[1])


def runCom1DFA(avaDir, overrides=None, res=False, rel=None):
    """Run com1DFA in its own process for one avalanche dir under runs/.

    avaDir/Inputs gets links to the shared inputs and its own CFGs/com1DFACfg.ini = the site baseline
    (Inputs/CFGs/com1DFACfg.ini) with `overrides` ({key: value}, all in [GENERAL]) applied, so no global
    avaframe/com1DFA/local_com1DFACfg.ini is touched. Output goes to avaDir/Outputs, log to avaDir/run.log.
    Inputs/RES (resistance area) is linked only when res=True; AvaFrame would otherwise use it in every run.
    rel: a REL directory to use instead of Inputs/REL (named release scenarios).
    """
    avaDir = pathlib.Path(avaDir)
    inputs = avaDir / "Inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    for src in INPUTS.iterdir():
        if src.name in ("CFGs", "remeshedRasters", "scenarios") or src.suffix == ".json" or (src.name == "RES" and not res):
            continue
        if src.name == "REL" and rel is not None:
            src = pathlib.Path(rel)
        if not (inputs / src.name).exists():
            (inputs / src.name).symlink_to(src)
    for prj in (inputs / "remeshedRasters").glob("*.prj"):  # AvaFrame would glob the .prj sidecar as a raster
        prj.unlink()
    shutil.rmtree(avaDir / "Outputs", ignore_errors=True)

    cfg = configparser.ConfigParser(interpolation=None)
    cfg.optionxform = str
    cfg.read(INPUTS / "CFGs/com1DFACfg.ini")
    for k, v in (overrides or {}).items():
        cfg["GENERAL"][k] = str(v)
    (inputs / "CFGs").mkdir(exist_ok=True)
    with open(inputs / "CFGs/com1DFACfg.ini", "w") as fh:
        cfg.write(fh)

    log = avaDir / "run.log"
    with open(log, "w") as fh:
        res = subprocess.run([sys.executable, str(REPO / "avaframe/runCom1DFA.py"), str(avaDir)], cwd=REPO,
                             stdout=fh, stderr=subprocess.STDOUT)
    if res.returncode:
        sys.exit(f"com1DFA failed for {avaDir.name}, see {log}")
