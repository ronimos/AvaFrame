"""Shared helpers for the Bethel scripts: paths, raster/terrain loading, and running com1DFA per avalanche dir.

Layout (everything relative to Bethel/):
    code/     these scripts
    data/     raw source data (DEM tif, release geojson, road lines)
    Inputs/   AvaFrame inputs (DEM, REL, CFGs/com1DFACfg.ini = baseline settings for this site)
    runs/     one AvaFrame avalanche dir per simulation (generated)
    Outputs/  figures, tables, animations (generated)
"""
import configparser
import pathlib
import shutil
import subprocess
import sys

import geopandas as gpd
import numpy as np
import rasterio
from matplotlib.colors import LightSource

BASE = pathlib.Path(__file__).resolve().parent.parent
REPO = BASE.parent
DATA = BASE / "data"
INPUTS = BASE / "Inputs"
RUNS = BASE / "runs"
OUTPUTS = BASE / "Outputs"
RELEASE = INPUTS / "REL/Bethel_W_Start_Zone.shp"
DEM = INPUTS / "Bethel_DEM.asc"
# source DEM that prepInputs.py turns into DEM (survey ground model, 4 cm cells; the older USGS 1 m tile is no longer used)
SOURCE_DEM = DATA / "DEM/261018_Bethel_Full_Ground-DSM.tiff"
# EPSG code of the DEM / release / road data (UTM 13N). Change for another site.
CRS = 26913


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
    """Road edge lines (hill side, far side) from data/Road as 2D shapely lines in CRS."""
    from shapely import force_2d
    get = lambda pat: force_2d(gpd.read_file(next((DATA / "Road").glob(pat))).to_crs(CRS).geometry.iloc[0])
    return get("*hill*.geojson"), get("*far*.geojson")


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


def runCom1DFA(avaDir, overrides=None):
    """Run com1DFA in its own process for one avalanche dir under runs/.

    avaDir/Inputs gets links to the shared inputs and its own CFGs/com1DFACfg.ini = the site baseline
    (Inputs/CFGs/com1DFACfg.ini) with `overrides` ({key: value}, all in [GENERAL]) applied, so no global
    avaframe/com1DFA/local_com1DFACfg.ini is touched. Output goes to avaDir/Outputs, log to avaDir/run.log.
    """
    avaDir = pathlib.Path(avaDir)
    inputs = avaDir / "Inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    for src in INPUTS.iterdir():
        if src.name in ("CFGs", "remeshedRasters"):
            continue
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
