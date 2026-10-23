# Bethel avalanche simulations (AvaFrame com1DFA)

Dense-flow simulations for the Bethel path with release-thickness scenarios, interactive/static plots, an animation and a
road-reach check. Built on the AvaFrame repo this folder lives in.

## Layout

```
code/     scripts; bethel.py holds the shared paths/helpers
data/     raw sources: DEM (not tracked), release_area/*.geojson, Road/*.geojson, original release shapefile
Inputs/   AvaFrame inputs; CFGs/com1DFACfg.ini = site baseline settings (rho, friction, ...)
runs/     one AvaFrame avalanche dir per simulation (generated, not tracked)
Outputs/  figures, tables, animations (generated, not tracked)
```

## Get it running

```bash
git clone <AvaFrame repo> && cd AvaFrame
uv sync                       # installs AvaFrame (compiles its extensions) + plotly etc. into .venv
```

The 1 m DEM is not in git (about 260 MB). Put the USGS tile
`USGS_1M_13_x42y440_CO_DRCOG_2020_B20.tif` in `Bethel/data/DEM/`. Optional: `ffmpeg` on the PATH adds mp4/webm to the animation.

## Quick start (run from the repo root)

```bash
uv run python Bethel/code/prepInputs.py                 # crop DEM -> Inputs/Bethel_DEM.asc, release geojson -> shapefile
uv run python Bethel/code/adjustRelease.py              # optional: slope-adjusted release polygon
uv run python Bethel/code/runThicknessSeries.py         # runs 0.5/0.75/1/1.5/2 m, writes Outputs/thickness_series.{png,html}
uv run python Bethel/code/makeAnimationSeries.py        # side-by-side animation -> Outputs/animation_<friction>/
uv run python Bethel/code/roadReach.py                  # does the flow reach the road? -> Outputs/road_reach.{csv,png}
```

Useful options: `runThicknessSeries.py --thickness 1 2 --rho 250 --friction samosATMedium`, `--plot-only` / `--render-only`
to reuse existing runs, `--clean` to delete a series and rerun. Run one simulation directly with
`uv run python avaframe/runCom1DFA.py Bethel` (uses `Inputs/CFGs/com1DFACfg.ini`).

## Settings

`Inputs/CFGs/com1DFACfg.ini` is read by AvaFrame automatically for this avalanche dir and filled up from
`avaframe/com1DFA/com1DFACfg.ini`. Current baseline: `rho = 300`, `frictModel = samosATAuto` (AvaFrame picks the friction class
from the release volume; `samosATSmall` is too small for thickness >= 0.75 m). The series scripts write the per-run values
(thickness, friction, ...) to `runs/<series>/<run>/Inputs/CFGs/com1DFACfg.ini` on top of this baseline, so
`avaframe/com1DFA/local_com1DFACfg.ini` is not used or modified.

## Adapting to another site / path

1. Copy `Bethel/` to `<Site>/`, clear `data/`, `Inputs/` (keep `CFGs/`) and put your DEM, release polygon (and road lines) in `data/`.
2. `code/bethel.py`: set `CRS` (EPSG of your DEM) and the `RELEASE` / `DEM` file names (and `BASE` follows the folder automatically).
3. `code/prepInputs.py`: DEM crop box `cropBounds`, release file name/thickness. Output names must match `bethel.py`.
4. `Inputs/CFGs/com1DFACfg.ini`: density, friction model, defaults for your site.
5. Hard-coded to Bethel and to review: map windows in `roadReach.py` (`xl, yl`), `plotRelease.py` and `checkStripes.py`;
   slope threshold in `adjustRelease.py`; the "Bethel" titles in `runThicknessSeries.py` / `makeAnimationSeries.py`;
   `THICK`, `FRICT`, `RHO` in `roadReach.py`.
6. Run the quick start with the `<Site>/code/` paths.
