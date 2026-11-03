# Bethel avalanche simulations (AvaFrame com1DFA)

Dense-flow simulations for the Bethel path with release-thickness scenarios, a self-contained interactive HTML report
(linked maps, closest-approach histogram), static plots, an animation and a road-reach check. Built on the AvaFrame repo this folder lives in.

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

The DEM is not in git (about 540 MB). Put the survey ground model `261018_Bethel_Full_Ground-DSM.tiff` (4 cm cells,
EPSG:26913) in `Bethel/data/DEM/`; `code/bethel.py` (`SOURCE_DEM`) names it. It covers only the path corridor, so the
simulation domain ends at its edge. (The earlier USGS 1 m tile is no longer used, except by the `checkStripes.py` diagnostic.) Optional: `ffmpeg` on the PATH adds mp4/webm to the animation.

## Quick start (run from the repo root)

```bash
uv run python Bethel/code/prepInputs.py                 # DEM -> 1 m Inputs/Bethel_DEM.asc, release geojson -> shapefile
                                                        # (--dem-only keeps Inputs/REL, e.g. the adjusted polygon)
uv run python Bethel/code/adjustRelease.py              # optional: slope-adjusted release polygon
uv run python Bethel/code/makeAnimationSeries.py        # side-by-side animation -> Outputs/animation_<friction>/ (needs ffmpeg for the video)
uv run python Bethel/code/roadReach.py                  # does the flow reach the road? -> Outputs/road_reach.{csv,png}
uv run python Bethel/code/runThicknessSeries.py         # last step: runs 0.25/0.5/0.75/1/1.5/2 m, writes Outputs/thickness_series.{png,html}
```

`runThicknessSeries.py` writes the HTML last and embeds the animation (`flow.webm`/`flow.mp4`), the road sensitivity paragraph
(`road_reach.csv`) and the `roadReach.py` runs in the closest-approach histogram only if they already exist at that moment. So run it after `makeAnimationSeries.py` and `roadReach.py`,
or rerun it with `--plot-only` once those are done.

The report's maps (`code/interactiveMaps.py`) are plain canvas/SVG with the data embedded, no plotly or web tiles, so the
HTML is one file that works offline. One panel per thickness; zoom, pan and the crosshair are linked; click pins a point
(the URL gets `#at=E,N`, so the link opens at that spot). Distances are to where the traced avalanche path meets the
hill-side road edge. Printing shows the static PNG instead.

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
2. `code/bethel.py`: set `CRS` (EPSG of your DEM) and the `RELEASE` / `DEM` / `SOURCE_DEM` file names (and `BASE` follows the folder automatically).
3. `code/prepInputs.py`: optional DEM crop box `cropBounds`, cell size, release file name/thickness. Output names must match `bethel.py`.
4. `Inputs/CFGs/com1DFACfg.ini`: density, friction model, defaults for your site.
5. Hard-coded to Bethel and to review: map windows in `roadReach.py` (`xl, yl`), `plotRelease.py` and `checkStripes.py`;
   slope threshold in `adjustRelease.py`; the "Bethel" titles in `runThicknessSeries.py` / `makeAnimationSeries.py`;
   `THICK`, `FRICT`, `RHO` in `roadReach.py`.
6. Run the quick start with the `<Site>/code/` paths.
