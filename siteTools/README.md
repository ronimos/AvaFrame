# siteTools: avalanche path runout scenarios with AvaFrame com1DFA

One set of scripts for every site. A site is a folder in the repo root (`McClure/`, `Bethel/`) with a `site.toml`
(see `McClure/site.toml` for all keys) and its raw data in `data/`. Pick the site with `--site <Name>` (or `AVA_SITE`).

## New site

1. Make `<Site>/data/` with the DEM (GeoTIFF), the release polygon and the two road edge lines (geojson, shp, kml or kmz).
2. Write `<Site>/site.toml`: `crs`, `sourceDem`, `release`, `[road] hillSide/farSide`; optional thicknesses, density,
   friction, `[vegetationFilter]`, `[shrubs]`, `[[scenario]]` (e.g. a D2 slab) and `[report]` text.
3. Copy `McClure/.gitignore` and rename the DEM line.
4. `uv run python siteTools/runSite.py --site <Site>` (inputs, animation, road matrix, series, shrub series, report).
5. Check `<Site>/Outputs/thickness_series*.html`, then `uv run python siteTools/publish.py --site <Site>` copies the
   reports to /home/www/html/ron/avalanche_paths/ as `<Site>[_shrubs]_runout_scenarios.html` and rebuilds the index.

## Scripts

| script | does |
|---|---|
| `avasite.py` | site settings and paths, raster/vector helpers, `runCom1DFA` (one avalanche dir per run under `runs/`) |
| `prepInputs.py` | DEM -> 1 m `.asc` (+ vegetation filter), release and scenario polygons -> shapefiles, shrub resistance raster |
| `groundFilter.py` | removes tree/shrub crowns from a surface model (progressive morphological filter) |
| `makeAnimationSeries.py` | flow-thickness animation, one panel per thickness |
| `roadReach.py` | thickness x friction class matrix: does the flow reach / cross the road |
| `runThicknessSeries.py` | the main series (+ `--shrubs`), static PNG and the self-contained HTML report |
| `interactiveMaps.py` | linked canvas maps and closest-approach histogram used in the report |
| `runSite.py` / `publish.py` | run everything in order / copy reports to the web folder |

Report sizes (D1-D5) are estimated from release mass only (typical masses 10 t ... 100,000 t).
`Bethel/code/` still holds the older Bethel-only copies of these scripts.
