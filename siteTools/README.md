# siteTools: avalanche path runout scenarios with AvaFrame com1DFA

One set of scripts for every site. A site is a folder in the repo root (`McClure/`, `Bethel/`) with a `site.toml`
(see `McClure/site.toml` for all keys) and its raw data in `data/`. Pick the site with `--site <Name>` (or `AVA_SITE`).

## New site

Drop the files into `<Name>/data/DEM/` (one GeoTIFF), `<Name>/data/release_area/` (the start zone polygon) and
`<Name>/data/Road/` (the two road edge lines; kml, kmz, geojson or shp), then:

```bash
uv run python siteTools/newSite.py <Name> --title "Path name" [--shrubs 1.5] [--d2 d2_outline.kml]
```

or give the files directly (`--dem`, `--release`, `--road EDGE1 EDGE2`); they are copied in. `newSite.py` checks the
inputs (CRS, coverage of start zone / road / runout, which road edge is the hill side, vegetation in the surface model,
travel angle, size by mass), writes `site.toml`, turns on the vegetation filter when the DEM still holds crowns, adds a D2
scenario (yours with `--d2`, otherwise an assumed ~100 t slab at the top of the start zone), and runs the pipeline.
A FAIL stops it; WARNs are printed and saved in `<Name>/checks.json`. Then review the draft
`<Name>/Outputs/thickness_series*.html` and publish:

```bash
uv run python siteTools/publish.py --site <Name>     # -> /home/www/html/ron/avalanche_paths/
```

To change settings later, edit `<Name>/site.toml` and rerun `uv run python siteTools/runSite.py --site <Name>`.

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
| `newSite.py` | set up and check a new site from its files, then run it |
| `runSite.py` / `publish.py` | run everything in order / copy reports to the web folder |

Report sizes (D1-D5) are estimated from release mass only (typical masses 10 t ... 100,000 t).
`Bethel/code/` still holds the older Bethel-only copies of these scripts.
