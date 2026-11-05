"""Build a site's full report in the right order: inputs, animation, road check, thickness series + HTML report.

Usage (from the repo root):
    uv run python siteTools/runSite.py --site McClure              # everything
    uv run python siteTools/runSite.py --site McClure --keep-rel   # prepInputs --dem-only (keep an adjusted Inputs/REL)
    uv run python siteTools/runSite.py --site McClure --skip-prep  # inputs already prepared
The HTML embeds the animation and road results only if they exist, so the series/report step runs last.
Output: <Site>/Outputs/thickness_series.html (+ .png, road_reach.{csv,png}, animation_<friction>/), and with [shrubs] in
site.toml also thickness_series_shrubs.{html,png}.
"""
import argparse
import subprocess
import sys

import avasite as b

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--keep-rel", action="store_true")
ap.add_argument("--skip-prep", action="store_true")
a = ap.parse_args()

here = b.REPO / "siteTools"
steps = [] if a.skip_prep else [["prepInputs.py"] + (["--dem-only"] if a.keep_rel else [])]
steps += [["makeAnimationSeries.py"]] + ([["roadReach.py"]] if b.hasRoad() else []) + [["runThicknessSeries.py"]]
if b.SHRUBS:  # shrub series, then redraw the plain report so its histogram includes the shrub runs
    steps += [["runThicknessSeries.py", "--shrubs"], ["runThicknessSeries.py", "--plot-only"]]
for s in steps:
    print(f"\n=== {b.SITE}: {' '.join(s)}", flush=True)
    subprocess.run([sys.executable, str(here / s[0]), "--site", b.SITE] + s[1:], check=True, cwd=b.REPO)
print(f"\nreport: {b.OUTPUTS / 'thickness_series.html'}")
