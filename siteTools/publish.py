"""Copy a site's report to the shared web folder and rebuild its index page.

Usage (from the repo root):  uv run python siteTools/publish.py --site McClure
Reports: <Site>/Outputs/thickness_series[_<variant>].html -> WEB/<Site>[_<variant>]_runout_scenarios.html (e.g. the
_shrubs variant). site.toml [report] publishVariants (e.g. [""] = plain only) limits which variants go online; this site's
other variants are then removed from WEB. WEB/index.html lists every report there,
with the site title and date read from each report. WEB defaults to /home/www/html/ron/avalanche_paths (env AVA_WEB overrides).
"""
import datetime
import html
import os
import pathlib
import re
import shutil

import avasite as b

WEB = pathlib.Path(os.environ.get("AVA_WEB", "/home/www/html/ron/avalanche_paths"))
WEB.mkdir(parents=True, exist_ok=True)
keep = b.CFG.get("report", {}).get("publishVariants")            # None = every variant
for src in sorted(b.OUTPUTS.glob("thickness_series*.html")):
    variant = src.stem.removeprefix("thickness_series")          # "" or "_shrubs"
    if keep is not None and variant not in keep:
        old = WEB / f"{b.SITE}{variant}_runout_scenarios.html"
        if old.exists():
            old.unlink()
            print("removed", old, "(not in publishVariants)")
        continue
    dst = WEB / f"{b.SITE}{variant}_runout_scenarios.html"
    shutil.copy2(src, dst)
    print("copied", src, "->", dst)

rows = []
for f in sorted(WEB.glob("*_runout_scenarios.html")):
    head = f.read_text(errors="ignore")[:4000]
    title = re.search(r"<title>(.*?)</title>", head, re.S)
    gen = re.search(r"Generated (\d{4}-\d{2}-\d{2})", f.read_text(errors="ignore")[-3000:])
    rows.append(f'<tr><td><a href="{f.name}">{f.name.removesuffix("_runout_scenarios.html")}</a></td>'
                f'<td>{title.group(1).strip() if title else ""}</td><td>{gen.group(1) if gen else ""}</td></tr>')
(WEB / "index.html").write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Avalanche path runouts</title>
<style>body {{ font-family: system-ui, sans-serif; max-width: 900px; margin: 0 auto; padding: 0 16px 40px; color: #222; }}
table {{ border-collapse: collapse; width: 100%; }} td, th {{ border-bottom: 1px solid #ddd; padding: 6px 8px; text-align: left; }}
.sub {{ color: #555; }}</style></head><body>
<h1>Avalanche path runout scenarios</h1>
<p class="sub">AvaFrame com1DFA simulations by release thickness, one interactive report per path. Uncalibrated model runs;
 read each report's methodology and limitations before using the numbers.</p>
<table><tr><th>Path</th><th>Report</th><th>Generated</th></tr>{''.join(rows)}</table>
<p class="sub">Index updated {datetime.date.today().isoformat()}.</p>
</body></html>""")
print("index:", WEB / "index.html", f"({len(rows)} reports)")
