"""Map of the original vs slope-adjusted release polygon over the slope >= 28 deg terrain -> Outputs/release_adjust.png.

Run from /home/ron/AvaFrame after adjustRelease.py:  uv run python Bethel/code/plotRelease.py
"""
import pathlib, numpy as np, rasterio, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import bethel as b
base = b.BASE
# slope (degrees) of the 5 m block-mean DEM, same recipe as adjustRelease.py
with rasterio.open(base / "Inputs/Bethel_DEM.asc") as r: z = r.read(1).astype(float); T = r.transform
f = 5; n0, n1 = z.shape[0] // f * f, z.shape[1] // f * f
z5 = z[:n0, :n1].reshape(n0 // f, f, n1 // f, f).mean(axis=(1, 3)); gy, gx = np.gradient(z5, f)
slope = np.degrees(np.arctan(np.hypot(gx, gy)))
# map extent of the cropped 5 m grid
ext = (T.c, T.c + T.a * n1, T.f + T.e * n0, T.f)
# background: terrain >= 28 deg in grey (darker = steeper), < 28 deg in blue
fig, ax = plt.subplots(figsize=(8, 9))
im = ax.imshow(np.where(slope >= 28, slope, np.nan), extent=ext, cmap="Greys", vmin=20, vmax=55, alpha=.55)
ax.imshow(np.where(slope < 28, 1, np.nan), extent=ext, cmap="Blues", vmin=0, vmax=3, alpha=.5)
# release polygons: original (dashed orange) and adjusted (green); the map window is Bethel-specific
o = gpd.read_file(base / "data/Bethel_W_orig.shp"); n = gpd.read_file(base / "Inputs/REL/Bethel_W_Start_Zone.shp")
o.boundary.plot(ax=ax, color="#d95f02", lw=2, ls="--"); n.boundary.plot(ax=ax, color="#1b9e77", lw=2)
ax.set_xlim(424750, 425250); ax.set_ylim(4395300, 4395900)
# legend with proxy lines, because the slope shading has no artists of its own
ax.legend(handles=[Line2D([0], [0], color="#d95f02", ls="--", lw=2, label="original"), Line2D([0], [0], color="#1b9e77", lw=2, label="adjusted"),
                   Line2D([0], [0], color="#9ecae1", lw=6, label="slope < 28°"), Line2D([0], [0], color="#777", lw=6, label="slope >= 28°")], loc="upper right")
ax.set_title("Release area: original vs adjusted (5 m block-mean slope)"); ax.ticklabel_format(style="plain", useOffset=False)
b.OUTPUTS.mkdir(exist_ok=True)
fig.savefig(base / "Outputs/release_adjust.png", dpi=90)
