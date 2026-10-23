"""Check raw 1 m DEM vs the 5 m prepared DEM for horizontal striping in the flow-path area."""
import pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import LightSource
from rasterio.windows import from_bounds
from scipy.ndimage import uniform_filter

import bethel as b
base = b.BASE
bb = (424700, 4394700, 425900, 4395900)  # flow-path area
# read a raster window for the box bb; returns (data, cell size, transform)
def win(path):
    with rasterio.open(path) as r:
        w = from_bounds(*bb, r.transform).round_offsets().round_lengths()
        return r.read(1, window=w).astype(float), abs(r.res[0]), r.transform
# raw 1 m tile vs the 5 m-prepared DEM Bethel_DEM.asc, same box
raw, c1, _ = win(next((base / "data/DEM").glob("*.tif")))
p5, c5, _ = win(base / "Inputs/Bethel_DEM.asc")

def stripe(a):
    """high-pass residual (value minus 9x9 local mean); compare row-mean vs column-mean of that residual"""
    hp = a - uniform_filter(a, 9)
    return hp.mean(axis=1), hp.mean(axis=0), hp
# numbers: std of the row-mean / column-mean residual [cm]; striping shows up as a larger row-mean (or column-mean) std
for name, a, c in (("raw 1 m", raw, c1), ("prepared 5 m", p5, c5)):
    rm, cm, _ = stripe(a)
    print(f"{name:14s} shape={a.shape}  std of row-mean residual={rm.std()*100:.2f} cm   column-mean residual={cm.std()*100:.2f} cm")

# figure per DEM: hillshade, high-pass residual, row/column mean profiles
fig, axs = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
for j, (name, a, c) in enumerate((("raw 1 m", raw, c1), ("prepared 5 m", p5, c5))):
    hs = LightSource(315, 45).hillshade(a, vert_exag=1, dx=c, dy=c)
    axs[j, 0].imshow(hs, cmap="gray"); axs[j, 0].set_title(f"Hillshade - {name}")
    rm, cm, hp = stripe(a)
    im = axs[j, 1].imshow(hp, cmap="RdBu_r", vmin=-.3, vmax=.3); axs[j, 1].set_title("High-pass residual (m)")
    fig.colorbar(im, ax=axs[j, 1], shrink=.7)
    axs[j, 2].plot(rm * 100, np.arange(len(rm)) * c, lw=.8, label="row means"); axs[j, 2].invert_yaxis()
    axs[j, 2].plot(cm * 100, np.arange(len(cm)) * c, lw=.8, label="column means (plotted vs. x idx)", alpha=.6)
    axs[j, 2].set_xlabel("mean residual (cm)"); axs[j, 2].set_ylabel("m along axis"); axs[j, 2].legend(fontsize=8)
    axs[j, 2].set_title("Row/column mean of residual")
b.OUTPUTS.mkdir(exist_ok=True)
fig.savefig(base / "Outputs/dem_stripe_check.png", dpi=100)
