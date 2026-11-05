"""Remove vegetation (tree / shrub crowns) from a surface model: progressive morphological filter + interpolation.

A grey-scale opening with a disk of diameter w removes raised objects narrower than w but leaves planar slopes unchanged,
so it works on steep terrain. Openings with growing windows (`windows`, m) flag a cell as vegetation when it sits more than
dh = dh0 + slope * w above the opened surface (the slope term keeps convex ground, e.g. gully banks and ridges, from being
cut). Flagged cells are filled by inverse-distance interpolation from the remaining ground cells, lightly smoothed, and
never set above the original surface. Same idea as Zhang et al. (2003) PMF / Pingel et al. (2013) SMRF, simplified.
"""
import numpy as np
from rasterio.fill import fillnodata
from scipy import ndimage


def disk(r):
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def groundFilter(z, cell, windows=(3, 5, 9, 15, 21), dh0=0.5, slope=0.15, maxFill=40):
    """z: 2D surface (nan = nodata), cell: cell size [m]. Returns (ground surface, vegetation mask)."""
    valid = ~np.isnan(z)
    # fill nodata with the nearest valid value so the filters do not see holes
    idx = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    zf = z[tuple(idx)]
    veg = np.zeros(z.shape, bool)
    for w in windows:
        r = max(1, int(round(w / cell / 2)))
        op = ndimage.grey_opening(zf, footprint=disk(r))
        veg |= (zf - op) > dh0 + slope * w
    veg &= valid
    g = np.where(veg | ~valid, np.nan, zf).astype("float32")
    mask = (~np.isnan(g)).astype("uint8")
    g = fillnodata(np.nan_to_num(g, nan=0), mask=mask, max_search_distance=maxFill / cell, smoothing_iterations=0)
    g = ndimage.uniform_filter(g, 3)                     # smooth the patched surface a little
    g = np.minimum(g, zf)                                # ground is never above the surface model
    return np.where(valid, g, np.nan), veg
