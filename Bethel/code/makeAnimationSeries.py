"""Run com1DFA with saved time steps for several release thicknesses and animate them side by side.

Layout: flow thickness only, one panel per release thickness in a grid of --ncols columns (default 2, so 6 runs = 3 rows),
one shared colour scale and a shared clock (a run that has stopped holds its last frame). Each thickness runs in its own avalanche dir
Bethel/runs/anim_<friction>/th_<t>/.

Usage (from /home/ron/AvaFrame):
    uv run python Bethel/code/makeAnimationSeries.py                          # run all + render (samosATAuto, rho 300)
    uv run python Bethel/code/makeAnimationSeries.py --render-only            # reuse existing runs
    uv run python Bethel/code/makeAnimationSeries.py --friction samosATSmall --thickness 1 2 --dt 1
Writes frames (PNG), flow.gif and, if ffmpeg is installed, flow.mp4 and flow.webm to Bethel/Outputs/animation_<friction>/.
"""
import argparse
import re
import shutil
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

import bethel as b

thr = 0.05  # m, flow thickness shown above this


def runOne(avaDir, relTh, rho, friction, dt):
    """com1DFA with flow fields saved every dt seconds, so the animation has time steps."""
    b.runCom1DFA(avaDir, dict(rho=rho, frictModel=friction, relThFromFile=False, relTh=relTh,
                              resType="ppr|pft|pfv|FT|FV", tSteps=f"0:{dt}"))
    print(f"relTh {relTh} m done")


def loadRun(avaDir):
    """Time-step rasters of one run: {t: (FT, FV)}, plus the friction model AvaFrame actually used."""
    steps = avaDir / "Outputs/com1DFA/peakFiles/timeSteps"
    data = {}
    for p in steps.glob("*_dfa_FT_t*.asc"):
        t = float(re.search(r"_t([\d.]+)\.asc$", p.name).group(1))
        ft, T, B = b.readRaster(p)
        fv = b.readRaster(steps / p.name.replace("_FT_", "_FV_"))[0]
        data[t] = (ft, fv)
    return data, b.runConfig(avaDir)["GENERAL"]["frictModel"], T, B


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--thickness", type=float, nargs="+", default=[0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
    ap.add_argument("--rho", type=float, default=300)
    ap.add_argument("--friction", default="samosATAuto")
    ap.add_argument("--dt", type=float, default=2, help="seconds between saved frames")
    ap.add_argument("--fps", type=int, default=6)
    ap.add_argument("--ncols", type=int, default=2, help="panels per row")
    ap.add_argument("--render-only", action="store_true")
    ap.add_argument("--skip-existing", action="store_true", help="only run thicknesses without saved time steps")
    a = ap.parse_args()
    # one avalanche dir per thickness: runs/anim_<friction>/th_<t>
    root = b.RUNS / f"anim_{a.friction}"
    outDir = b.OUTPUTS / f"animation_{a.friction}"
    dirs = {t: root / f"th_{t:g}" for t in a.thickness}

    # run the simulations (skippable), then render from whatever is in runs/
    if not a.render_only:
        for t, d in dirs.items():
            if a.skip_existing and any((d / "Outputs/com1DFA/peakFiles/timeSteps").glob("*.asc")):
                print(f"relTh {t:g} m: existing run kept")
                continue
            runOne(d, f"{t:g}", f"{a.rho:g}", a.friction, f"{a.dt:g}")

    # load all runs: {thickness: (time steps, friction class used, transform, bounds)}; all runs share one grid
    runs = {t: loadRun(d) for t, d in dirs.items()}
    T, B = runs[a.thickness[0]][2:]
    ext = b.extent(B)
    for t in a.thickness:
        print(f"relTh {t:g} m: friction used = {runs[t][1]}, {len(runs[t][0])} steps, ends at {max(runs[t][0]):.1f} s")
    # shared clock = union of all saved times; colour scales (max over all runs/steps) are fixed so frames are comparable
    clock = sorted(set(x for r in runs.values() for x in r[0]))
    vmaxT = max(np.nanmax(ft) for r in runs.values() for ft, _ in r[0].values())

    # fixed map window around everything the flow ever covers + the release area
    terrain = b.loadTerrain()
    rel, hill, dExt = terrain["rel"], terrain["hill"], terrain["ext"]
    anyF = [ft > thr for r in runs.values() for ft, _ in r[0].values()]
    xl, yl = b.zoomWindow(anyF, T, rel, pad=100)
    anyF = np.any(anyF, axis=0)

    # figure: one flow-thickness panel per release thickness, filled row by row. The artists are created once
    # and only their data is swapped per frame, which is much faster than redrawing.
    outDir.mkdir(parents=True, exist_ok=True)
    for old in outDir.glob("frame_*.png"):
        old.unlink()
    n = len(a.thickness)
    ncol = min(a.ncols, n); nrow = -(-n // ncol)
    pw = 5.0                                                     # panel width [in]; height follows the map window
    ph = pw * (yl[1] - yl[0]) / (xl[1] - xl[0])
    fig, axs = plt.subplots(nrow, ncol, figsize=(pw * ncol + 1.2, ph * nrow + .6), constrained_layout=True, squeeze=False)
    ims = {}
    for j, t in enumerate(a.thickness):
        ax = axs.flat[j]
        ax.imshow(hill, cmap="gray", extent=dExt, alpha=.8)
        ims[j] = ax.imshow(np.full(anyF.shape, np.nan), cmap="Blues", vmin=0, vmax=vmaxT, extent=ext)
        rel.boundary.plot(ax=ax, color="#d95f02", lw=1.2)
        ax.set_xlim(xl); ax.set_ylim(yl); ax.set_aspect("equal"); ax.tick_params(labelsize=6)
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.set_title(f"Release {t:g} m ({runs[t][1]})", loc="left", fontsize=10)
    for ax in axs.flat[n:]:                                      # unused grid cells
        ax.set_visible(False)
    fig.colorbar(ims[0], ax=axs, shrink=.5, pad=.01, label="flow thickness (m)")
    # render one PNG per time step
    title = fig.suptitle("")
    frames = []
    for k, tc in enumerate(clock):
        for j, t in enumerate(a.thickness):
            data = runs[t][0]
            ts = max(x for x in data if x <= tc)  # hold last saved step (and final state once stopped)
            ft = data[ts][0]
            ims[j].set_data(np.where(ft > thr, ft, np.nan))
        title.set_text(f"Bethel avalanche: flow thickness, density {a.rho:g} kg/m3, t = {tc:5.1f} s")
        f = outDir / f"frame_{k:03d}_t{tc:06.1f}s.png"
        fig.savefig(f, dpi=70); frames.append(f)
    # assemble the GIF (last frame held longer)
    imgs = [Image.open(f).convert("P", palette=Image.ADAPTIVE) for f in frames]
    dur = [int(1000 / a.fps)] * len(imgs); dur[-1] = 2500
    imgs[0].save(outDir / "flow.gif", save_all=True, append_images=imgs[1:], duration=dur, loop=0, optimize=True)
    print(f"{len(frames)} frames -> {outDir}\ngif: {outDir / 'flow.gif'} ({(outDir / 'flow.gif').stat().st_size / 1e6:.1f} MB)")
    # optional video encodes: mp4 (H.264) and webm (VP9, plays in browsers without H.264)
    if shutil.which("ffmpeg"):
        subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(a.fps), "-pattern_type", "glob", "-i",
                        str(outDir / "frame_*.png"), "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p",
                        str(outDir / "flow.mp4")], check=True)
        print("mp4:", outDir / "flow.mp4")
        subprocess.run(["ffmpeg", "-loglevel", "error", "-y", "-framerate", str(a.fps), "-pattern_type", "glob", "-i",
                        str(outDir / "frame_*.png"), "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-c:v", "libvpx-vp9", "-crf", "34", "-b:v", "0",
                        str(outDir / "flow.webm")], check=False)  # webm plays in browsers without H.264
        print("webm:", outDir / "flow.webm")


if __name__ == "__main__":
    main()
