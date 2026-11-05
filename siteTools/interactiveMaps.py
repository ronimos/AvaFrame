"""Interactive linked maps + closest-approach histogram for the self-contained site HTML report.

Python builds one JSON payload; plain JavaScript draws it with <canvas> and inline SVG (no server, no tiles, no JS
libraries), so the report stays one file that can be emailed. Adapted from the avachain runout-forecast report:
there one panel per forecast day shows a probability map; here one panel per release-thickness run shows a peak
field (thickness / velocity / pressure, switchable).

Pieces:
  * linked maps: shared zoom/pan/crosshair, click to pin (URL #at=E,N), exact cursor readout from gzipped uint16 fields
  * avalanche path: D8 steepest descent on the pit-filled DEM from the top of the flow core, stopped 40 m from the
    road and projected onto it; the hit point ("road crossing") is the reference for all distances
  * histogram: closest approach of every run to that crossing, 25 m bins, one bar per friction setting
  * print: the static matplotlib figure replaces the interactive block

Entry point: section(...) returns the HTML (CSS + markup + payload + JS) to paste into the report body.
"""
import base64
import gzip
import heapq
import io
import json

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import LightSource
from PIL import Image
from rasterio.features import rasterize
from rasterio.warp import Resampling, reproject
from shapely import distance, points
from shapely.geometry import LineString, Point

STOP_M = 40        # stop following the terrain this far from the road (nearer, D8 turns along the road-cut ditch)
DIR_M = 100        # track direction is taken over the last DIR_M of the traced path
PROJECT_M = 300    # ... and projected at most this far to meet the road
GROW_MAX_M = 150   # the view box may grow at most this much to include nearby road
BIN_M = 25         # histogram bin width


# ---------------------------------------------------------------- encoding helpers
def _gz(a):
    """Array -> gzipped bytes -> base64 text (browser: DecompressionStream('gzip'))."""
    return base64.b64encode(gzip.compress(np.ascontiguousarray(a).tobytes(), 9)).decode()


def _png(rgba):
    buf = io.BytesIO()
    Image.fromarray(rgba).save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _viridis(v, vmax, show):
    """RGBA viridis image of v/vmax, alpha 0.8, transparent where not `show`."""
    rgba = (colormaps["viridis"](np.clip(np.nan_to_num(v) / vmax, 0, 1)) * 255).astype("uint8")
    rgba[..., 3] = np.where(show, 204, 0)
    return rgba


# ---------------------------------------------------------------- avalanche path and road crossing
def priorityFlood(z):
    """Pit-filled copy of z (nan = outside): priority flood + epsilon (Barnes et al. 2014), so every cell drains.

    Seeds are the valid cells on the grid edge or next to nan; cells are raised to just above their lowest
    already-processed neighbour, so D8 descent on the result never gets stuck.
    """
    h, w = z.shape
    valid = ~np.isnan(z)
    out = z.copy()
    done = ~valid
    pq = []
    pad = np.pad(valid, 1, constant_values=False)
    edge = valid & ~(pad[:-2, 1:-1] & pad[2:, 1:-1] & pad[1:-1, :-2] & pad[1:-1, 2:]
                     & pad[:-2, :-2] & pad[:-2, 2:] & pad[2:, :-2] & pad[2:, 2:])
    for r, c in zip(*np.nonzero(edge)):
        heapq.heappush(pq, (out[r, c], r, c)); done[r, c] = True
    while pq:
        zc, r, c = heapq.heappop(pq)
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                rr, cc = r + dr, c + dc
                if (dr or dc) and 0 <= rr < h and 0 <= cc < w and not done[rr, cc]:
                    done[rr, cc] = True
                    out[rr, cc] = max(out[rr, cc], np.nextafter(zc, np.inf) + 1e-4)
                    heapq.heappush(pq, (out[rr, cc], rr, cc))
    return out


def avalanchePath(dem, T, core, road):
    """Trace the avalanche path and find where it meets the road line.

    dem, core: arrays on transform T (core = cells reached by at least half the runs). Start = highest core cell,
    follow D8 steepest descent on the filled DEM until STOP_M from `road`, then project the direction of the
    last DIR_M up to PROJECT_M onto the road. Returns dict(line=[[x, y], ...], crossing=[x, y], rule=...).
    """
    filled = priorityFlood(dem)
    res = abs(T.a)
    zc = np.where(core, dem, -np.inf)
    r, c = np.unravel_index(np.argmax(zc), zc.shape)
    xy = lambda r, c: tuple(T * (c + .5, r + .5))
    line = [xy(r, c)]
    for _ in range(dem.size):
        if road.distance(Point(line[-1])) <= STOP_M:
            break
        best, bd = None, 0.0
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                rr, cc = r + dr, c + dc
                if (dr or dc) and 0 <= rr < dem.shape[0] and 0 <= cc < dem.shape[1] and not np.isnan(filled[rr, cc]):
                    s = (filled[r, c] - filled[rr, cc]) / np.hypot(dr, dc)
                    if s > bd:
                        best, bd = (rr, cc), s
        if best is None:
            break
        r, c = best
        line.append(xy(r, c))
    pts = np.array(line)
    end = pts[-1]
    # direction over the last DIR_M of path length
    seg = np.r_[0, np.cumsum(np.hypot(*np.diff(pts, axis=0).T))]
    back = pts[np.searchsorted(seg, seg[-1] - DIR_M, side="right") - 1] if seg[-1] > 0 else end
    v = end - back
    rule = "nearest road point to the path end"
    crossing = None
    if np.hypot(*v) > 0:
        hit = LineString([end, end + v / np.hypot(*v) * PROJECT_M]).intersection(road)
        if not hit.is_empty:
            # first hit along the projection
            cand = [hit] if hit.geom_type == "Point" else list(getattr(hit, "geoms", [hit]))
            p = min((Point(q.coords[0]) for q in cand), key=lambda q: q.distance(Point(end)))
            crossing, rule = (p.x, p.y), f"path direction over the last {DIR_M} m projected onto the road"
            line.append(crossing)
    if crossing is None:
        p = road.interpolate(road.project(Point(end)))
        crossing = (p.x, p.y)
    return dict(line=[[round(x, 1), round(y, 1)] for x, y in line], crossing=[round(crossing[0], 1), round(crossing[1], 1)],
                rule=rule, res=res)


def closestApproach(flow, T, crossing, roadMask):
    """Min distance [m] from any flow cell centre to the crossing; 0 if the flow touches the road anywhere."""
    if (flow & roadMask).any():
        return 0.0
    rr, cc = np.nonzero(flow)
    if not len(rr):
        return None
    x, y = T * (cc + .5, rr + .5)
    return float(np.min(np.hypot(x - crossing[0], y - crossing[1])))


# ---------------------------------------------------------------- payload
def payload(panels, T, fields, thr, roads, roadRef, release, dem1, demT, groups):
    """Build the JSON payload.

    panels: [{id, title, sub, fields: {key: array}}] all on transform T (same shape, checked)
    fields: [{key, name, unit, scale, dec}] (scale = factor to store the value as uint16)
    roads: list of shapely lines to draw; roadRef: the line that counts as "the road" (crossing / reached)
    release: list of polygons; dem1, demT: hillshade DEM and its transform
    groups: [{name, runs: [{id, pft}]}] for the histogram, pft on transform T
    """
    shape = panels[0]["fields"]["pft"].shape
    for p in panels + [r for g in groups for r in g["runs"]]:
        arrs = p["fields"].values() if "fields" in p else [p["pft"]]
        if any(a.shape != shape for a in arrs):
            raise ValueError(f"run {p['id']} is not on the shared grid {shape}; the linked crosshair needs one grid")
    res = abs(T.a)
    flows = [p["fields"]["pft"] > thr for p in panels]
    count = np.sum(flows, axis=0).astype("uint8")

    # view extent: footprint + release, padded, grown to include the road nearby (by at most GROW_MAX_M)
    rr, cc = np.nonzero(count)
    xs, ys = T * (cc + .5, rr + .5)
    rb = np.array([g.bounds for g in release])
    x0, x1 = min(xs.min(), rb[:, 0].min()) - 60, max(xs.max(), rb[:, 2].max()) + 60
    y0, y1 = min(ys.min(), rb[:, 1].min()) - 60, max(ys.max(), rb[:, 3].max()) + 60
    lb = roadRef.bounds
    x0, y0 = max(min(x0, lb[0]), x0 - GROW_MAX_M), max(min(y0, lb[1]), y0 - GROW_MAX_M)
    x1, y1 = min(max(x1, lb[2]), x1 + GROW_MAX_M), min(max(y1, lb[3]), y1 + GROW_MAX_M)
    # crop the simulation grid to the extent (whole cells)
    c0 = max(int((x0 - T.c) // res), 0); c1 = min(int(np.ceil((x1 - T.c) / res)), shape[1])
    r0 = max(int((T.f - y1) // res), 0); r1 = min(int(np.ceil((T.f - y0) / res)), shape[0])
    crop = lambda a: a[r0:r1, c0:c1]
    grid = dict(left=T.c + c0 * res, top=T.f - r0 * res, res=res, w=c1 - c0, h=r1 - r0)
    x0, x1 = max(x0, grid["left"]), min(x1, grid["left"] + grid["w"] * res)
    y1, y0 = min(y1, grid["top"]), max(y0, grid["top"] - grid["h"] * res)

    # path and crossing from the 5 m DEM on the simulation grid (DEM block-averaged onto T)
    dem5 = np.full(shape, np.nan, "float32")
    reproject(dem1.astype("float32"), dem5, src_transform=demT, dst_transform=T, src_crs="EPSG:26913",
              dst_crs="EPSG:26913", src_nodata=np.nan, dst_nodata=np.nan, resampling=Resampling.average)
    core = count >= max(1, int(np.ceil(len(panels) / 2)))
    path = avalanchePath(dem5.astype(float), T, core, roadRef)
    roadMask = rasterize([roadRef], out_shape=shape, transform=T, all_touched=True).astype(bool)

    # fields: shared colour range per field over all panels (footprint cells only)
    fmeta = []
    for f in fields:
        vmax = max(float(np.nanmax(np.where(fl, p["fields"][f["key"]], np.nan))) for p, fl in zip(panels, flows))
        fmeta.append(dict(key=f["key"], name=f["name"], unit=f["unit"], scale=f["scale"], dec=f["dec"], vmax=vmax))
    outPanels = []
    for p, fl in zip(panels, flows):
        d = closestApproach(fl, T, path["crossing"], roadMask)
        status = p.get("status", "") + ("" if not d else f"; closest {d:.0f} m from the road crossing")
        q = dict(id=p["id"], title=p["title"], sub=p["sub"], status=status, d=d, data={}, png={})
        for f in fmeta:
            v = crop(np.where(fl, np.nan_to_num(p["fields"][f["key"]]), 0))
            q["data"][f["key"]] = _gz(np.clip(np.round(v * f["scale"]), 0, 65535).astype("<u2"))
            q["png"][f["key"]] = _png(_viridis(v, f["vmax"], crop(fl)))
        outPanels.append(q)

    # hillshade (light 315/45, z x 1.3) of the 1 m DEM inside the extent, transparent outside the DEM
    hc0, hr0 = int((x0 - demT.c) // demT.a), int((demT.f - y1) // -demT.e)
    hc1, hr1 = int(np.ceil((x1 - demT.c) / demT.a)), int(np.ceil((demT.f - y0) / -demT.e))
    hc0, hr0 = max(hc0, 0), max(hr0, 0)
    sub = dem1[hr0:hr1, hc0:hc1]
    ok = ~np.isnan(sub)
    hs = LightSource(315, 45).hillshade(np.where(ok, sub, np.nanmean(sub)), vert_exag=1.3, dx=demT.a, dy=-demT.e)
    hill = np.dstack([(hs * 255).astype("uint8")] * 3 + [np.where(ok, 255, 0).astype("uint8")])

    hist = [dict(name=g["name"], scen=[dict(id=r["id"], d=closestApproach(r["pft"] > thr, T, path["crossing"], roadMask))
                                       for r in g["runs"]]) for g in groups]
    xyList = lambda g: [[round(x, 1), round(y, 1)] for x, y, *_ in g.coords]
    return dict(
        grid=grid, extent=dict(x0=x0, x1=x1, y0=y0, y1=y1), thr=thr, fields=fmeta, panels=outPanels,
        count=_gz(crop(count)), countN=len(panels), bin=BIN_M,
        hill=dict(png=_png(hill), left=demT.c + hc0 * demT.a, top=demT.f + hr0 * demT.e,
                  w=sub.shape[1] * demT.a, h=sub.shape[0] * -demT.e),
        roads=[xyList(g) for g in roads], release=[xyList(g.exterior) for g in release], path=path, hist=hist,
        ramp=[f"#{''.join(f'{int(c * 255):02x}' for c in colormaps['viridis'](i / 8)[:3])}" for i in range(9)])


# ---------------------------------------------------------------- HTML / CSS / JS
CSS = """
.imap { margin: 8px 0 18px; }
.imap-bar { display: flex; flex-wrap: wrap; gap: 6px; align-items: center; margin: 4px 0 8px; }
.imap-bar button { font: inherit; padding: 3px 10px; border: 1px solid #999; background: #fff; border-radius: 3px; cursor: pointer; }
.imap-bar button.on { background: #222; color: #fff; border-color: #222; }
.imap-read { border: 1px solid #ccc; padding: 6px 10px; margin: 0 0 8px; font-size: .9rem; font-variant-numeric: tabular-nums; overflow-x: auto; }
.imap-read table { margin: 0 0 4px; border-collapse: collapse; } .imap-read td, .imap-read th { border: none; border-bottom: 1px solid #eee; padding: 1px 12px 1px 0; text-align: right; white-space: nowrap; }
.imap-read td:first-child, .imap-read th:first-child, .imap-read td.l, .imap-read th.l { text-align: left; }
/* equal-width panels that wrap on narrow screens (a flex row would stretch a lone last panel) */
.imap-row { display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 10px; }
.imap-panel { min-width: 0; }
.imap-title { font-size: .85rem; margin: 0 0 3px; line-height: 1.3; min-height: 3.9em; }
.imap-box { position: relative; border: 1px solid #888; background: #eee; overflow: hidden; touch-action: none; }
.imap-box canvas { position: absolute; inset: 0; width: 100%; height: 100%; cursor: crosshair; }
.imap-tag { position: absolute; left: 4px; top: 4px; background: rgba(255,255,255,.85); padding: 0 5px; font-size: .8rem;
  font-variant-numeric: tabular-nums; pointer-events: none; border-radius: 2px; }
.imap-legend { display: flex; flex-wrap: wrap; gap: 6px 22px; align-items: center; font-size: .8rem; margin: 8px 0 2px; }
.imap-ramp { display: inline-block; width: 180px; height: 10px; vertical-align: middle; border: 1px solid #888; }
.imap-help { color: #555; font-size: .8rem; }
.imap-hist { position: relative; margin: 10px 0; max-width: 900px; }
.imap-hist svg { width: 100%; height: auto; font-size: 11px; }
.imap-tip { position: absolute; display: none; background: #fff; border: 1px solid #888; padding: 4px 8px; font-size: .8rem;
  pointer-events: none; box-shadow: 0 1px 4px rgba(0,0,0,.2); white-space: nowrap; font-variant-numeric: tabular-nums; }
.imap-print { display: none; }
@media print { .imap, .imap-hist-wrap { display: none; } .imap-print { display: block; max-width: 100%; } }
"""

MARKUP = """
<div class="imap" id="imap">
 <div class="imap-bar" id="imap-fields"><span>Show:</span></div>
 <div class="imap-read" id="imap-read"></div>
 <div class="imap-row" id="imap-row"></div>
 <div class="imap-legend" id="imap-legend"></div>
</div>
<img class="imap-print" src="{printPng}" alt="Static maps of the same runs">
"""

HIST = """
<div class="imap-hist-wrap">
<h2>How close does each run come to the road?</h2>
<div class="imap-hist" id="imap-hist"><div class="imap-tip" id="imap-tip"></div></div>
<p class="sub">Closest approach of each run's flow (peak thickness &gt; {thr:g}&nbsp;m) to where the avalanche path crosses the
 road; each run is counted once. The first bar holds runs whose flow touches the road anywhere. Bars are grouped by friction
 setting.</p>
</div>
"""

JS = r"""
(async function () {
const D = JSON.parse(document.getElementById('imap-data').textContent);
const G = D.grid, E = D.extent, $ = id => document.getElementById(id);
const fmt = (v, d) => v == null ? '-' : v.toFixed(d);
async function unz(b64, T) {           // base64 gzip -> typed array
  const bin = Uint8Array.from(atob(b64), c => c.charCodeAt(0));
  const buf = await new Response(new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
  return new T(buf);
}
const img = src => new Promise(r => { const i = new Image(); i.onload = () => r(i); i.src = src; });
// decode everything once
const hill = await img(D.hill.png);
const count = await unz(D.count, Uint8Array);
for (const p of D.panels) {
  p.vals = {}; p.imgs = {};
  for (const f of D.fields) { p.vals[f.key] = await unz(p.data[f.key], Uint16Array); p.imgs[f.key] = await img(p.png[f.key]); }
}
let field = D.fields[0].key;
const F = () => D.fields.find(f => f.key === field);
const cell = (x, y) => {                // grid index under map point, or -1
  const c = Math.floor((x - G.left) / G.res), r = Math.floor((G.top - y) / G.res);
  return (c < 0 || r < 0 || c >= G.w || r >= G.h) ? -1 : r * G.w + c;
};
// distance reference: where the avalanche path crosses the road (fallback: nearest road point)
function nearestRoad(x, y) {
  let best = null, bd = Infinity;
  for (const L of D.roads) for (let i = 1; i < L.length; i++) {
    const [ax, ay] = L[i - 1], [bx, by] = L[i], dx = bx - ax, dy = by - ay;
    const t = Math.max(0, Math.min(1, ((x - ax) * dx + (y - ay) * dy) / (dx * dx + dy * dy || 1)));
    const px = ax + t * dx, py = ay + t * dy, d = Math.hypot(x - px, y - py);
    if (d < bd) { bd = d; best = [px, py]; }
  }
  return best;
}
const ref = (x, y) => D.path ? D.path.crossing : nearestRoad(x, y);

// ---- shared state: view {cx, cy, k} (k = px per metre, null = fit) and cursor in map metres
const view = { cx: (E.x0 + E.x1) / 2, cy: (E.y0 + E.y1) / 2, k: null };
let cursor = null, pinned = false;
const m = location.hash.match(/at=([\d.]+),([\d.]+)/);
if (m) { cursor = [+m[1], +m[2]]; pinned = true; }

// ---- field buttons
for (const f of D.fields) {
  const b = document.createElement('button'); b.textContent = f.name; b.dataset.k = f.key;
  b.onclick = () => { field = f.key; document.querySelectorAll('#imap-fields button').forEach(x => x.classList.toggle('on', x.dataset.k === field)); legend(); drawAll(); readout(); };
  $('imap-fields').appendChild(b);
}
$('imap-fields').querySelector('button').classList.add('on');

// ---- panels
const panels = D.panels.map(p => {
  const el = document.createElement('div'); el.className = 'imap-panel';
  el.innerHTML = `<div class="imap-title"><b>${p.title}</b> ${p.sub}<br>${p.status}</div><div class="imap-box"><canvas></canvas><div class="imap-tag"></div></div>`;
  $('imap-row').appendChild(el);
  const box = el.querySelector('.imap-box'), cv = el.querySelector('canvas');
  const P = { p, box, cv, tag: el.querySelector('.imap-tag'), w: 0, h: 0 };
  new ResizeObserver(() => { box.style.height = (box.clientWidth * (E.y1 - E.y0) / (E.x1 - E.x0)) + 'px'; draw(P); }).observe(box);
  hook(P);
  return P;
});
const fitK = P => Math.min(P.w / (E.x1 - E.x0), P.h / (E.y1 - E.y0));
const K = P => view.k ?? fitK(P);
const toPx = (P, x, y) => [(x - view.cx) * K(P) + P.w / 2, (view.cy - y) * K(P) + P.h / 2];
const toMap = (P, px, py) => [(px - P.w / 2) / K(P) + view.cx, view.cy - (py - P.h / 2) / K(P)];

function draw(P) {
  const dpr = window.devicePixelRatio || 1, w = P.cv.clientWidth, h = P.cv.clientHeight;
  if (!w || !h) return;
  if (P.cv.width !== Math.round(w * dpr) || P.cv.height !== Math.round(h * dpr)) { P.cv.width = Math.round(w * dpr); P.cv.height = Math.round(h * dpr); }
  P.w = w; P.h = h;
  const c = P.cv.getContext('2d'), k = K(P);
  c.setTransform(dpr, 0, 0, dpr, 0, 0); c.clearRect(0, 0, w, h);
  const place = (im, left, top, ww, hh) => { const [a, b] = toPx(P, left, top); c.drawImage(im, a, b, ww * k, hh * k); };
  c.imageSmoothingEnabled = true; place(hill, D.hill.left, D.hill.top, D.hill.w, D.hill.h);
  c.imageSmoothingEnabled = k * G.res < 2; place(P.p.imgs[field], G.left, G.top, G.w * G.res, G.h * G.res);
  const poly = (L, close) => { c.beginPath(); L.forEach(([x, y], i) => { const [a, b] = toPx(P, x, y); i ? c.lineTo(a, b) : c.moveTo(a, b); }); if (close) c.closePath(); c.stroke(); };
  c.setLineDash([]); c.lineJoin = 'round';
  c.strokeStyle = '#FF4FA0'; c.lineWidth = 1.8; D.release.forEach(L => poly(L, true));
  c.strokeStyle = '#E87722'; c.lineWidth = 2.5; D.roads.forEach(L => poly(L));
  if (D.path) {
    c.strokeStyle = '#fff'; c.lineWidth = 3; poly(D.path.line);
    c.strokeStyle = '#000'; c.lineWidth = 1.2; c.setLineDash([5, 4]); poly(D.path.line); c.setLineDash([]);
    const [a, b] = toPx(P, ...D.path.crossing);
    c.fillStyle = '#fff'; c.strokeStyle = '#E87722'; c.lineWidth = 2;
    c.beginPath(); c.arc(a, b, 6, 0, 7); c.fill();
    c.beginPath(); c.moveTo(a - 3.5, b - 3.5); c.lineTo(a + 3.5, b + 3.5); c.moveTo(a + 3.5, b - 3.5); c.lineTo(a - 3.5, b + 3.5); c.stroke();
  }
  // cursor: distance line to the road reference, crosshair, value tag
  const f = F();
  if (cursor) {
    const [px, py] = toPx(P, ...cursor), R = ref(...cursor);
    if (R) {
      const [rx, ry] = toPx(P, ...R), d = Math.hypot(cursor[0] - R[0], cursor[1] - R[1]);
      c.strokeStyle = '#fff'; c.lineWidth = 3; c.setLineDash([]); c.beginPath(); c.moveTo(px, py); c.lineTo(rx, ry); c.stroke();
      c.strokeStyle = '#000'; c.lineWidth = 1.2; c.setLineDash([4, 3]); c.beginPath(); c.moveTo(px, py); c.lineTo(rx, ry); c.stroke(); c.setLineDash([]);
      c.fillStyle = '#000'; c.beginPath(); c.arc(rx, ry, 3, 0, 7); c.fill();
      const lab = Math.round(d) + ' m'; c.font = '12px system-ui, sans-serif';
      const mx = (px + rx) / 2, my = (py + ry) / 2, tw = c.measureText(lab).width;
      c.fillStyle = 'rgba(255,255,255,.85)'; c.fillRect(mx - tw / 2 - 3, my - 9, tw + 6, 16); c.fillStyle = '#000'; c.fillText(lab, mx - tw / 2, my + 3);
    }
    c.strokeStyle = pinned ? '#A6192E' : '#000'; c.lineWidth = pinned ? 1.6 : 1;
    c.beginPath(); c.moveTo(px - 10, py); c.lineTo(px + 10, py); c.moveTo(px, py - 10); c.lineTo(px, py + 10); c.stroke();
    const i = cell(...cursor), v = i < 0 ? null : P.p.vals[field][i] / f.scale;
    P.tag.textContent = `${f.name.replace('Peak flow ', '')} = ${v == null ? '-' : fmt(v, f.dec) + ' ' + f.unit}`;
  } else P.tag.textContent = '';
}
const drawAll = () => panels.forEach(draw);

// ---- interaction: wheel zoom about the cursor, drag pan, double-click reset, click pin
function hook(P) {
  const cv = P.cv; let drag = null;
  const ev = e => { const r = cv.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
  cv.addEventListener('wheel', e => {
    e.preventDefault();
    const [px, py] = ev(e), [x, y] = toMap(P, px, py), f0 = fitK(P);
    const dy = e.deltaY * (e.deltaMode === 1 ? 33 : e.deltaMode === 2 ? 400 : 1);   // lines/pages -> px (Firefox)
    const k = Math.min(40, Math.max(0.8 * f0, K(P) * Math.pow(1.0015, -dy)));
    view.k = k; view.cx = x - (px - P.w / 2) / k; view.cy = y + (py - P.h / 2) / k; drawAll();
  }, { passive: false });
  cv.addEventListener('pointerdown', e => { drag = { p: ev(e), c: [view.cx, view.cy], moved: false }; cv.setPointerCapture(e.pointerId); });
  cv.addEventListener('pointermove', e => {
    const [px, py] = ev(e);
    if (drag) {
      const dx = px - drag.p[0], dy = py - drag.p[1];
      if (Math.hypot(dx, dy) > 3) drag.moved = true;
      if (drag.moved) { view.k = K(P); view.cx = drag.c[0] - dx / view.k; view.cy = drag.c[1] + dy / view.k; drawAll(); return; }
    }
    if (!pinned) { cursor = toMap(P, px, py); drawAll(); readout(); }
  });
  cv.addEventListener('pointerup', e => {
    if (drag && !drag.moved) {
      if (pinned) { pinned = false; history.replaceState(null, '', location.pathname + location.search); }
      else { cursor = toMap(P, ...ev(e)); pinned = true; history.replaceState(null, '', '#at=' + cursor.map(v => v.toFixed(0)).join(',')); }
      drawAll(); readout();
    }
    drag = null;
  });
  cv.addEventListener('pointerleave', () => { if (!pinned && !drag) { cursor = null; drawAll(); readout(); } });
  cv.addEventListener('dblclick', () => { view.k = null; view.cx = (E.x0 + E.x1) / 2; view.cy = (E.y0 + E.y1) / 2; drawAll(); });
}

// ---- readout table above the maps
function readout() {
  const i = cursor ? cell(...cursor) : -1;
  let h = '<table><tr><th class="l">Run</th>' + D.fields.map(f => `<th>${f.name.replace('Peak flow ', 'peak ')} (${f.unit})</th>`).join('') + '<th class="l">Road</th></tr>';
  for (const p of D.panels) {
    h += `<tr><td>${p.title}</td>` + D.fields.map(f => `<td>${i < 0 ? '-' : fmt(p.vals[f.key][i] / f.scale, f.dec)}</td>`).join('') + `<td class="l">${p.status}</td></tr>`;
  }
  h += '</table>';
  if (!cursor) h += '<span class="imap-help">Move over a map to read each run\'s peak values at that point and its distance to where the avalanche path crosses the road. Click to pin.</span>';
  else {
    const R = ref(...cursor), d = R ? Math.hypot(cursor[0] - R[0], cursor[1] - R[1]) : null;
    h += `<b>${d == null ? '-' : Math.round(d) + ' m'}</b> from ${D.path ? 'where the avalanche path crosses the road' : 'the nearest road point'}` +
         ` | E ${cursor[0].toFixed(0)}, N ${cursor[1].toFixed(0)} | flow reaches this cell in ${i < 0 ? '-' : count[i]} of ${D.countN} runs` + (pinned ? ' | <b style="color:#A6192E">pinned</b>' : '');
  }
  $('imap-read').innerHTML = h;
}

// ---- legend
function legend() {
  const f = F(), ticks = [0, .5, 1].map(t => fmt(t * f.vmax, f.dec)).join(' &nbsp; &hellip; &nbsp; ');
  $('imap-legend').innerHTML =
    `<span>${f.name} (${f.unit}): <span class="imap-ramp" style="background:linear-gradient(90deg,${D.ramp.join(',')})"></span> ${ticks}</span>` +
    `<span><svg width="26" height="10"><line x1="0" y1="5" x2="26" y2="5" stroke="#E87722" stroke-width="2.5"/></svg> road edges</span>` +
    `<span><svg width="26" height="10"><line x1="0" y1="5" x2="26" y2="5" stroke="#FF4FA0" stroke-width="2"/></svg> release area</span>` +
    `<span><svg width="44" height="14"><line x1="0" y1="7" x2="30" y2="7" stroke="#fff" stroke-width="3"/><line x1="0" y1="7" x2="30" y2="7" stroke="#000" stroke-width="1.2" stroke-dasharray="5 4"/>` +
    `<circle cx="37" cy="7" r="5.5" fill="#fff" stroke="#E87722" stroke-width="2"/></svg> avalanche path, road crossing</span>` +
    `<span class="imap-help">Scroll to zoom, drag to pan, double-click to reset, click to pin (the link then points to that spot).` +
    ` Values are cells with peak thickness &gt; ${D.thr} m.</span>`;
}

// ---- closest-approach histogram (inline SVG)
function histogram() {
  const host = $('imap-hist'); if (!host || !D.hist.length) return;
  const B = D.bin, all = D.hist.flatMap(g => g.scen).filter(s => s.d != null);
  const nb = Math.max(1, Math.ceil(Math.max(0, ...all.map(s => s.d)) / B));
  const slots = [{ lo: null }].concat([...Array(nb)].map((_, i) => ({ lo: i * B, hi: (i + 1) * B })));
  const inSlot = (s, b) => b.lo == null ? s.d === 0 : s.d > b.lo && s.d <= b.hi;   // 'reached' = 0, bins (lo, hi]
  const greys = ['#222', '#7A8893', '#C9CED2', '#4D5A63', '#A3ADB4'];
  const cnt = slots.map(b => D.hist.map(g => g.scen.filter(s => s.d != null && inSlot(s, b))));
  const ymax = Math.max(1, ...cnt.flat().map(a => a.length));
  const W = 860, H = 300, L = 40, Rm = 10, T = 34, Bm = 42, pw = (W - L - Rm) / slots.length, bw = Math.min(26, (pw - 8) / D.hist.length);
  const y = v => T + (H - T - Bm) * (1 - v / ymax);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Closest approach to the road crossing">`;
  for (let v = 0; v <= ymax; v++) s += `<line x1="${L}" x2="${W - Rm}" y1="${y(v)}" y2="${y(v)}" stroke="#e3e3e3"/><text x="${L - 6}" y="${y(v) + 4}" text-anchor="end">${v}</text>`;
  s += `<text x="12" y="${(T + H - Bm) / 2}" transform="rotate(-90 12 ${(T + H - Bm) / 2})" text-anchor="middle">number of runs</text>`;
  slots.forEach((b, i) => {
    const x0 = L + i * pw, g0 = x0 + (pw - bw * D.hist.length) / 2;
    cnt[i].forEach((arr, j) => { if (arr.length) s += `<rect x="${g0 + j * bw}" y="${y(arr.length)}" width="${bw - 2}" height="${y(0) - y(arr.length)}" fill="${b.lo == null ? '#A6192E' : greys[j % greys.length]}" ${b.lo == null && j ? `fill-opacity="${1 - j * .25}"` : ''}/>`; });
    s += `<text x="${x0 + pw / 2}" y="${H - Bm + 16}" text-anchor="middle">${b.lo == null ? 'reached' : b.hi}</text>`;
    s += `<rect class="hit" data-i="${i}" x="${x0}" y="${T}" width="${pw}" height="${H - T - Bm}" fill="transparent"/>`;
  });
  s += `<line x1="${L + pw}" x2="${L + pw}" y1="${T - 4}" y2="${H - Bm}" stroke="#E87722" stroke-width="2"/>`;
  s += `<text x="${(L + pw + W - Rm) / 2}" y="${H - 6}" text-anchor="middle">closest approach to the road crossing (m, bins labelled by upper edge)</text>`;
  D.hist.forEach((g, j) => { s += `<rect x="${L + 10 + j * 190}" y="6" width="12" height="12" fill="${greys[j % greys.length]}"/><text x="${L + 27 + j * 190}" y="16">${g.name}</text>`; });
  s += '</svg>';
  host.insertAdjacentHTML('afterbegin', s);
  const tip = $('imap-tip');
  host.querySelectorAll('rect.hit').forEach(r => {
    r.addEventListener('mousemove', e => {
      const i = +r.dataset.i, b = slots[i];
      const rng = b.lo == null ? 'flow reaches the road' : `${b.lo}&ndash;${b.hi} m short of the road crossing`;
      tip.innerHTML = `<b>${rng}</b><br>` + D.hist.map((g, j) => {
        const a = cnt[i][j], n = g.scen.length;
        return `${g.name}: ${a.length} of ${n} (${Math.round(100 * a.length / n)}%)${a.length ? ': ' + a.map(x => x.id).join(', ') : ''}`;
      }).join('<br>');
      const hr = host.getBoundingClientRect();
      tip.style.display = 'block'; tip.style.left = Math.min(e.clientX - hr.left + 12, hr.width - tip.offsetWidth - 4) + 'px'; tip.style.top = (e.clientY - hr.top + 12) + 'px';
    });
    r.addEventListener('mouseleave', () => tip.style.display = 'none');
  });
}

legend(); readout(); drawAll(); histogram();
window.addEventListener('resize', drawAll);
})();
"""


def section(printPng, **kw):
    """CSS, markup, JSON payload and script for the report body. kw: see payload(). Returns (css, html, histHtml)."""
    data = json.dumps(payload(**kw), separators=(",", ":")).replace("</", "<\\/")
    html = MARKUP.format(printPng=printPng) + f'<script type="application/json" id="imap-data">{data}</script>'
    return CSS, html, HIST.format(thr=kw["thr"]) + f"<script>{JS}</script>"
