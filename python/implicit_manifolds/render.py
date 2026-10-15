"""Rendering of combinatorial-skeleton (.pol) files with matplotlib.

polygons(...)  : 2-cells of a 2-dimensional approximation (ordered vertex cycles)
slice3(...)    : slice of a 3-dimensional approximation by a hyperplane x_axis = c
draw_polys(...) : shaded 3D drawing after a linear projection R^n -> R^3
"""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
sys.path.insert(0, REPO + "/python/implicit_manifolds")
from polfile import read_pol, polygons as cell_polygons


def _key(p, nd=6):
    return tuple(np.round(p, nd))


def polygons(cells, dedupe=True):
    out, seen = [], set()
    for c in cells:
        for cyc in cell_polygons(c):
            P = c.coords[cyc]
            if dedupe:
                key = frozenset(_key(p) for p in P)
                if key in seen:
                    continue
                seen.add(key)
            out.append(P)
    return out


def segments(cells, dedupe=True):
    out, seen = [], set()
    for c in cells:
        if len(c.faces) < 2:
            continue
        for e in c.faces[1]:
            P = c.coords[e]
            key = frozenset(_key(p) for p in P)
            if dedupe and key in seen:
                continue
            seen.add(key)
            out.append(P)
    return out


def _tri_slice(T, axis, c):
    """intersection segment of triangle T (3 x n) with x_axis = c"""
    s = T[:, axis] - c
    pts = []
    for i in range(3):
        a, b = T[i], T[(i + 1) % 3]
        sa, sb = s[i], s[(i + 1) % 3]
        if (sa < 0) != (sb < 0):
            t = sa / (sa - sb)
            pts.append(a + t * (b - a))
    return pts if len(pts) == 2 else None


def slice3(cells, axis, c):
    """slice every 3-cell by the hyperplane x_axis = c; return polygons"""
    out = []
    for cell in cells:
        if len(cell.faces) < 4:
            continue
        polys = cell_polygons(cell)
        segs = []
        for cyc in polys:
            P = cell.coords[cyc]
            g = P.mean(axis=0)
            for i in range(len(P)):
                seg = _tri_slice(np.array([g, P[i], P[(i + 1) % len(P)]]), axis, c)
                if seg:
                    segs.append(seg)
        if len(segs) < 3:
            continue
        # chain segments into loops by matching endpoints
        pts = [segs[0][0], segs[0][1]]
        used = {0}
        while len(used) < len(segs):
            last = pts[-1]
            best, bj, bk = None, None, None
            for j, sg in enumerate(segs):
                if j in used:
                    continue
                for kk in (0, 1):
                    d = np.linalg.norm(sg[kk] - last)
                    if best is None or d < best:
                        best, bj, bk = d, j, kk
            used.add(bj)
            pts.append(segs[bj][1 - bk])
        out.append(np.array(pts[:-1]))
    return out


def fan(P):
    g = P.mean(axis=0)
    return [np.array([g, P[i], P[(i + 1) % len(P)]]) for i in range(len(P))]


def newell(Q):
    n = np.zeros(3)
    for i in range(len(Q)):
        a, b = Q[i], Q[(i + 1) % len(Q)]
        n += np.array([(a[1]-b[1])*(a[2]+b[2]), (a[2]-b[2])*(a[0]+b[0]), (a[0]-b[0])*(a[1]+b[1])])
    return n


def draw_polys(ax, polys, proj, color_values=None, cmap="viridis", edge=True,
               edgecolor=(0.12, 0.12, 0.12, 0.7), lw=0.3, alpha=1.0,
               light=(0.4, 0.3, 0.85), base_color=(0.36, 0.55, 0.85), vmin=None, vmax=None,
               zoom=1.35):
    proj = np.asarray(proj, float)
    L = np.asarray(light, float); L /= np.linalg.norm(L)
    faces, cols = [], []
    cm = plt.get_cmap(cmap)
    if color_values is not None:
        vmin = min(color_values) if vmin is None else vmin
        vmax = max(color_values) if vmax is None else vmax
    for i, P in enumerate(polys):
        Q = P @ proj.T
        nrm = newell(Q); nn = np.linalg.norm(nrm)
        shade = 0.5 + 0.5 * abs(nrm @ L) / nn if nn > 0 else 0.7
        if color_values is not None:
            c = np.array(cm((color_values[i] - vmin) / (vmax - vmin + 1e-12)))[:3]
        else:
            c = np.array(base_color)
        faces.append(Q)
        cols.append(tuple(np.clip(c * shade, 0, 1)) + (alpha,))
    pc = Poly3DCollection(faces, facecolors=cols,
                          edgecolors=edgecolor if edge else "none", linewidths=lw if edge else 0)
    ax.add_collection3d(pc)
    allp = np.vstack(faces)
    _equal(ax, allp, zoom)


def draw_segments(ax, segs, proj, color="C0", lw=0.6, color_values=None, cmap="viridis"):
    proj = np.asarray(proj, float)
    S = [np.array(s) @ proj.T for s in segs]
    if color_values is not None:
        cm = plt.get_cmap(cmap)
        v = np.asarray(color_values); v = (v - v.min()) / (np.ptp(v) + 1e-12)
        cols = [cm(x) for x in v]
    else:
        cols = [color]
    ax.add_collection3d(Line3DCollection(S, colors=cols, linewidths=lw))
    _equal(ax, np.vstack(S))


def _equal(ax, P, zoom=1.35):
    lo, hi = P.min(0), P.max(0)
    c, r = (lo + hi) / 2, (hi - lo).max() / 2
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)


def new_ax(fig, pos=111, elev=22, azim=-60):
    ax = fig.add_subplot(pos, projection="3d")
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    return ax


def axes_proj(n, idx=(0, 1, 2)):
    P = np.zeros((3, n))
    for r, i in enumerate(idx):
        P[r, i] = 1
    return P
