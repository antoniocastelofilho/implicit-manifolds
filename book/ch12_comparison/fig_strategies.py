"""Marching vs continuation, simplices vs squares, on a curve with two components in the plane."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from collections import deque
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
def F(x, y): return np.minimum((x + 0.9) ** 2 / 0.55 + (y - 0.15) ** 2 / 0.35 - 1, ((x - 1.05) ** 2 + (y + 0.25) ** 2) / 0.28 - 1)
N = 14; lo, hi = -2.0, 2.0; h = (hi - lo) / N
xs = lo + h * np.arange(N + 1); V = F(xs[:, None], xs[None, :])
def tris(i, j):   # Kuhn triangles of square (i,j)
    a, b, c, d = (i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)
    return [(a, b, c), (a, c, d)]
def crossed(vs): s = [V[v] > 0 for v in vs]; return any(s) and not all(s)
def zero_pts(vs):
    P = []
    for p, q in [(vs[k], vs[(k + 1) % len(vs)]) for k in range(len(vs))]:
        fp, fq = V[p], V[q]
        if (fp > 0) != (fq > 0):
            t = fp / (fp - fq); P.append((xs[p[0]] + t * (xs[q[0]] - xs[p[0]]), xs[p[1]] + t * (xs[q[1]] - xs[p[1]])))
    return P
def sq(i, j): return [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
seed = (int((-0.9 + 0.74 - lo) / h), int((0.15 - lo) / h))
# continuation over squares (BFS through crossed edges) and over triangles
def cont_squares():
    vis = set(); crossedset = set()
    s = None
    for di in range(-2, 3):
        for dj in range(-2, 3):
            c = (seed[0] + di, seed[1] + dj)
            if crossed(sq(*c)): s = c; break
        if s: break
    Q = deque([s]); vis.add(s)
    while Q:
        c = Q.popleft(); crossedset.add(c)
        for d in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            nb = (c[0] + d[0], c[1] + d[1])
            if 0 <= nb[0] < N and 0 <= nb[1] < N and nb not in vis:
                vis.add(nb)
                if crossed(sq(*nb)): Q.append(nb)
    return vis, crossedset
def cont_tris():
    alltri = {}
    for i in range(N):
        for j in range(N):
            for t in tris(i, j): alltri[t] = (i, j)
    edges = {}
    for t in alltri:
        for e in [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])]:
            edges.setdefault(tuple(sorted(e)), []).append(t)
    s = next(t for t in alltri if crossed(t) and abs(alltri[t][0] - seed[0]) <= 2 and abs(alltri[t][1] - seed[1]) <= 2)
    vis = {s}; Q = deque([s]); cr = set()
    while Q:
        t = Q.popleft(); cr.add(t)
        for e in [(t[0], t[1]), (t[1], t[2]), (t[0], t[2])]:
            if (V[e[0]] > 0) == (V[e[1]] > 0): continue      # continuation only crosses transverse edges
            for u in edges[tuple(sorted(e))]:
                if u not in vis: vis.add(u); Q.append(u)
    return vis, cr
fig, axs = plt.subplots(1, 4, figsize=(10, 2.8))
titles = ['(a) Marching Simplex', '(b) Continuation Simplex', '(c) Marching Hypercubes', '(d) Continuation Hypercubes']
ctv, ctc = cont_tris(); csv, csc = cont_squares()
for a, tt in zip(axs, titles):
    a.set_title(tt, fontsize=8.5); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_xlim(lo, hi); a.set_ylim(lo, hi)
for i in range(N):
    for j in range(N):
        S = sq(i, j); P = [(xs[p[0]], xs[p[1]]) for p in S]
        # marching squares: every square visited
        axs[2].add_patch(Polygon(P, fc='#ececec' if not crossed(S) else '#f4b183', ec='#999999', lw=0.4))
        axs[3].add_patch(Polygon(P, fc=('#f4b183' if (i, j) in csc else ('#ececec' if (i, j) in csv else 'white')), ec='#999999', lw=0.4))
        for t in tris(i, j):
            Pt = [(xs[p[0]], xs[p[1]]) for p in t]
            axs[0].add_patch(Polygon(Pt, fc='#ececec' if not crossed(t) else '#9dc3e6', ec='#999999', lw=0.3))
            axs[1].add_patch(Polygon(Pt, fc=('#9dc3e6' if t in ctc else 'white'), ec='#999999', lw=0.3))
            for k, sel in ((0, crossed(t)), (1, t in ctc)):
                if sel:
                    Z = zero_pts(t)
                    if len(Z) == 2: axs[k].plot(*zip(*Z), 'k-', lw=1.0)
        for k, sel in ((2, crossed(S)), (3, (i, j) in csc)):
            if sel:
                Z = zero_pts(S)
                for q in range(0, len(Z) - 1, 2): axs[k].plot(*zip(Z[q], Z[q + 1]), 'k-', lw=1.0)
for a in (axs[1], axs[3]):
    a.plot(xs[seed[0]] + h / 2, xs[seed[1]] + h / 2, 'r*', ms=8)
plt.tight_layout(); plt.savefig(REPO + '/book/figures/cmp_strategies.pdf'); plt.close()
print('marching squares', N * N, 'continuation visited', len(csv), 'crossed', len(csc), '| triangles', 2 * N * N, 'cont visited', len(ctv), 'crossed', len(ctc))
