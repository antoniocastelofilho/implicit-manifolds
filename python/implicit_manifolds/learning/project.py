"""Projection of points of R^n onto the polygons of a piecewise linear surface (m = 2) or the
segments of a curve (m = 1), given as cells of a .pol file: closest point by brute force over the
simplices whose bucket is within one bucket of the point."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np, itertools
from numba import njit
from polfile import polygons


def simplices_of(cells, m):
    S = []
    for c in cells:
        if m == 1:
            for e in c.faces[1]: S.append(c.coords[list(e)])
        else:
            for cyc in polygons(c):
                Q = c.coords[cyc]; g = Q.mean(0)
                for i in range(len(Q)): S.append(np.array([g, Q[i], Q[(i + 1) % len(Q)]]))
    return np.array(S)


@njit(cache=True)
def _closest_seg(p, a, b):
    d = b - a; t = np.dot(p - a, d) / max(np.dot(d, d), 1e-300)
    t = min(1.0, max(0.0, t)); return a + t * d


@njit(cache=True)
def _closest_tri(p, a, b, c):
    e0 = b - a; e1 = c - a; v = p - a
    A00 = np.dot(e0, e0); A01 = np.dot(e0, e1); A11 = np.dot(e1, e1)
    b0 = np.dot(v, e0); b1 = np.dot(v, e1); det = A00 * A11 - A01 * A01
    if det > 1e-300:
        s = (A11 * b0 - A01 * b1) / det; t = (A00 * b1 - A01 * b0) / det
        if s >= 0 and t >= 0 and s + t <= 1: return a + s * e0 + t * e1
    q1 = _closest_seg(p, a, b); q2 = _closest_seg(p, b, c); q3 = _closest_seg(p, a, c)
    d1 = np.sum((p - q1) ** 2); d2 = np.sum((p - q2) ** 2); d3 = np.sum((p - q3) ** 2)
    if d1 <= d2 and d1 <= d3: return q1
    if d2 <= d3: return q2
    return q3


@njit(cache=True)
def _project(X, S, order, start, keys_sorted, pkeys_nb, m):
    out = np.empty_like(X); dist = np.empty(X.shape[0])
    for i in range(X.shape[0]):
        best = 1e300; bq = X[i].copy()
        for kk in range(pkeys_nb.shape[1]):
            key = pkeys_nb[i, kk]
            lo = np.searchsorted(keys_sorted, key); hi = np.searchsorted(keys_sorted, key, 'right')
            for t in range(lo, hi):
                s = order[t]
                if m == 1: q = _closest_seg(X[i], S[s, 0], S[s, 1])
                else: q = _closest_tri(X[i], S[s, 0], S[s, 1], S[s, 2])
                d = np.sum((X[i] - q) ** 2)
                if d < best: best = d; bq = q
        out[i] = bq; dist[i] = np.sqrt(best)
    return out, dist


def project(X, cells, m, dmin, h):
    S = simplices_of(cells, m)
    n = X.shape[1]
    big = np.int64(1 << 16)
    def key(I): return np.sum(I.astype(np.int64) * big ** np.arange(n), axis=-1)
    skeys = key(np.floor((S.mean(1) - dmin) / h))
    order = np.argsort(skeys); ks = skeys[order]
    offs = np.array(list(itertools.product((-1, 0, 1), repeat=n)))
    I = np.floor((X - dmin) / h).astype(np.int64)
    pk = key(I[:, None, :] + offs[None, :, :])
    return _project(np.ascontiguousarray(X), np.ascontiguousarray(S), order, 0, ks, pk, m)
