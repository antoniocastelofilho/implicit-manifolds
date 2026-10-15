"""Compare polynomial, RBF and Hermite RBF local approximations in the partition of unity.
Metrics: distance of the vertices of the FCH approximation to the exact manifold (rms, max), Euler
characteristic, construction and evaluation times."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, pickle
sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np
from pui_fast import PUIFast
from pui_rbf import PUIRBF
from runfch import run_fch
from polfile import read_pol, euler
from scipy.spatial import cKDTree
R = REPO + '/book/runs/ml/'


def run(tag, X, td, m, dmin, dmax, div, kind, refine=1, **kw):
    n = X.shape[1]; k = n - m
    t0 = time.time()
    M = PUIFast(X, dmin, dmax, div, m) if kind == 'poly' else PUIRBF(X, dmin, dmax, div, m, kind=kind, **kw)
    tb = time.time() - t0
    shim, ev = M.cfunc(); y = np.zeros(k); ev(X[0], y)
    Z = X[:20000]; t0 = time.time()
    for x in Z: ev(x, y)
    te = (time.time() - t0) / len(Z) * 1e6
    vals = []
    for x in X[:200]:
        ev(x, y); vals.append(np.abs(y).max())
    pol, tm = run_fch(shim, n, k, dmin, dmax, np.asarray(div) * refine, X[int(np.argmin(vals))], 'rbf_' + tag + '_' + kind)
    _, _, _, cells = read_pol(pol)
    V = np.unique(np.round(np.concatenate([c.coords for c in cells]), 9), axis=0); d = td(V)
    W = newton(ev, V, k); dz = td(W)
    r = dict(kind=kind, cells=len(cells), chi=euler(cells)[1], rms=np.sqrt(np.mean(d ** 2)), mx=d.max(), rms_z=np.sqrt(np.mean(dz ** 2)), mx_z=dz.max(), t_build=tb,
             t_eval=te, t_fch=tm['fch'], rbf_cells=getattr(M, 'rbf_cells', 0), data=M.ndata)
    print(tag, r, flush=True); return r


def newton(ev, V, k, it=4, e=1e-6):
    """Gauss-Newton projection of the vertices onto F^{-1}(0): measures the error of F itself"""
    W = V.copy(); n = V.shape[1]; y = np.zeros(k); y2 = np.zeros(k)
    for i in range(len(W)):
        x = W[i].copy()
        for _ in range(it):
            ev(x, y); J = np.empty((k, n))
            for a in range(n):
                xa = x.copy(); xa[a] += e; ev(xa, y2); J[:, a] = (y2 - y) / e
            x = x - np.linalg.pinv(J) @ y
        W[i] = x
    return W


def param_dist(S, ngrid=800):
    """exact distance to a parametric surface S(u,v) (2*pi periodic): nearest grid point, then Gauss-Newton"""
    tt = np.linspace(0, 2 * np.pi, ngrid, endpoint=False); U0, V0 = np.meshgrid(tt, tt); U0 = U0.ravel(); V0 = V0.ravel()
    tree = cKDTree(S(U0, V0))
    def td(Y):
        _, i = tree.query(Y); u = U0[i].copy(); v = V0[i].copy(); e = 1e-6
        for _ in range(8):
            P = S(u, v); Su = (S(u + e, v) - P) / e; Sv = (S(u, v + e) - P) / e; r = Y - P
            a11 = (Su * Su).sum(1); a12 = (Su * Sv).sum(1); a22 = (Sv * Sv).sum(1); b1 = (Su * r).sum(1); b2 = (Sv * r).sum(1)
            det = a11 * a22 - a12 ** 2; u += (a22 * b1 - a12 * b2) / det; v += (a11 * b2 - a12 * b1) / det
        return np.linalg.norm(Y - S(u, v), axis=1)
    return td


rng = np.random.default_rng(7)
out = {}
which = sys.argv[1]; sig = float(sys.argv[2]) if len(sys.argv) > 2 else 0.002
kinds = sys.argv[3].split(',') if len(sys.argv) > 3 else ['poly', 'rbf', 'hrbf']
if which == 'torus4':
    N = 20000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    X = np.c_[np.cos(u), np.sin(u), 0.6 * np.cos(v), 0.6 * np.sin(v)] + sig * rng.normal(size=(N, 4))
    td = lambda Y: np.hypot(np.hypot(Y[:, 0], Y[:, 1]) - 1, np.hypot(Y[:, 2], Y[:, 3]) - 0.6)
    for g in [(6, 6, 4, 4), (8, 8, 6, 6), (12, 12, 8, 8)]:
        for kind in kinds:
            out[(g, kind)] = run('t4_%d' % g[0], X, td, 2, [-1.27, -1.26, -0.83, -0.82], [1.25, 1.26, 0.81, 0.82], g, kind, refine=2)
if which == 'surface':
    N = 30000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    rr = 0.35 * (1 + 0.25 * np.sin(5 * u))
    X = np.c_[(1 + rr * np.cos(v)) * np.cos(u), (1 + rr * np.cos(v)) * np.sin(u), rr * np.sin(v)] + sig * rng.normal(size=(N, 3))
    tt = np.linspace(0, 2 * np.pi, 1500); U, V = np.meshgrid(tt, tt); Rr = 0.35 * (1 + 0.25 * np.sin(5 * U))
    G = np.c_[((1 + Rr * np.cos(V)) * np.cos(U)).ravel(), ((1 + Rr * np.cos(V)) * np.sin(U)).ravel(), (Rr * np.sin(V)).ravel()]
    TG = cKDTree(G)
    for g in [(18, 18, 6), (24, 24, 8), (36, 36, 12)]:
        for kind in kinds:
            out[(g, kind)] = run('s_%d' % g[0], X, lambda Y: TG.query(Y)[0], 2, [-1.53, -1.52, -0.51], [1.51, 1.52, 0.52], g, kind, refine=2)
if which in ('ripple3', 'ripple4'):
    N = 30000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    if which == 'ripple3':      # torus of R^3 whose tube radius oscillates 12 times around the core circle
        def S(u, v):
            rr = 0.35 * (1 + 0.08 * np.sin(12 * u))
            return np.c_[(1 + rr * np.cos(v)) * np.cos(u), (1 + rr * np.cos(v)) * np.sin(u), rr * np.sin(v)]
        box = ([-1.53, -1.52, -0.51], [1.51, 1.52, 0.52]); grids = [(16, 16, 6), (24, 24, 8), (36, 36, 12)]
    else:                        # torus of R^4 whose second radius oscillates 9 times
        def S(u, v):
            rho = 0.6 + 0.04 * np.sin(9 * u)
            return np.c_[np.cos(u), np.sin(u), rho * np.cos(v), rho * np.sin(v)]
        box = ([-1.27, -1.26, -0.86, -0.85], [1.25, 1.26, 0.85, 0.86]); grids = [(8, 8, 6, 6), (12, 12, 8, 8), (16, 16, 10, 10)]
    X = S(u, v) + sig * rng.normal(size=(N, S(u[:1], v[:1]).shape[1])); td = param_dist(S)
    for g in grids:
        for kind in kinds:
            out[(g, kind)] = run(which + '_%d' % g[0], X, td, 2, box[0], box[1], g, kind, refine=2)
pickle.dump(out, open(R + 'rbf_%s_%g.pkl' % (which, sig), 'wb'))
