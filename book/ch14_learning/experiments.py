import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, os, time, pickle
sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np
from pui import PUI
from pui_fast import PUIFast
from runfch import run_fch
from polfile import read_pol, euler
from project import project
R = REPO + '/book/runs/ml/'; os.makedirs(R, exist_ok=True)
rng = np.random.default_rng(2026)
res = {}


def experiment(tag, X, truth_dist, m, dmin, dmax, div_pui, refine=2, trim=None):
    n = X.shape[1]; k = n - m
    t0 = time.time(); M = PUIFast(X, dmin, dmax, div_pui, m); t1 = time.time()
    shim, ev = M.cfunc()
    y = np.zeros(k); ev(X[0], y)                       # compile
    # starting point: the sample point with the smallest |F|
    vals = []
    for x in X[:200]:
        ev(x, y); vals.append(np.abs(y).max())
    x0 = X[int(np.argmin(vals))]
    div_fch = np.asarray(div_pui) * refine
    pol, tm = run_fch(shim, n, k, dmin, dmax, div_fch, x0, tag)
    n_, k_, D, cells = read_pol(pol)
    if trim is not None:                               # discard cells far from the data (open manifolds)
        from scipy.spatial import cKDTree
        tree = cKDTree(X); d, _ = tree.query(np.array([c.coords.mean(0) for c in cells]))
        cells = [c for c, dd in zip(cells, d) if dd < trim]
    chi = euler(cells)
    h = (np.asarray(dmax) - np.asarray(dmin)) / div_fch
    t2 = time.time(); P, dist = project(X, cells, m, np.asarray(dmin, float), h); t3 = time.time()
    r = dict(n=n, m=m, k=k, npts=len(X), pui_cells=len(M.cells), div_pui=list(div_pui), div_fch=list(div_fch),
             t_pui=t1 - t0, t_fch=tm, t_proj=t3 - t2, cells=len(cells), chi=chi,
             err_before=truth_dist(X), err_after=truth_dist(P), dist_proj=dist, eig=M.eig)
    res[tag] = r
    pickle.dump(dict(X=X, P=P, eig=M.eig, cells_path=pol, trim=trim, r={kk: v for kk, v in r.items()}),
                open(R + tag + '.pkl', 'wb'))
    print(tag, 'pui cells', len(M.cells), 'cells', len(cells), 'chi', chi[1],
          'err before %.4f after %.4f (rms)' % (np.sqrt(np.mean(r['err_before'] ** 2)), np.sqrt(np.mean(r['err_after'] ** 2))),
          'times pui %.1f fch %.2f proj %.1f' % (t1 - t0, tm['fch'], t3 - t2), flush=True)
    return r


which = sys.argv[1:]
sig = 0.01
if 'curve' in which:
    th = rng.uniform(0, 2 * np.pi, 1500); rr = 1 + 0.25 * np.cos(3 * th)
    X = np.c_[rr * np.cos(th), rr * np.sin(th)] + 0.015 * rng.normal(size=(1500, 2))
    def td(X):
        a = np.arctan2(X[:, 1], X[:, 0]); best = np.full(len(X), 1e9)
        tt = np.linspace(0, 2 * np.pi, 4000); C = np.c_[(1 + 0.25 * np.cos(3 * tt)) * np.cos(tt), (1 + 0.25 * np.cos(3 * tt)) * np.sin(tt)]
        from scipy.spatial import cKDTree
        return cKDTree(C).query(X)[0]
    experiment('curve', X, td, 1, [-1.63, -1.61], [1.62, 1.64], [20, 20], refine=2)
if 'surface' in which:
    N = 30000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    rr = 0.35 * (1 + 0.25 * np.sin(5 * u))
    X = np.c_[(1 + rr * np.cos(v)) * np.cos(u), (1 + rr * np.cos(v)) * np.sin(u), rr * np.sin(v)] + sig * rng.normal(size=(N, 3))
    tt = np.linspace(0, 2 * np.pi, 600); U, V = np.meshgrid(tt, tt); Rr = 0.35 * (1 + 0.25 * np.sin(5 * U))
    G = np.c_[((1 + Rr * np.cos(V)) * np.cos(U)).ravel(), ((1 + Rr * np.cos(V)) * np.sin(U)).ravel(), (Rr * np.sin(V)).ravel()]
    from scipy.spatial import cKDTree; TG = cKDTree(G)
    experiment('surface', X, lambda Y: TG.query(Y)[0], 2, [-1.53, -1.52, -0.51], [1.51, 1.52, 0.52], [36, 36, 12], refine=1)
if 'torus4' in which:
    N = 20000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    X = np.c_[np.cos(u), np.sin(u), 0.6 * np.cos(v), 0.6 * np.sin(v)] + sig * rng.normal(size=(N, 4))
    td = lambda Y: np.hypot(np.hypot(Y[:, 0], Y[:, 1]) - 1, np.hypot(Y[:, 2], Y[:, 3]) - 0.6)
    experiment('torus4', X, td, 2, [-1.27, -1.26, -0.83, -0.82], [1.25, 1.26, 0.81, 0.82], [12, 12, 8, 8], refine=2)
if 'torus6' in which:
    N = 40000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    def emb(u, v): return np.c_[np.cos(u), np.sin(u), np.cos(v), np.sin(v), 0.5 * np.cos(u + v), 0.5 * np.sin(u + v)]
    X = emb(u, v) + sig * rng.normal(size=(N, 6))
    tt = np.linspace(0, 2 * np.pi, 500); U, V = np.meshgrid(tt, tt)
    from scipy.spatial import cKDTree; TG = cKDTree(emb(U.ravel(), V.ravel()))
    experiment('torus6', X, lambda Y: TG.query(Y)[0], 2, [-1.17, -1.16, -1.15, -1.14, -0.63, -0.62], [1.13, 1.14, 1.15, 1.16, 0.61, 0.62],
               [12, 12, 12, 12, 6, 6], refine=1)
if 'swiss' in which:
    N = 12000; t = rng.uniform(1.5 * np.pi, 4.5 * np.pi, N); yy = rng.uniform(0, 1, N)
    def emb(t, yy): return np.c_[t * np.cos(t) / 15, 1.6 * yy - 0.8, t * np.sin(t) / 15]
    X = emb(t, yy) + 0.008 * rng.normal(size=(N, 3))
    tt = np.linspace(1.5 * np.pi, 4.5 * np.pi, 1500); Y_ = np.linspace(0, 1, 200); T_, YY = np.meshgrid(tt, Y_)
    from scipy.spatial import cKDTree; TG = cKDTree(emb(T_.ravel(), YY.ravel()))
    np.save(R + 'swiss_param.npy', np.c_[t, yy])
    experiment('swiss', X, lambda Y: TG.query(Y)[0], 2, [-0.93, -0.92, -0.91], [0.91, 0.92, 0.93], [22, 20, 22], refine=2, trim=0.08)
if 'swiss2' in which:
    N = 12000; t = rng.uniform(1.5 * np.pi, 4.5 * np.pi, N); yy = rng.uniform(0, 1, N)
    def emb(t, yy): return np.c_[t * np.cos(t) / 15, 1.6 * yy - 0.8, t * np.sin(t) / 15]
    X = emb(t, yy) + 0.03 * rng.normal(size=(N, 3))
    tt = np.linspace(1.5 * np.pi, 4.5 * np.pi, 1500); Y_ = np.linspace(0, 1, 200); T_, YY = np.meshgrid(tt, Y_)
    from scipy.spatial import cKDTree; TG = cKDTree(emb(T_.ravel(), YY.ravel()))
    np.save(R + 'swiss2_param.npy', np.c_[t, yy])
    experiment('swiss2', X, lambda Y: TG.query(Y)[0], 2, [-0.93, -0.92, -0.91], [0.91, 0.92, 0.93], [22, 20, 22], refine=2, trim=0.1)
pickle.dump(res, open(R + 'res_' + '_'.join(which) + '.pkl', 'wb'))
