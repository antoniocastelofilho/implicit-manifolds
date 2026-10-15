"""Dimensionality reduction experiments for Chapter 14.
(1) Swiss roll: embeddings by PCA, Isomap, LLE, Laplacian eigenmaps, t-SNE, Isomap on the points projected
onto the reconstructed polygons, and MDS of geodesics on the polygons; trustworthiness, continuity,
Spearman correlation of distances, relative error with respect to the true coordinates, time.
(2) A torus of R^6 hidden in R^20: PCA finds the 6-dimensional span, the partition of unity and FCH
reconstruct the torus in the reduced space."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, pickle; sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
from sklearn.manifold import trustworthiness, Isomap, LocallyLinearEmbedding, SpectralEmbedding, TSNE
from sklearn.decomposition import PCA
R = REPO + '/book/runs/ml/'
rng = np.random.default_rng(11)
which = sys.argv[1:]

def metrics(X, Y, truth, k=12):
    T = trustworthiness(X, Y, n_neighbors=k); C = trustworthiness(Y, X, n_neighbors=k)
    rho = spearmanr(pdist(X), pdist(Y)).statistic
    Zc = Y - Y.mean(0); Tc = truth - truth.mean(0); U, S, Vt = np.linalg.svd(Zc.T @ Tc)
    Zr = Zc @ (U @ Vt); s = np.sum(Zr * Tc) / np.sum(Zr * Zr)          # best similarity transform
    err = np.sqrt(np.mean(np.sum((s * Zr - Tc) ** 2, 1))) / np.sqrt(np.mean(np.sum(Tc ** 2, 1)))
    return T, C, rho, err

if 'swiss' in which:
    D = pickle.load(open(R + 'swiss.pkl', 'rb')); X = D['X']; P = D['P']; par = np.load(R + 'swiss_param.npy')
    t = par[:, 0]; s_true = 0.5 * (t * np.sqrt(1 + t * t) + np.arcsinh(t)) / 15; truth = np.c_[s_true, 1.6 * par[:, 1]]
    sub = rng.choice(len(X), 2000, replace=False); Xs = X[sub]; Ps = P[sub]; Ts = truth[sub]
    U = pickle.load(open(R + 'swiss_unroll_emb.pkl', 'rb')) if False else None
    meths = {
        'PCA': lambda Z: PCA(2).fit_transform(Z),
        'Isomap': lambda Z: Isomap(n_neighbors=10, n_components=2).fit_transform(Z),
        'LLE': lambda Z: LocallyLinearEmbedding(n_neighbors=12, n_components=2, random_state=0).fit_transform(Z),
        'Laplacian eigenmaps': lambda Z: SpectralEmbedding(n_components=2, n_neighbors=12, random_state=0).fit_transform(Z),
        't-SNE': lambda Z: TSNE(2, perplexity=30, random_state=0, init='pca').fit_transform(Z),
    }
    res = {}; emb = {}
    for name, f in meths.items():
        t0 = time.time(); Y = f(Xs); dt = time.time() - t0
        res[name] = metrics(Xs, Y, Ts) + (dt,); emb[name] = Y; print(name, res[name], flush=True)
    t0 = time.time(); Y = Isomap(n_neighbors=10, n_components=2).fit_transform(Ps); dt = time.time() - t0
    res['Isomap, projected points'] = metrics(Xs, Y, Ts) + (dt,); emb['Isomap, projected points'] = Y; print(res['Isomap, projected points'])
    pickle.dump(dict(res=res, emb=emb, t=t[sub]), open(R + 'dimred_swiss.pkl', 'wb'))

if 'hidden' in which:
    from pui_fast import PUIFast
    from runfch import run_fch
    from polfile import read_pol, euler
    from tda import betti_of_cells
    N = 40000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    X6 = np.c_[np.cos(u), np.sin(u), np.cos(v), np.sin(v), 0.5 * np.cos(u + v), 0.5 * np.sin(u + v)]
    Q, _ = np.linalg.qr(rng.normal(size=(20, 20)))
    X20 = np.c_[X6, np.zeros((N, 14))] @ Q.T + 0.01 * rng.normal(size=(N, 20)) + rng.normal(size=20)
    t0 = time.time(); pca = PCA(20).fit(X20); tp = time.time() - t0
    ev = pca.explained_variance_
    Z = pca.transform(X20)[:, :6]
    lo = Z.min(0) - 0.15; hi = Z.max(0) + 0.15
    ext = hi - lo; h0 = 2.3 / 12
    div = np.maximum(np.round(ext / h0).astype(int), 4)
    lo = lo - 0.0137 * np.arange(1, 7) / 6; hi = lo + div * h0                # asymmetric box, uniform spacing
    M = PUIFast(Z, lo, hi, div, 2); evf = M.njit_eval(); shim, _ = M.cfunc()
    y = np.zeros(4); vals = []
    for x in Z[:300]:
        evf(x, y); vals.append(np.linalg.norm(y))
    pol, tm = run_fch(shim, 6, 4, lo, hi, div, Z[int(np.argmin(vals))], 'hidden20')
    n_, k_, Dv, cells = read_pol(pol); chi = euler(cells); b = betti_of_cells(cells, 2)
    out = dict(ev=ev, t_pca=tp, div=div, cells=len(cells), chi=chi, betti=b, fch=tm, pui=M.timings)
    print(out); pickle.dump(out, open(R + 'dimred_hidden.pkl', 'wb'))
