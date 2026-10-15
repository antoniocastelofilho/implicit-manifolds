import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle, os; sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from polfile import read_pol, polygons
from zbuf import render, crop
plt.rcParams.update({'font.size': 9, 'mathtext.fontset': 'cm'})
R = REPO + '/book/runs/ml/'; FG = REPO + '/book/figures/'
BLUE = '#2a78d6'; ORANGE = '#eb6834'; GREEN = '#1baf7a'


def load(tag):
    D = pickle.load(open(R + tag + '.pkl', 'rb'))
    n, k, Dv, cells = read_pol(D['cells_path'])
    if D['trim'] is not None:
        from scipy.spatial import cKDTree
        dd, _ = cKDTree(D['X']).query(np.array([c.coords.mean(0) for c in cells])); cells = [c for c, x in zip(cells, dd) if x < D['trim']]
    return D, cells


def polys(cells): return [c.coords[cyc] for c in cells for cyc in polygons(c)]
def cols(v, cmap, lo=None, hi=None):
    v = np.asarray(v); lo = v.min() if lo is None else lo; hi = v.max() if hi is None else hi
    cm = plt.get_cmap(cmap); return [cm(np.clip((x - lo) / (hi - lo + 1e-12), 0, 1))[:3] for x in v]


def spheres_render(Pts, color, size=(900, 800), elev=30, azim=-60, r=0.012, R_=None, bbox=None):
    """render points as small octahedra (cheap 'splats') with the z-buffer"""
    o = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float) * r
    F = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4), (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    P = []; C = []
    for p, c in zip(Pts, color):
        for f in F: P.append(p + o[list(f)]); C.append(c)
    return render(P, C, elev=elev, azim=azim, size=size, edges=False, R=R_)


which = sys.argv[1:]
# ---------------------------------------------------------------- 2D illustration
if 'curve' in which:
    from pui_fast import PUIFast
    D, cells = load('curve'); X = D['X']
    dmin = np.array([-1.63, -1.61]); dmax = np.array([1.62, 1.64]); div = np.array([20, 20])
    M = PUIFast(X, dmin, dmax, div, 1); ev = M.njit_eval()
    xs = np.linspace(dmin[0], dmax[0], 400); ys = np.linspace(dmin[1], dmax[1], 400); XX, YY = np.meshgrid(xs, ys)
    Z = np.zeros_like(XX); y = np.zeros(1)
    for i in range(XX.shape[0]):
        for j in range(XX.shape[1]):
            ev(np.array([XX[i, j], YY[i, j]]), y); Z[i, j] = y[0]
    h = M.h
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.9))
    ax = axs[0]
    for i in range(int(div[0]) + 1): ax.axvline(dmin[0] + i * h[0], color='0.88', lw=0.5)
    for j in range(int(div[1]) + 1): ax.axhline(dmin[1] + j * h[1], color='0.88', lw=0.5)
    for c in M.cells[:M.ndata]:
        x0 = dmin + np.array(c) * h; ax.add_patch(plt.Rectangle(x0, h[0], h[1], color=BLUE, alpha=0.22, lw=0))
    ax.plot(X[:, 0], X[:, 1], '.', ms=1.2, color='0.25')
    # a few local approximations, drawn over their supports
    rng = np.random.default_rng(3)
    for i in rng.choice(M.ndata, 9, replace=False):
        p, T, N, cf = M.P[i], M.T[i][:, 0], M.N[i][:, 0], M.coef[i, 0]
        u = np.linspace(-1.5 * h[0], 1.5 * h[0], 40); q = cf[0] + cf[1] * u + cf[2] * u * u
        pts = p + np.outer(u, T) + np.outer(q, N); ax.plot(pts[:, 0], pts[:, 1], color=ORANGE, lw=1.6)
        ctr = dmin + (np.array(M.cells[i]) + 0.5) * h
        ax.add_patch(plt.Rectangle(ctr - 1.5 * h, 3 * h[0], 3 * h[1], fill=False, ec=ORANGE, lw=0.6, ls='--'))
    ax.set_title('(a) sample, data cells and local fits', fontsize=9)
    ax = axs[1]
    lev = np.linspace(-0.5, 0.5, 21)
    cs = ax.contourf(XX, YY, Z, levels=lev, cmap='RdBu_r', extend='both'); ax.contour(XX, YY, Z, [0], colors='k', linewidths=1.2)
    ax.set_title('(b) partition of unity implicit $F$', fontsize=9)
    ax = axs[2]
    P = D['P']
    for c in cells: ax.plot(c.coords[:, 0], c.coords[:, 1], color=BLUE, lw=1.2)
    sel = slice(0, None, 4)
    ax.plot(X[sel, 0], X[sel, 1], '.', ms=1.5, color='0.55')
    for a, b in zip(X[sel], P[sel]): ax.plot([a[0], b[0]], [a[1], b[1]], color=ORANGE, lw=0.5)
    ax.set_title(f'(c) approximation ({len(cells)} segments) and projections', fontsize=9)
    for a in axs: a.set_aspect('equal'); a.set_xlim(-1.5, 1.5); a.set_ylim(-1.5, 1.5); a.set_xticks([]); a.set_yticks([])
    axs[2].set_xlim(0.2, 1.4); axs[2].set_ylim(-0.6, 0.6)
    plt.subplots_adjust(0.01, 0.01, 0.99, 0.93, 0.05, 0); plt.savefig(FG + 'ml_pui2d.pdf'); plt.close()
# ---------------------------------------------------------------- surface in R^3
if 'surface' in which:
    D, cells = load('surface'); X = D['X']; P = D['P']
    rng = np.random.default_rng(0); s = rng.choice(len(X), 9000, replace=False)
    im1 = crop(spheres_render(X[s], cols(X[s, 2], 'viridis'), elev=35, azim=-60, r=0.011))
    Q = polys(cells)
    im2 = crop(render(Q, cols([q[:, 2].mean() for q in Q], 'viridis'), elev=35, azim=-60, size=(900, 800), edge_alpha=0.5))
    im3 = crop(spheres_render(P[s], cols(P[s, 2], 'viridis'), elev=35, azim=-60, r=0.011))
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.6))
    tt = [f'(a) sample ({len(X)} points, noise $0.01$)', f'(b) PUI and FCH: {len(cells)} polygons, $\\chi=0$', '(c) sample projected on the polygons']
    for a, im, t in zip(axs, [im1, im2, im3], tt): a.imshow(im); a.axis('off'); a.set_title(t, fontsize=9)
    plt.subplots_adjust(0.01, 0.01, 0.99, 0.9, 0.03, 0); plt.savefig(FG + 'ml_surface.png', dpi=200); plt.close()
# ---------------------------------------------------------------- tori in R^4 and R^6 + spectra
if 'tori' in which:
    D4, c4 = load('torus4'); D6, c6 = load('torus6')
    def stereo(Q, r): return r * Q[:, :3] / (r - Q[:, 3:4])
    Q = polys(c4); r0 = np.sqrt(1 + 0.36)
    im1 = crop(render([stereo(q, r0) for q in Q], cols([np.arctan2(q[:, 1].mean(), q[:, 0].mean()) for q in Q], 'twilight'), elev=40, azim=-60, size=(900, 800), edge_alpha=0.5))
    Q = polys(c6)
    def tmap(q): return np.c_[(2 + 0.6 * q[:, 2]) * q[:, 0], (2 + 0.6 * q[:, 2]) * q[:, 1], 0.6 * q[:, 3]]
    im2 = crop(render([tmap(q) for q in Q], cols([np.arctan2(q[:, 5].mean(), q[:, 4].mean()) for q in Q], 'twilight'), elev=40, azim=-60, size=(900, 800), edge_alpha=0.5))
    fig = plt.figure(figsize=(12, 3.7)); gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.9])
    a = fig.add_subplot(gs[0]); a.imshow(im1); a.axis('off'); a.set_title(f'(a) $T^2\\subset\\mathbb{{R}}^4$: {len(c4)} polygons, $\\chi=0$', fontsize=9)
    a = fig.add_subplot(gs[1]); a.imshow(im2); a.axis('off'); a.set_title(f'(b) $T^2\\subset\\mathbb{{R}}^6$: {len(c6)} polygons, $\\chi=0$', fontsize=9)
    a = fig.add_subplot(gs[2])
    for D, lab, c in [(D4, '$\\mathbb{R}^4$', BLUE), (D6, '$\\mathbb{R}^6$', ORANGE)]:
        E = D['eig'] / D['eig'][:, :1]
        med = np.median(E, 0); lo = np.percentile(E, 10, 0); hi = np.percentile(E, 90, 0)
        x = np.arange(1, len(med) + 1)
        a.semilogy(x, med, '-o', color=c, ms=4, lw=2, label='sample in ' + lab); a.fill_between(x, lo, hi, color=c, alpha=0.15, lw=0)
    a.set_xlabel('index $i$ of the eigenvalue'); a.set_ylabel('$\\lambda_i/\\lambda_1$ (median, 10–90%)'); a.legend(frameon=False, fontsize=8)
    a.set_title('(c) local PCA spectra: two large eigenvalues', fontsize=9); a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
    a.axvline(2.5, color='0.6', lw=0.8, ls='--')
    plt.subplots_adjust(0.01, 0.13, 0.99, 0.9, 0.3, 0); plt.savefig(FG + 'ml_tori.png', dpi=200); plt.close()
# ---------------------------------------------------------------- Swiss roll: unrolling by geodesics on the polygons
if 'swiss' in which:
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra
    from scipy.spatial import cKDTree
    D, cells = load('swiss'); X = D['X']; P = D['P']; par = np.load(R + 'swiss_param.npy')
    key = {}; V = []; E = set()
    for c in cells:
        ids = []
        for p in c.coords:
            kk = tuple(np.round(p, 6))
            if kk not in key: key[kk] = len(V); V.append(p)
            ids.append(key[kk])
        # all pairs of vertices of each polygon: the graph metric then approaches the intrinsic one
        for cyc in polygons(c):
            for a_ in range(len(cyc)):
                for b_ in range(a_ + 1, len(cyc)): E.add(tuple(sorted((ids[cyc[a_]], ids[cyc[b_]]))))
    V = np.array(V); E = np.array(sorted(E)); w = np.linalg.norm(V[E[:, 0]] - V[E[:, 1]], axis=1)
    G = coo_matrix((np.r_[w, w], (np.r_[E[:, 0], E[:, 1]], np.r_[E[:, 1], E[:, 0]])), shape=(len(V), len(V))).tocsr()
    rng = np.random.default_rng(0); L = rng.choice(len(V), 250, replace=False)
    Dl = dijkstra(G, indices=L)                                   # landmarks x vertices
    Dl[~np.isfinite(Dl)] = Dl[np.isfinite(Dl)].max()
    # landmark MDS
    D2 = Dl[:, L] ** 2; nL = len(L); J = np.eye(nL) - 1 / nL
    B = -0.5 * J @ D2 @ J; lam, U = np.linalg.eigh(B); lam = lam[::-1][:2]; U = U[:, ::-1][:, :2]
    Lh = U / np.sqrt(lam); mu = D2.mean(1)
    Y = -0.5 * ((Dl ** 2).T - mu) @ Lh                           # all vertices
    _, nv = cKDTree(V).query(P); emb_mesh = Y[nv]
    # ground truth: arc length of the spiral and height
    t = par[:, 0]; s_true = 0.5 * (t * np.sqrt(1 + t * t) + np.arcsinh(t)) / 15
    truth = np.c_[s_true, 1.6 * par[:, 1]]
    def fit_err(Z):
        Zc = Z - Z.mean(0); Tc = truth - truth.mean(0)
        U_, S_, Vt_ = np.linalg.svd(Zc.T @ Tc); Rr = U_ @ Vt_
        return np.sqrt(np.mean(np.sum((Zc @ Rr - Tc) ** 2, 1))) / np.sqrt(np.mean(np.sum(Tc ** 2, 1)))
    from sklearn.manifold import Isomap
    sub = rng.choice(len(X), 3000, replace=False)
    iso = Isomap(n_neighbors=10, n_components=2).fit_transform(X[sub])
    emb_iso = np.zeros((len(X), 2)); emb_iso[:] = np.nan; emb_iso[sub] = iso
    def fit_err_sub(Z, idx):
        Zc = Z[idx] - Z[idx].mean(0); Tc = truth[idx] - truth[idx].mean(0)
        U_, S_, Vt_ = np.linalg.svd(Zc.T @ Tc); Rr = U_ @ Vt_
        return np.sqrt(np.mean(np.sum((Zc @ Rr - Tc) ** 2, 1))) / np.sqrt(np.mean(np.sum(Tc ** 2, 1)))
    iso_p = Isomap(n_neighbors=10, n_components=2).fit_transform(P[sub])
    emb_isop = np.zeros((len(X), 2)); emb_isop[sub] = iso_p
    e_mesh = fit_err(emb_mesh); e_iso = fit_err_sub(emb_iso, sub); e_isop = fit_err_sub(emb_isop, sub)
    # stronger noise: the same comparison when the sample is noisier
    print('isomap on projected points %.3f' % e_isop)
    pickle.dump(dict(e_mesh=e_mesh, e_iso=e_iso, e_isop=e_isop, nV=len(V), nE=len(E)), open(R + 'swiss_unroll.pkl', 'wb'))
    print('relative error: mesh geodesics %.3f, isomap %.3f' % (e_mesh, e_iso))
    Q = polys(cells)
    im = crop(render(Q, cols([np.arctan2(q[:, 2].mean(), q[:, 0].mean()) for q in Q], 'Spectral'), elev=15, azim=-70, size=(900, 800), edge_alpha=0.4))
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.4))
    axs[0].imshow(im); axs[0].axis('off'); axs[0].set_title(f'(a) PUI + FCH: {len(cells)} polygons, $\\chi=1$', fontsize=9)
    cmap = plt.get_cmap('Spectral')
    axs[1].scatter(emb_mesh[::3, 0], emb_mesh[::3, 1], c=t[::3], s=1.2, cmap=cmap); axs[1].set_title(f'(b) geodesics on the polygons ({e_mesh:.3f})', fontsize=9)
    axs[2].scatter(iso_p[:, 0], iso_p[:, 1], c=t[sub], s=1.5, cmap=cmap); axs[2].set_title(f'(c) Isomap, projected points ({e_isop:.3f})', fontsize=9)
    axs[3].scatter(iso[:, 0], iso[:, 1], c=t[sub], s=1.5, cmap=cmap); axs[3].set_title(f'(d) Isomap, sample ({e_iso:.3f})', fontsize=9)
    for a in axs[1:]: a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([])
    plt.subplots_adjust(0.01, 0.02, 0.99, 0.9, 0.06, 0); plt.savefig(FG + 'ml_swiss.png', dpi=200); plt.close()
# ---------------------------------------------------------------- persistence
if 'tda' in which:
    T = pickle.load(open(R + 'tda.pkl', 'rb'))
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.3))
    colors = {0: BLUE, 1: ORANGE, 2: GREEN}
    for a, (tag, title) in zip(axs, [('curve', 'closed curve in $\\mathbb{R}^2$'), ('torus4', 'torus in $\\mathbb{R}^4$'), ('swiss', 'Swiss roll in $\\mathbb{R}^3$')]):
        dg = T[tag]['dgm']; me = T[tag]['max_edge']
        for d in (0, 1, 2):
            pts = np.array([(b_, (e_ if np.isfinite(e_) else me * 1.05)) for dd, (b_, e_) in dg if dd == d])
            if len(pts): a.scatter(pts[:, 0], pts[:, 1], s=9 if d else 5, color=colors[d], label=f'$H_{d}$', zorder=3 - d * 0)
        a.plot([0, me * 1.05], [0, me * 1.05], color='0.6', lw=0.8); a.axhline(me * 1.05, color='0.8', lw=0.6, ls='--')
        a.set_xlim(-0.02 * me, me * 1.08); a.set_ylim(-0.02 * me, me * 1.1); a.set_xlabel('birth'); a.set_ylabel('death')
        a.set_title(f'{title}: Betti of the cells {tuple(T[tag]["betti"])}', fontsize=9)
        a.spines['top'].set_visible(False); a.spines['right'].set_visible(False); a.set_aspect('equal')
    axs[1].legend(frameon=False, fontsize=8, loc='lower right')
    plt.tight_layout(); plt.savefig(FG + 'ml_tda.pdf'); plt.close()
