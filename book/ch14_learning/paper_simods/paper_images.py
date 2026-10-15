"""Rotating images: PCA to d coordinates, partition of unity implicit of a curve (m=1) in R^d, Fast Continuation
Hypercubes, Betti numbers of the reconstruction, and recovery of the rotation angle by arc length along the reconstructed
closed curve, compared with the angle read from two-dimensional embeddings."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, pickle; sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap, TSNE
import umap
from pui_fast import PUIFast
from runfch import run_fch
from polfile import read_pol, euler
from project import project
from tda import betti_of_cells
R = REPO + '/book/runs/ml/'
D = pickle.load(open(R + 'paper_data.pkl', 'rb')); X, ang = D['images']
d = int(sys.argv[1]) if len(sys.argv) > 1 else 4; ndiv = int(sys.argv[2]) if len(sys.argv) > 2 else 16


def angle_error(theta, ang):
    """mean absolute circular error (degrees) after the best rotation and orientation"""
    best = 1e9
    for s in (1, -1):
        dif = (s * theta - ang + 180) % 360 - 180
        c = np.angle(np.mean(np.exp(1j * np.radians(dif))), deg=True)
        e = np.mean(np.abs((dif - c + 180) % 360 - 180)); best = min(best, e)
    return best


t0 = time.time(); pca = PCA(d).fit(X); Z = pca.transform(X); tp = time.time() - t0
lo = Z.min(0) - 0.25 * np.ptp(Z, 0) / ndiv * 3; hi = Z.max(0) + 0.25 * np.ptp(Z, 0) / ndiv * 3
lo = lo - 0.0137 * np.arange(1, d + 1) / d
div = np.full(d, ndiv)
t0 = time.time(); M = PUIFast(Z, lo, hi, div, 1); tb = time.time() - t0
shim, ev = M.cfunc(); y = np.zeros(d - 1); vals = []
for z in Z:
    ev(z, y); vals.append(np.abs(y).max())
pol, tm = run_fch(shim, d, d - 1, lo, hi, div * 2, Z[int(np.argmin(vals))], 'images_d%d' % d)
_, _, _, cells = read_pol(pol)
chi = euler(cells); b = betti_of_cells(cells, 1)
h = (hi - lo) / (div * 2)
P, dist = project(Z, cells, 1, lo, h)
# arc length along the reconstructed cycle
segs = [c.coords for c in cells]
V, inv = np.unique(np.round(np.concatenate(segs), 9), axis=0, return_inverse=True); E = inv.reshape(-1, 2)
adj = {}
for a_, b_ in E: adj.setdefault(a_, []).append(b_); adj.setdefault(b_, []).append(a_)
order = [E[0, 0]]; prev = None
while True:
    cur = order[-1]; nx = [w for w in adj[cur] if w != prev]
    if not nx or nx[0] == order[0]: break
    prev = cur; order.append(nx[0])
C = V[order]; L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(np.r_[C, C[:1]], axis=0), axis=1))]
s = np.empty(len(P))
for i, p in enumerate(P):
    A = C; B = np.r_[C[1:], C[:1]]; AB = B - A; tt = np.clip(np.sum((p - A) * AB, 1) / np.sum(AB * AB, 1), 0, 1)
    dd = np.linalg.norm(A + tt[:, None] * AB - p, axis=1); j = int(np.argmin(dd)); s[i] = L[j] + tt[j] * (L[j + 1] - L[j])
theta = 360 * s / L[-1]
err = {'ours (arc length)': angle_error(theta, ang)}
for name, Y in [('PCA', PCA(2).fit_transform(X)), ('Isomap', Isomap(n_neighbors=10, n_components=2).fit_transform(X)),
                ('t-SNE', TSNE(2, perplexity=30, random_state=0, init='pca').fit_transform(X)),
                ('UMAP', umap.UMAP(n_neighbors=15, random_state=0).fit_transform(X))]:
    Yc = Y - Y.mean(0); err[name] = angle_error(np.degrees(np.arctan2(Yc[:, 1], Yc[:, 0])), ang)
out = dict(d=d, ndiv=ndiv, ev=pca.explained_variance_ratio_, cells=len(cells), loops=len(order), chi=chi, betti=b, t_pca=tp,
           t_pui=tb, t_fch=tm['fch'], dist=dist, err=err, theta=theta, ang=ang, C=C, Z=Z, P=P)
pickle.dump(out, open(R + 'paper_images_d%d.pkl' % d, 'wb'))
print(d, ndiv, 'cells', len(cells), 'cycle', len(order), 'chi', chi, 'betti', b, 'err', {k: round(v, 2) for k, v in err.items()},
      'times', round(tp, 2), round(tb, 2), round(tm['fch'], 3), 'proj dist', np.sqrt(np.mean(dist ** 2)))
