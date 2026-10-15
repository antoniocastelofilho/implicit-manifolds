"""Where dimensionality reduction fails: linear methods on nonlinear data, and every 2-D embedding when the
topology of the data manifold does not fit in the plane."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, pickle; sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import rotate, gaussian_filter
from sklearn.manifold import trustworthiness, Isomap, TSNE
from sklearn.decomposition import PCA
from sklearn.datasets import load_digits
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
import gudhi
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
R = REPO + '/book/runs/ml/'; FG = REPO + '/book/figures/'
rng = np.random.default_rng(5)
data = {}
# Swiss roll
t = rng.uniform(1.5 * np.pi, 4.5 * np.pi, 2000); yy = rng.uniform(0, 1, 2000)
data['swiss'] = (np.c_[t * np.cos(t) / 15, 1.6 * yy - 0.8, t * np.sin(t) / 15] + 0.008 * rng.normal(size=(2000, 3)), t, 'Swiss roll, $\\mathbb{R}^3$')
# torus of R^6
u = rng.uniform(0, 2 * np.pi, 2000); v = rng.uniform(0, 2 * np.pi, 2000)
data['torus'] = (np.c_[np.cos(u), np.sin(u), np.cos(v), np.sin(v), .5 * np.cos(u + v), .5 * np.sin(u + v)] + 0.01 * rng.normal(size=(2000, 6)), u, 'torus, $\\mathbb{R}^6$')
# images of a rotating object: a closed curve in R^1024
base = np.zeros((32, 32)); base[8:24, 14:18] = 1; base[8:12, 14:24] = 1; base = gaussian_filter(base, 1.0)
ang = np.sort(rng.uniform(0, 360, 720))
imgs = np.array([rotate(base, a, reshape=False, order=1).ravel() for a in ang]) + 0.01 * rng.normal(size=(720, 1024))
data['images'] = (imgs, ang, 'rotating images, $\\mathbb{R}^{1024}$')
# handwritten digits
dg = load_digits(); data['digits'] = (dg.data / 16.0, dg.target, 'digits, $\\mathbb{R}^{64}$')
meths = {'PCA': lambda Z: PCA(2).fit_transform(Z),
         'Isomap': lambda Z: Isomap(n_neighbors=10, n_components=2).fit_transform(Z),
         't-SNE': lambda Z: TSNE(2, perplexity=30, random_state=0, init='pca').fit_transform(Z)}
res = {}; emb = {}
for dn, (X, c, title) in data.items():
    for mn, f in meths.items():
        Y = f(X); emb[(dn, mn)] = Y
        T = trustworthiness(X, Y, n_neighbors=10); C = trustworthiness(Y, X, n_neighbors=10)
        extra = None
        if dn == 'digits':
            extra = cross_val_score(KNeighborsClassifier(10), Y, dg.target, cv=5).mean()
        res[(dn, mn)] = (T, C, extra); print(dn, mn, res[(dn, mn)], flush=True)
    if dn == 'digits':
        res[(dn, 'full')] = (1, 1, cross_val_score(KNeighborsClassifier(10), X, dg.target, cv=5).mean())
# global PCA spectra (number of components needed) and persistent H1 of the images
spec = {dn: PCA().fit(X).explained_variance_ratio_ for dn, (X, c, t_) in data.items()}
Xi = data['images'][0]; idx = rng.choice(len(Xi), 300, replace=False)
st = gudhi.RipsComplex(points=Xi[idx], max_edge_length=12.0).create_simplex_tree(max_dimension=2)
dgm = st.persistence(); h1 = sorted([d - b for k, (b, d) in dgm if k == 1], reverse=True)[:3]
info = dict(res=res, spec={k: v[:10] for k, v in spec.items()}, h1=h1)
pickle.dump(info, open(R + 'dimred_fail.pkl', 'wb')); print('H1 bars', h1); print({k: np.round(v[:8], 3) for k, v in spec.items()})
# figure
fig, axs = plt.subplots(4, 3, figsize=(9, 11.5))
for i, (dn, (X, c, title)) in enumerate(data.items()):
    for j, mn in enumerate(meths):
        a = axs[i, j]; Y = emb[(dn, mn)]
        cm = 'tab10' if dn == 'digits' else ('twilight' if dn in ('torus', 'images') else 'Spectral')
        a.scatter(Y[:, 0], Y[:, 1], c=c, s=1.5 if dn != 'images' else 3, cmap=cm)
        T, C, e = res[(dn, mn)]
        lab = f'T={T:.2f}, C={C:.2f}' + (f', kNN acc.={e:.2f}' if e is not None else '')
        a.set_title((title + '\n' if j == 0 else '\n') + f'{mn}: {lab}', fontsize=8)
        a.set_xticks([]); a.set_yticks([])
plt.tight_layout(); plt.savefig(FG + 'ml_dimred_fail.png', dpi=180); plt.close()
