"""Baselines for the paper: two-dimensional reductions of four data sets by PCA, kernel PCA, ms-gLPCA (a simplification
of kNN-tPCA: graph-Laplacian PCA with the average Laplacian of kNN graphs over a range of k), Isomap, LLE,
Laplacian eigenmaps, t-SNE and UMAP.  Metrics: trustworthiness T and continuity C (10 neighbours), Spearman correlation of
pairwise distances, and for the digits the 5-fold accuracy of a 10-NN classifier on the two coordinates."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, pickle; sys.path.insert(0, REPO + '/python/implicit_manifolds/learning')
import numpy as np, warnings; warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import rotate, gaussian_filter
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
from scipy.sparse import csgraph
from sklearn.manifold import trustworthiness, Isomap, TSNE, LocallyLinearEmbedding, SpectralEmbedding
from sklearn.decomposition import PCA, KernelPCA
from sklearn.datasets import load_digits
from sklearn.neighbors import KNeighborsClassifier, kneighbors_graph
from sklearn.model_selection import cross_val_score
import umap
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
R = REPO + '/book/runs/ml/'; FG = REPO + '/book/figures/paper/'
rng = np.random.default_rng(5)


def tpca_f(X, d=2, ks=(5, 10, 15, 20, 30), beta=0.5):
    """min ||X^T - U Q^T||_F^2 + alpha Tr(Q^T PL Q), Q^T Q = I, with PL = mean_t L_t over kNN graphs (Jiang et al. 2013
    closed form, normalized as G = (1-beta)(I - XX^T/lambda_max) + beta PL/xi_max): Q = eigenvectors of the d smallest eigenvalues"""
    Xc = X - X.mean(0); Gm = Xc @ Xc.T; lam = np.linalg.eigvalsh(Gm)[-1]
    PL = np.zeros_like(Gm)
    for k in ks:
        A = kneighbors_graph(X, k, include_self=False); A = 0.5 * (A + A.T)
        PL += csgraph.laplacian(A).toarray() / len(ks)
    xi = np.linalg.eigvalsh(PL)[-1]
    G = (1 - beta) * (np.eye(len(X)) - Gm / lam) + beta * PL / xi
    w, V = np.linalg.eigh(G)
    return V[:, :d]


data = {}
t = rng.uniform(1.5 * np.pi, 4.5 * np.pi, 2000); yy = rng.uniform(0, 1, 2000)
data['swiss'] = (np.c_[t * np.cos(t) / 15, 1.6 * yy - 0.8, t * np.sin(t) / 15] + 0.008 * rng.normal(size=(2000, 3)), t, 'Swiss roll, $\\mathbb{R}^3$')
u = rng.uniform(0, 2 * np.pi, 2000); v = rng.uniform(0, 2 * np.pi, 2000)
data['torus'] = (np.c_[np.cos(u), np.sin(u), np.cos(v), np.sin(v), .5 * np.cos(u + v), .5 * np.sin(u + v)] + 0.01 * rng.normal(size=(2000, 6)), u, 'torus, $\\mathbb{R}^6$')
base = np.zeros((32, 32)); base[8:24, 14:18] = 1; base[8:12, 14:24] = 1; base = gaussian_filter(base, 1.0)
ang = np.sort(rng.uniform(0, 360, 720))
imgs = np.array([rotate(base, a, reshape=False, order=1).ravel() for a in ang]) + 0.01 * rng.normal(size=(720, 1024))
data['images'] = (imgs, ang, 'rotating images, $\\mathbb{R}^{1024}$')
dg = load_digits(); data['digits'] = (dg.data / 16.0, dg.target, 'digits, $\\mathbb{R}^{64}$')
pickle.dump({k: (v[0], v[1]) for k, v in data.items()}, open(R + 'paper_data.pkl', 'wb'))

meths = {'PCA': lambda Z: PCA(2).fit_transform(Z),
         'kernel PCA': lambda Z: KernelPCA(2, kernel='rbf', gamma=1.0 / (np.median(pdist(Z)) ** 2)).fit_transform(Z),
         'ms-gLPCA': lambda Z: tpca_f(Z),
         'Isomap': lambda Z: Isomap(n_neighbors=10, n_components=2).fit_transform(Z),
         'LLE': lambda Z: LocallyLinearEmbedding(n_neighbors=12, n_components=2, random_state=0).fit_transform(Z),
         'Lapl. eigenmaps': lambda Z: SpectralEmbedding(n_components=2, n_neighbors=12, random_state=0).fit_transform(Z),
         't-SNE': lambda Z: TSNE(2, perplexity=30, random_state=0, init='pca').fit_transform(Z),
         'UMAP': lambda Z: umap.UMAP(n_neighbors=15, random_state=0).fit_transform(Z)}
res = {}; emb = {}
for dn, (X, c, title) in data.items():
    D0 = pdist(X)
    for mn, f in meths.items():
        t0 = time.time(); Y = f(X); dt = time.time() - t0; emb[(dn, mn)] = Y
        T = trustworthiness(X, Y, n_neighbors=10); C = trustworthiness(Y, X, n_neighbors=10)
        rho = spearmanr(D0, pdist(Y)).statistic
        acc = cross_val_score(KNeighborsClassifier(10), Y, dg.target, cv=5).mean() if dn == 'digits' else None
        res[(dn, mn)] = dict(T=T, C=C, rho=rho, acc=acc, time=dt); print(dn, mn, res[(dn, mn)], flush=True)
    if dn == 'digits':
        res[(dn, 'full')] = dict(acc=cross_val_score(KNeighborsClassifier(10), X, dg.target, cv=5).mean())
pickle.dump(dict(res=res, emb=emb), open(R + 'paper_baselines.pkl', 'wb'))
show = ['PCA', 'ms-gLPCA', 'Isomap', 't-SNE', 'UMAP']
fig, axs = plt.subplots(4, 5, figsize=(10, 8.6))
for i, (dn, (X, c, title)) in enumerate(data.items()):
    for j, mn in enumerate(show):
        a = axs[i, j]; Y = emb[(dn, mn)]
        cm = 'tab10' if dn == 'digits' else ('twilight' if dn in ('torus', 'images') else 'Spectral')
        a.scatter(Y[:, 0], Y[:, 1], c=c, s=1.2 if dn != 'images' else 3, cmap=cm, rasterized=True)
        r = res[(dn, mn)]; lab = f'T={r["T"]:.2f} C={r["C"]:.2f}' + (f' acc={r["acc"]:.2f}' if r['acc'] is not None else '')
        a.set_title(lab, fontsize=7.5); a.set_xticks([]); a.set_yticks([])
        if i == 0: a.text(0.5, 1.25, mn, transform=a.transAxes, ha='center', fontsize=9.5, weight='bold')
        if j == 0: a.set_ylabel(title, fontsize=8.5)
plt.tight_layout(); plt.savefig(FG + 'baselines.pdf', dpi=200); plt.close()
