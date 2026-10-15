import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
R=REPO + '/book/runs/ml/'; FG=REPO + '/book/figures/paper/'
P=pickle.load(open(R+'paper_baselines.pkl','rb')); res=P['res']; emb=P['emb']
for k in list(res):
    if k[1]=='tPCA-F': res[(k[0],'ms-gLPCA')]=res[k]
for k in list(emb):
    if k[1]=='tPCA-F': emb[(k[0],'ms-gLPCA')]=emb[k]
from sklearn.datasets import load_digits
D=pickle.load(open(R+'paper_data.pkl','rb'))
titles={'swiss':'Swiss roll, $\\mathbb{R}^3$','torus':'torus, $\\mathbb{R}^6$','images':'rotating images, $\\mathbb{R}^{1024}$','digits':'digits, $\\mathbb{R}^{64}$'}
data={k:(D[k][0],D[k][1],titles[k]) for k in ['swiss','torus','images','digits']}
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
