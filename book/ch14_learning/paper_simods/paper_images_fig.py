import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import pickle, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import rotate, gaussian_filter
from sklearn.decomposition import PCA
import umap
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
R = REPO + '/book/runs/ml/'; o = pickle.load(open(R + 'paper_images_d4.pkl', 'rb'))
D = pickle.load(open(R + 'paper_data.pkl', 'rb')); X, ang = D['images']
fig = plt.figure(figsize=(10, 3.1))
base = np.zeros((32, 32)); base[8:24, 14:18] = 1; base[8:12, 14:24] = 1; base = gaussian_filter(base, 1.0)
for i, a in enumerate([0, 60, 120, 180, 240, 300]):
    ax = fig.add_axes([0.005 + 0.07 * (i % 2), 0.68 - 0.31 * (i // 2), 0.065, 0.27]); ax.imshow(rotate(base, a, reshape=False, order=1), cmap='gray_r'); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f'{a}°', fontsize=7, pad=1)
fig.text(0.075, 0.02, '(a)', ha='center')
ax = fig.add_axes([0.17, 0.12, 0.2, 0.78]); Y = PCA(2).fit_transform(X)
ax.scatter(Y[:, 0], Y[:, 1], c=ang, cmap='twilight', s=3); ax.set_xticks([]); ax.set_yticks([]); ax.set_title('(b) PCA, two coordinates', fontsize=8.5)
ax = fig.add_axes([0.38, 0.05, 0.27, 0.9], projection='3d'); C = o['C']; Z = o['Z']
ax.scatter(Z[:, 0], Z[:, 2], Z[:, 3], c=ang, cmap='twilight', s=2, alpha=.5)
CC = np.r_[C, C[:1]]; ax.plot(CC[:, 0], CC[:, 2], CC[:, 3], 'k-', lw=.8)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.set_title('(c) reconstruction in $\\mathbb{R}^4$ (coords 1, 3, 4)', fontsize=8.5)
ax = fig.add_axes([0.71, 0.14, 0.27, 0.74])
def align(theta):
    best = None
    for s in (1, -1):
        dif = (s * theta - ang + 180) % 360 - 180; c = np.angle(np.mean(np.exp(1j * np.radians(dif))), deg=True)
        th = (s * theta - c) % 360; e = np.mean(np.abs((th - ang + 180) % 360 - 180))
        if best is None or e < best[0]: best = (e, th)
    return best[1]
U = umap.UMAP(n_neighbors=15, random_state=0).fit_transform(X); U = U - U.mean(0)
ax.plot(ang, align(np.degrees(np.arctan2(U[:, 1], U[:, 0]))), '.', ms=2, color='tab:orange', label='UMAP (%.1f°)' % o['err']['UMAP'])
ax.plot(ang, align(o['theta']), '.', ms=2, color='k', label='arc length (%.2f°)' % o['err']['ours (arc length)'])
ax.set_xlabel('true angle'); ax.set_ylabel('recovered angle'); ax.legend(fontsize=7, markerscale=4, loc='upper left'); ax.set_title('(d) recovered rotation angle', fontsize=8.5)
plt.savefig(REPO + '/book/figures/paper/images.pdf', dpi=200); plt.close()
