import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0, REPO + '/python/implicit_manifolds'); sys.path.insert(0, REPO + '/book/ch12_comparison')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from polfile import read_pol, polygons
from zbuf import render, crop
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix
from bratu import natural, arclength
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
B = REPO + '/book/runs/bench/'
def polys_comp(path):
    _, _, _, cells = read_pol(path)
    P = []; key = {}; owner = []
    for ci, c in enumerate(cells):
        for cyc in polygons(c):
            P.append(c.coords[cyc]); owner.append(ci)
    # components by shared rounded vertices
    ids = []; rows = []; cols = []
    for i, p in enumerate(P):
        for v in np.round(p, 7):
            t = tuple(v); j = key.setdefault(t, len(key)); rows.append(i); cols.append(j)
    A = coo_matrix((np.ones(len(rows)), (rows, cols))).tocsr(); G = A @ A.T
    nc, lab = connected_components(G, directed=False)
    return P, lab, nc
Pm, lm, nm = polys_comp(B + 'BOOK_BLOBS3/ms.pol'); Pc, lc, ncc = polys_comp(B + 'BOOK_BLOBS3/cs.pol')
print('components', nm, ncc)
cm = plt.get_cmap('tab20')
colm = [cm(l % 20)[:3] for l in lm]
# continuation: the traced component in colour, the others grey (not found)
cskeys = set(tuple(v) for p in Pc for v in np.round(p, 7))
colc = []
for p, l in zip(Pm, lm):
    colc.append((0.85, 0.85, 0.85))
traced = set(l for p, l in zip(Pm, lm) if any(tuple(v) in cskeys for v in np.round(p, 7)))
colc = [cm(l % 20)[:3] if l in traced else (0.88, 0.88, 0.88) for l in lm]
im1 = crop(render(Pm, colm, elev=28, azim=-58, size=(900, 800), edge_alpha=0.25))
im2 = crop(render(Pm, colc, elev=28, azim=-58, size=(900, 800), edge_alpha=0.25))
fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
axs[0].imshow(im1); axs[0].set_title('(a) Marching Simplex: %d components' % nm)
axs[1].imshow(im2); axs[1].set_title('(b) Continuation Simplex from one point: 1 component')
for a in axs: a.axis('off')
plt.tight_layout(); plt.savefig(REPO + '/book/figures/cmp_blobs.png', dpi=200); plt.close()
# Bratu
fig, axs = plt.subplots(1, 2, figsize=(9, 3.2))
for a, N, mid in ((axs[0], 3, 2), (axs[1], 5, 3)):
    ref = arclength(N, umax=9); a.plot(ref[:, 0], ref[:, mid], '-', color='0.75', lw=3, label='pseudo-arclength (reference)')
    _, _, _, cells = read_pol(B + 'BOOK_BRATU%d/cs.pol' % N)
    for c in cells:
        a.plot(c.coords[:, 0], c.coords[:, mid], 'k-', lw=0.9)
    a.plot([], [], 'k-', lw=0.9, label='Continuation Simplex')
    nat, fail = natural(N, np.arange(0, 4.0, 0.05))
    a.plot(nat[:, 0], nat[:, mid], 'o', ms=2.5, color='tab:red', label='natural parameter (fails at $\\lambda=%.2f$)' % fail)
    a.set_xlabel('$\\lambda$'); a.set_ylabel('$u_{%d}$' % mid); a.set_title('discrete Bratu problem, $N=%d$: a curve in $\\mathbb{R}^{%d}$' % (N, N + 1))
    a.set_xlim(-0.1, 3.75); a.set_ylim(-0.2, 6.8)
axs[0].legend(fontsize=7, loc='upper right')
plt.tight_layout(); plt.savefig(REPO + '/book/figures/cmp_bratu.pdf'); plt.close()
