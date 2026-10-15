import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
from render import *
fig = plt.figure(figsize=(9,4.4))
for i,(f,t) in enumerate([('sph3_MS.pol','Marching Simplex ($K_1$)'),('sph3_GMHC.pol','Combinatorial Marching Hypercubes')]):
    n,k,D,cells = read_pol(REPO + '/book/runs/hypercubes/'+f)
    P = polygons(cells)
    ax = new_ax(fig, 121+i, elev=20, azim=-55)
    draw_polys(ax, P, axes_proj(3), color_values=[p[:,2].mean() for p in P], cmap='Blues', vmin=-2, vmax=1.3, lw=0.35)
    ax.set_title(t+f'\n{len(P)} cells', fontsize=10)
plt.subplots_adjust(0,0,1,0.9,0,0)
plt.savefig(REPO + '/book/figures/sphere_ms_cmh.png', dpi=220)
