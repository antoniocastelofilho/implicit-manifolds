import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
from render import *
base=REPO + '/python/fch_lib/out/'
# Cahn-Hilliard 100D
n,k,D,cells=read_pol(base+'hilliard100/results.pol')
P=polygons(cells)
fig=plt.figure(figsize=(10,4.6))
for j,(idx,el,az,ttl) in enumerate([((0,2,4),20,-60,r'$(\lambda,a_0,a_2)$'),((0,1,4),20,-50,r'$(\lambda,c,a_2)$')]):
    ax=new_ax(fig,121+j,elev=el,azim=az)
    pr=axes_proj(n,idx)
    # rescale axes for visibility
    S=np.diag([1/4.0,1/0.6,1/0.27]) if idx==(0,2,4) else np.diag([1/4.0,1/0.55,1/0.27])
    draw_polys(ax,P,S@pr,color_values=[p[:,1 if idx==(0,2,4) else 2].mean() for p in P],cmap='viridis',lw=0.15,zoom=1.25)
    ax.set_title(ttl,fontsize=10)
plt.subplots_adjust(0,0,1,0.95,0,0)
plt.savefig(REPO + '/book/figures/cahn_hilliard100.png',dpi=200)
