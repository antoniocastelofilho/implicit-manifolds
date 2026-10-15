import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from polfile import read_pol, polygons
plt.rcParams.update({'font.size':9,'mathtext.fontset':'stix','font.family':'STIXGeneral'})
cols={0:np.array([0.80,0.82,0.86]),1:np.array([0.55,0.68,0.85]),2:np.array([0.30,0.50,0.80]),3:np.array([0.85,0.35,0.30])}
L=np.array([0.3,-0.5,0.8]);L/=np.linalg.norm(L)
fig=plt.figure(figsize=(10,4.2))
for i,lev in enumerate([1,3]):
    n,k,D,cells=read_pol(REPO + f'/book/runs/j1a/torus_L{lev}.pol')
    P=[];C=[]
    for c in cells:
        lv=c.labels[0][0]
        for cyc in polygons(c):
            Q=c.coords[cyc]; P.append(Q)
            nrm=np.cross(Q[1]-Q[0],Q[2]-Q[0]); sh=0.55+0.45*abs(nrm@L)/(np.linalg.norm(nrm)+1e-15)
            C.append(tuple(np.clip(cols[lv]*sh,0,1)))
    ax=fig.add_subplot(1,2,i+1,projection='3d')
    ax.add_collection3d(Poly3DCollection(P,facecolors=C,edgecolors=(0,0,0,0.25),linewidths=0.08))
    ax.set_xlim(-1.15,1.15);ax.set_ylim(-1.15,1.15);ax.set_zlim(-0.55,0.55);ax.set_box_aspect((2.8,2.8,1.4))
    ax.view_init(elev=38,azim=-60);ax.axis('off')
    ax.set_title(f'({"ab"[i]}) {lev} level{"s" if lev>1 else ""} of refinement, {len(cells)} cells',fontsize=10)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=cols[l],label=f'cells in {l}-blocks') for l in range(4)],loc='lower center',ncol=4,frameon=False)
plt.subplots_adjust(0,0.07,1,0.95,0,0)
plt.savefig(REPO + '/book/figures/j1a_torus.png',dpi=220)
