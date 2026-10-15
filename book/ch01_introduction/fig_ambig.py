import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from polfile import read_pol, polygons
OUT=REPO + '/book/figures/'
V=np.array([[(b>>0)&1,(b>>1)&1,(b>>2)&1] for b in range(8)],float)
E=[(a,b) for a in range(8) for b in range(a+1,8) if bin(a^b).count('1')==1]
def panel(ax,fn,pos,title):
    n,k,D,cells=read_pol(fn)
    P=[]
    for c in cells:
        for cyc in polygons(c): P.append(c.coords[cyc])
    for a,b in E: ax.plot(*zip(V[a],V[b]),color='k',lw=1.4)
    ax.plot(*zip(V[0],V[7]),color='0.55',lw=0.9,ls='--')
    pc=Poly3DCollection(P,facecolors=(0.55,0.62,0.78,0.9),edgecolors=(0.15,0.15,0.2,0.8),linewidths=0.7)
    ax.add_collection3d(pc)
    for b in range(8):
        ax.scatter(*V[b],s=70 if b in pos else 18,color='k' if b in pos else 'w',edgecolors='k',depthshade=False,zorder=5)
    ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_zlim(0,1);ax.set_box_aspect((1,1,1))
    ax.view_init(elev=18,azim=-62);ax.axis('off');ax.set_title(title,fontsize=10)
    return len(P)
fig=plt.figure(figsize=(8.4,4))
a1=fig.add_subplot(121,projection='3d'); a2=fig.add_subplot(122,projection='3d')
print(panel(a1,REPO + '/book/runs/simplicial/ambig_split.pol',(1,6),'(a) two pieces'))
print(panel(a2,REPO + '/book/runs/simplicial/ambig_tube.pol',(0,7),'(b) a tube'))
plt.subplots_adjust(0,0,1,0.92,0.02,0)
plt.savefig(OUT+'intro_ambiguity.png',dpi=220)
