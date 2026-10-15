import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from polfile import read_pol, polygons
OUT=REPO + '/book/figures/'
lo,h=-1.1,2.2/12
def sel(fn):
    n,k,D,cells=read_pol(fn); P=[];cubes=set()
    for c in cells:
        g=tuple(np.floor((c.coords.mean(0)-lo)/h).astype(int))
        if 10<=g[0]<=11 and 4<=g[1]<=7 and 4<=g[2]<=7:
            cubes.add(g)
            for cyc in polygons(c): P.append(c.coords[cyc])
    return P,cubes
def draw(ax,P,cubes,title,col):
    for g in cubes:
        o=lo+h*np.array(g)
        for a in range(8):
            for b in range(a+1,8):
                if bin(a^b).count('1')==1:
                    A=o+h*np.array([(a>>i)&1 for i in range(3)]);B=o+h*np.array([(b>>i)&1 for i in range(3)])
                    ax.plot(*zip(A,B),color='0.6',lw=0.5)
    L=np.array([0.5,0.3,0.8]);L/=np.linalg.norm(L)
    fc=[]
    for p in P:
        nrm=np.cross(p[1]-p[0],p[2]-p[0]); s=0.55+0.45*abs(nrm@L)/(np.linalg.norm(nrm)+1e-15)
        fc.append(tuple(np.array(col)*s)+(0.95,))
    ax.add_collection3d(Poly3DCollection(P,facecolors=fc,edgecolors=(0.1,0.1,0.1,0.9),linewidths=0.6))
    ax.set_xlim(lo+10*h,lo+12*h);ax.set_ylim(lo+4*h,lo+8*h);ax.set_zlim(lo+4*h,lo+8*h)
    ax.set_box_aspect((2,4,4));ax.view_init(elev=15,azim=-25);ax.axis('off')
    ax.set_title(title+f' ({len(P)} cells)',fontsize=10)
fig=plt.figure(figsize=(9,4.6))
P,c=sel(REPO + '/book/runs/simplicial/sph3_ms.pol'); draw(fig.add_subplot(121,projection='3d'),P,c,'(a) Marching Simplex',(0.85,0.45,0.35)); print(len(P),len(c))
P,c=sel(REPO + '/book/runs/hypercubes/sph3_MHC.pol'); draw(fig.add_subplot(122,projection='3d'),P,c,'(b) Combinatorial Marching Hypercubes',(0.35,0.5,0.85)); print(len(P),len(c))
plt.subplots_adjust(0,0,1,0.93,0,0)
plt.savefig(OUT+'intro_cells.png',dpi=220)
