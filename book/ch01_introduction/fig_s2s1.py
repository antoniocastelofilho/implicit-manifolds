import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from polfile import read_pol
from render import slice3, draw_polys, draw_segments, segments, new_ax
OUT=REPO + '/book/figures/'
# global view: coarse approximation, edges, through (x,y,z,u,v) -> ((2+x)2u,(2+x)2v,y)
n,k,D,cc=read_pol(REPO + '/book/runs/simplicial/s2s1_coarse.pol')
S=segments(cc)
phi=lambda P: np.c_[(2+P[:,0])*2*P[:,3],(2+P[:,0])*2*P[:,4],P[:,1]]
fig=plt.figure(figsize=(5.2,4.6))
ax=new_ax(fig,111,elev=35,azim=-60)
segs=[phi(np.asarray(s)) for s in S]
ang=[np.arctan2(np.mean(np.asarray(s)[:,4]),np.mean(np.asarray(s)[:,3])) for s in S]
draw_segments(ax,segs,np.eye(3),lw=0.12,color_values=ang,cmap='hsv')
plt.subplots_adjust(0,0,1,1); plt.savefig(OUT+'intro_s2s1_global.png',dpi=220); plt.close()
print('global',len(S))
del cc
# sections of the fine approximation
n,k,D,cells=read_pol(REPO + '/book/runs/simplicial/s2s1_cs.pol')
X0=[c.coords.copy() for c in cells]
def section(w, c):
    """slice by w.x=c: replace coordinate 2 (z) by w.x, slice, and restore z"""
    sub=[]
    for cl,X in zip(cells,X0):
        h=X@w
        if h.min()<=c<=h.max():
            Y=X.copy(); Y[:,2]=h; cl.coords=Y; sub.append(cl)
    P=slice3(sub,2,c)
    for cl,X in zip(cells,X0): cl.coords=X
    out=[]
    for p in P:
        q=p.copy(); zz=w[2]
        q[:,2]=(c-(p@w-p[:,2]*w[2]))/w[2] if w[2]!=0 else p[:,2]
        out.append(q)
    return out
def ez(i): e=np.zeros(5); e[i]=1; return e
panels=[
 ('torus: $z=0$', ez(2), 0.0, lambda P: np.c_[(1+P[:,3])*P[:,0],(1+P[:,3])*P[:,1],P[:,4]]),
 ('two spheres: $u=0.2$', ez(3), 0.2, lambda P: np.c_[P[:,0]+3*P[:,4],P[:,1],P[:,2]]),
 ('sphere: $z+u=1.2$', ez(2)+ez(3), 1.2, lambda P: np.c_[P[:,0],P[:,1],2*P[:,4]]),
]
fig=plt.figure(figsize=(12,4.2))
for i,(t,w,c,m) in enumerate(panels):
    if w[2]==0:   # hyperplane u=c: slice along coordinate 3 directly
        P=slice3([cl for cl in cells if cl.coords[:,3].min()<=c<=cl.coords[:,3].max()],3,c)
    else:
        P=section(w,c)
    ax=new_ax(fig,131+i,elev=28,azim=-58)
    Q=[m(p) for p in P]
    cv=[q[:,2].mean() for q in Q]
    draw_polys(ax,Q,np.eye(3),color_values=cv,cmap='viridis',lw=0.12,zoom=1.25)
    ax.set_title(t+f' ({len(P)} polygons)',fontsize=10)
    print(t,len(P))
plt.subplots_adjust(0,0,1,0.93,0,0)
plt.savefig(OUT+'intro_s2s1_sections.png',dpi=200)
