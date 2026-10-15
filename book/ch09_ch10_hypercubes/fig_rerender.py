"""z-buffer re-renders of the hypercube figures of chapters 8 and 9"""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from polfile import read_pol, polygons
from slicer import slice_all, axis
from zbuf import render, crop
FG=REPO + '/book/figures/'; R=REPO + '/book/runs/hypercubes/'; CB=REPO + '/book/runs/fch_book/'
def polys_of(cells): return [x.coords[cyc] for x in cells for cyc in polygons(x)]
def cols(v,cmap,lo=None,hi=None):
    v=np.asarray(v); lo=v.min() if lo is None else lo; hi=v.max() if hi is None else hi
    cm=plt.get_cmap(cmap); return [cm(np.clip((x-lo)/(hi-lo+1e-12),0,1))[:3] for x in v]
def panel(ims,titles,out,figsize):
    fig,axs=plt.subplots(1,len(ims),figsize=figsize)
    for a,im,t in zip(np.atleast_1d(axs),ims,titles): a.imshow(im); a.axis('off'); a.set_title(t,fontsize=10)
    plt.subplots_adjust(0.01,0.01,0.99,0.9,0.04,0); plt.savefig(FG+out,dpi=200); plt.close()
def stereo(P,r): return r*P[:,:3]/(r-P[:,3:4]+1e-9)
if "s3only" not in sys.argv:
    # sphere: MS vs CMH
    ims=[];tt=[]
    for f,t in [('sph3_MS.pol','Marching Simplex ($K_1$)'),('sph3_GMHC.pol','Combinatorial Marching Hypercubes')]:
        n,k,D,c=read_pol(R+f); P=polys_of(c)
        ims.append(crop(render(P,cols([p[:,2].mean() for p in P],'Blues',-2,1.3),elev=20,azim=-55,size=(900,900),edge_alpha=0.6)))
        tt.append(t+f': {len(P)} cells')
    panel(ims,tt,'sphere_ms_cmh.png',(9,4.4))

# S^3: exploded cells of a few cubes and a slice
n,k,D,cc=read_pol(R+'sph4_GMHC.pol')
G={}
for c in cc: G.setdefault(tuple(c.grid),[]).append(c)
best=None
for g in G:
    blk=[tuple(np.add(g,(i,j,l,0))) for i in (0,1) for j in (0,1) for l in (0,1)]
    if all(b in G for b in blk):
        sc=sum(len(x.coords) for b in blk for x in G[b])
        if best is None or sc>best[0]: best=(sc,blk)
sel=set(best[1]); S=[x for b in best[1] for x in G[b]]
ctr=np.vstack([c.coords for c in S]).mean(0)
P=[];cv=[]
tab=plt.get_cmap('tab10')
for i,c in enumerate(S):
    g=c.coords.mean(0); sh=0.6*(g-ctr)
    for cyc in polygons(c):
        Q=c.coords[cyc]+sh; P.append(Q[:,:3]); cv.append(tab(i%10)[:3])
im1=crop(render(P,cv,elev=18,azim=-50,size=(900,900),edge_alpha=0.8))
SL=slice_all(cc,axis(4,3),0.3); Q=[p[:,:3] for p in polys_of(SL)]
im2=crop(render(Q,cols([q[:,2].mean() for q in Q],'Blues',-2.5,1.2),elev=18,azim=-50,size=(900,900),edge_alpha=0.6))
panel([im1,im2],[f'cells of {len(sel)} cubes ({len(S)} polytopes, exploded)',f'slice $x_4=0.3$ ({len(Q)} polygons)'],'s3_cmh.png',(10,4.6))
import sys as _s
if "s3only" in _s.argv: raise SystemExit
# torus in R^4 by GCCH
n,k,D,c=read_pol(R+'tor4_GCHC.pol'); P=polys_of(c)
Q=[stereo(p,np.sqrt(1.25)) for p in P]
im=crop(render(Q,cols([p[:,2].mean() for p in P],'viridis'),elev=35,azim=-55,size=(1300,950),edge_alpha=0.6))
plt.imsave(FG+'torus4_hypercube.png',im); print('torus',len(P))
# complex curves
ims=[]
for cs in ('ZCOSW','Z2COSW','Z2W5'):
    n,k,D,c=read_pol(CB+f'BOOK_{cs}/fch.pol'); P=polys_of(c)
    ims.append(crop(render([p[:,:3] for p in P],cols([p[:,3].mean() for p in P],'viridis'),elev=22,azim=-60,size=(900,900),edge_alpha=0.5)))
panel(ims,['$z=\\cos w$','$z^2=\\cos w$','$z^2=w^5+1+i$'],'complex_gcch.png',(11,4))
