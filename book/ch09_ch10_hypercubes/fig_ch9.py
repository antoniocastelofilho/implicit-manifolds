import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from slicer import iter_pol, slice_all, axis
from polfile import polygons, read_pol
from zbuf import render, crop
OUT=REPO + '/book/figures/'; R=REPO + '/book/runs/fch_book/'
def polys_of(cells): return [x.coords[cyc] for x in cells for cyc in polygons(x)]
def section(path, a, c):
    a=np.asarray(a,float); keep=[]; tot=0
    for x in iter_pol(path):
        tot+=1; h=x.coords@a
        if h.min()<=c<=h.max(): keep.append(x)
    return tot, polys_of(slice_all(keep,a,c))
def stereo4(P, rad):
    # stereographic projection of S^3(rad) from (0,0,0,rad)
    P=np.asarray(P); return rad*P[:,:3]/(rad-P[:,3:4]+1e-9)
def colors(vals, cmap, lo=None, hi=None):
    v=np.asarray(vals); lo=v.min() if lo is None else lo; hi=v.max() if hi is None else hi
    cm=plt.get_cmap(cmap); return [cm((x-lo)/(hi-lo+1e-12))[:3] for x in v]
def save(img, name):
    plt.imsave(OUT+name, crop(img)); print('saved', name)
info={}
# ---------- T3 = S1xS1xS1 in R^6, section x6 = 0.9: two Clifford tori (one for each sign of x5)
tot,P=section(R+'BOOK_T3_R6/fch.pol', axis(6,5), 0.9)
rad=np.sqrt(1.25); ims=[]
for s in (-1,1):
    PP=[p for p in P if np.sign(p[:,4].mean())==s]
    Q=[stereo4(p[:,:4],rad) for p in PP]
    cols=[plt.get_cmap('hsv')((np.arctan2(p[:,1].mean(),p[:,0].mean())+np.pi)/(2*np.pi))[:3] for p in PP]
    ims.append(crop(render(Q,cols,elev=32,azim=-60,size=(1000,800),edge_alpha=0.5)))
fig,axs=plt.subplots(1,2,figsize=(10,4))
for a_,im,t in zip(axs,ims,['$x_5<0$','$x_5>0$']): a_.imshow(im); a_.axis('off'); a_.set_title(t,fontsize=12)
plt.subplots_adjust(0.01,0.01,0.99,0.92,0.05,0); plt.savefig(OUT+'hc_t3_section.png',dpi=200); plt.close()
info['T3']=(tot,len(P))
# ---------- torus in R^10, stereographic view of coordinates 1..4
n,k,D,cells=read_pol(R+'BOOK_TORUS_R10/fch.pol'); P=polys_of(cells)
Q=[stereo4(p[:,:4],np.sqrt(2)) for p in P]
cols=colors([np.arctan2(p[:,1].mean(),p[:,0].mean()) for p in P],'hsv')
save(render(Q,cols,elev=40,azim=-60,size=(1300,950),edge_alpha=0.5),'hc_torus_r10.png'); info['T10']=(len(cells),len(P))
# ---------- S2 in R^10: three coordinate projections
n,k,D,cells=read_pol(R+'BOOK_S2_R10/fch.pol'); P=polys_of(cells); info['S2']=(len(cells),len(P))
views=[((0,1,2),9,'x1x2x3'),((0,1,3),5,'x1x2x4'),((3,4,5),0,'x4x5x6')]
imgs=[]
for ax3,cc,nm in views:
    Q=[p[:,list(ax3)] for p in P]; cols=colors([p[:,cc].mean() for p in P],'viridis')
    imgs.append(crop(render(Q,cols,elev=25,azim=-55,size=(900,900),edge_alpha=0.45)))
fig,axs=plt.subplots(1,3,figsize=(12,4.3))
for a,im,(ax3,cc,nm) in zip(axs,imgs,views):
    a.imshow(im); a.axis('off'); a.set_title('$(x_{%d},x_{%d},x_{%d})$, colour $x_{%d}$'%(ax3[0]+1,ax3[1]+1,ax3[2]+1,cc+1),fontsize=11)
plt.subplots_adjust(0.01,0.01,0.99,0.92,0.03,0); plt.savefig(OUT+'hc_s2_r10.png',dpi=200); plt.close()
# ---------- S2 x S1 in R^5, the three sections of the Introduction
path=R+'BOOK_S2XS1_R5/fch.pol'
def e(i): v=np.zeros(5); v[i]=1; return v
panels=[('torus: $z=0$',e(2),0.0,lambda P: np.c_[(1+P[:,3])*P[:,0],(1+P[:,3])*P[:,1],P[:,4]]),
        ('two spheres: $u=0.2$',e(3),0.2,lambda P: np.c_[P[:,0]+3*P[:,4],P[:,1],P[:,2]]),
        ('sphere: $z+u=1.2$',e(2)+e(3),1.2,lambda P: np.c_[P[:,0],P[:,1],2*P[:,4]])]
imgs=[]; ns=[]
for t,a,c,m in panels:
    tot,P=section(path,a,c); Q=[m(p) for p in P]; ns.append(len(P))
    cols=colors([q[:,2].mean() for q in Q],'viridis')
    imgs.append(crop(render(Q,cols,elev=28,azim=-58,size=(1000,900),edge_alpha=0.45)))
info['S2S1']=(tot,ns)
fig,axs=plt.subplots(1,3,figsize=(12,4.2))
for a,im,(t,*_),nn in zip(axs,imgs,panels,ns):
    a.imshow(im); a.axis('off'); a.set_title(t+f' ({nn} polygons)',fontsize=11)
plt.subplots_adjust(0.01,0.01,0.99,0.92,0.03,0); plt.savefig(OUT+'hc_s2s1_sections.png',dpi=200); plt.close()
# ---------- Klein bottle (article example), z-buffer re-render
n,k,D,cells=read_pol(REPO + '/book/runs/fch/KLEIN_5TO3/fch.pol'); P=polys_of(cells)
cols=colors([(p[:,3].mean()/(2*np.pi))%1 for p in P],'twilight',0,1)
save(render([p[:,:3] for p in P],cols,elev=30,azim=-60,size=(1500,1050),edge_alpha=0.5),'klein.png'); info['Klein']=(len(cells),len(P))
print(info); pickle.dump(info,open(REPO + '/book/ch09_ch10_hypercubes/ch9_info.pkl','wb'))
