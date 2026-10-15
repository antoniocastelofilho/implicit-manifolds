import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle, os; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from polfile import read_pol, polygons, euler
from zbuf import render, crop
FG=REPO + '/book/figures/'; R=REPO + '/book/runs/modeling/'
info=pickle.load(open(R+'info.pkl','rb')) if os.path.exists(R+'info.pkl') else {}
def polys_of(cells): return [x.coords[cyc] for x in cells for cyc in polygons(x)]
def cols(v,cmap,lo=None,hi=None):
    v=np.asarray(v); lo=v.min() if lo is None else lo; hi=v.max() if hi is None else hi
    cm=plt.get_cmap(cmap); return [cm(np.clip((x-lo)/(hi-lo+1e-12),0,1))[:3] for x in v]
def pair(case, out, elev, azim, cmap='Blues', cv=2, lo=None, hi=None, size=(1000,900), figsize=(10,4.6), titles=('Marching Simplex ($K_1$)','Fast Continuation Hypercubes'), extra=None, ea=0.45):
    ims=[]; tt=[]
    for alg,t in zip(('ms','fch'),titles):
        n,k,D,c=read_pol(R+case+f'/{alg}.pol'); P=polys_of(c); chi=euler(c)[1]
        info[(case,alg)]=(len(c),chi)
        ims.append(crop(render([p[:,:3] for p in P],cols([p[:,cv].mean() for p in P],cmap,lo,hi),elev=elev,azim=azim,size=size,edge_alpha=ea)))
        tt.append(f'{t}: {len(c)} cells, $\\chi={chi}$')
    if extra is not None:
        for case2,alg,t in extra:
            n,k,D,c=read_pol(R+case2+f'/{alg}.pol'); P=polys_of(c); chi=euler(c)[1]; info[(case2,alg)]=(len(c),chi)
            ims.append(crop(render([p[:,:3] for p in P],cols([p[:,cv].mean() for p in P],cmap,lo,hi),elev=elev,azim=azim,size=size,edge_alpha=ea)))
            tt.append(f'{t}: {len(c)} cells, $\\chi={chi}$')
    fig,axs=plt.subplots(1,len(ims),figsize=figsize)
    for a,im,t in zip(axs,ims,tt): a.imshow(im); a.axis('off'); a.set_title(t,fontsize=9)
    plt.subplots_adjust(0.01,0.01,0.99,0.9,0.03,0); plt.savefig(FG+out,dpi=200); plt.close(); print(out,tt,flush=True)
which=sys.argv[1:]
if 'spheres' in which:
    pair('MOD_TWO_SPHERES','two_spheres.png',12,-5,cmap='YlGnBu',cv=2,lo=-1.8,hi=1.3,size=(900,800),figsize=(12,3.9),
         titles=('$\\mathcal{M}_\\cup$, Marching Simplex','$\\mathcal{M}_\\cup$, FCH'),extra=[('MOD_TWO_SPHERES_P8','fch','$\\mathcal{M}_{\\cup,8}$, FCH')])
if 'tori' in which:
    pair('MOD_EIGHT_TORI','eight_tori.png',40,-60,cmap='viridis',cv=0,size=(1100,800),figsize=(11,4.3))
if 'knot' in which:
    pair('MOD_KNOT','knot.png',35,-60,cmap='plasma',cv=2,size=(1000,900),figsize=(10,4.6))
if 'inter' in which:
    pair('MOD_INTERSECT3','intersection.png',55,-120,cmap='YlOrRd',cv=1,size=(900,900),figsize=(9,4.4))
if 'pipebox' in which:
    pair('MOD_PIPEBOX','pipebox.png',32,-55,cmap='PuBu',cv=2,lo=-9,hi=5,ea=0.18,size=(800,1100),figsize=(9,5.6),
         titles=('Marching Simplex','FCH'),extra=[('MOD_PIPEBOX_P8','fch','$p=8$, FCH')])
pickle.dump(info,open(R+'info.pkl','wb'))
