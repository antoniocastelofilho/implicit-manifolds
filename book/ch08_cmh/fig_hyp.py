import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from polfile import read_pol, polygons, euler
from slicer import slice_all, axis
from zbuf import render, crop
FG=REPO + '/book/figures/'
def components(cells, nd=6):
    key={}; par=[]
    def f(i):
        while par[i]!=i: par[i]=par[par[i]]; i=par[i]
        return i
    for c in cells:
        ids=[]
        for p in c.coords:
            k=tuple(np.round(p,nd))
            if k not in key: key[k]=len(par); par.append(len(par))
            ids.append(key[k])
        for a in ids[1:]: par[f(a)]=f(ids[0])
    return len({f(i) for i in range(len(par))})
ims=[];tt=[]; res={}
for f,t in [(REPO + '/book/runs/hypercubes/hyp_MHC.pol','CMH, positive corners separated'),(REPO + '/book/runs/modeling/BOOK_HYPERBOLA4/ms.pol','Freudenthal–Kuhn (Marching Simplex)')]:
    n,k,D,cells=read_pol(f)
    comp=components(cells); chi=None
    S=slice_all([c for c in cells if c.coords[:,3].min()<=0.35<=c.coords[:,3].max()],axis(4,3),0.35)
    P=[x.coords[cyc][:,:3] for x in S for cyc in polygons(x)]
    sc=components(S); se=euler(S)[1]
    res[f]=(len(cells),comp,len(P),sc,se)
    # colour the section polygons by connected component
    key={}; par=[]
    def fd(i):
        while par[i]!=i: par[i]=par[par[i]]; i=par[i]
        return i
    pid=[]
    for p in P:
        ids=[]
        for q in p:
            kk=tuple(np.round(q,6))
            if kk not in key: key[kk]=len(par); par.append(len(par))
            ids.append(key[kk])
        for a in ids[1:]: par[fd(a)]=fd(ids[0])
        pid.append(ids[0])
    roots={}; cmap=plt.get_cmap('tab10')
    cols=[cmap(roots.setdefault(fd(i),len(roots))%10)[:3] for i in pid]
    ims.append(crop(render(P,cols,elev=25,azim=-50,size=(1000,1000),edge_alpha=0.5)))
    tt.append(f'{t}: {len(cells)} cells, {comp} component{'s' if comp>1 else ''}\nsection $x_4=0.35$: {sc} component{'s' if sc>1 else ''}, $\\chi={se}$')
fig,axs=plt.subplots(1,2,figsize=(10,5.2))
for a,im,t in zip(axs,ims,tt): a.imshow(im); a.axis('off'); a.set_title(t,fontsize=9)
plt.subplots_adjust(0.01,0.01,0.99,0.88,0.03,0); plt.savefig(FG+'cmh_hyperbola.png',dpi=200); print(res)
