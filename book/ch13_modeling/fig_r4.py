import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle, os; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from polfile import polygons, euler
from slicer import iter_pol, slice_all
from zbuf import render, crop
FG=REPO + '/book/figures/'; R=REPO + '/book/runs/modeling/'
def cols(v,cmap,lo,hi):
    cm=plt.get_cmap(cmap); return [cm(np.clip((x-lo)/(hi-lo+1e-12),0,1))[:3] for x in v]
def sections(path, cuts):
    """cuts: list of (axis index, value); returns list of sliced cell lists, total cells"""
    keep=[[] for _ in cuts]; tot=0
    for c in iter_pol(path):
        tot+=1
        for j,(i,v) in enumerate(cuts):
            if c.coords[:,i].min()<=v<=c.coords[:,i].max(): keep[j].append(c)
    out=[]
    for (i,v),K in zip(cuts,keep):
        a=np.zeros(4); a[i]=1; out.append(slice_all(K,a,v))
    return tot,out
def drop(P,i): return P[:,[j for j in range(4) if j!=i]]
def figure(case, cuts, out, cmap, elev=25, azim=-60, alg='fch', figsize=None, cv=2):
    tot,S=sections(R+case+f'/{alg}.pol', cuts)
    ims=[]; tt=[]; res=[]
    for (i,v),cells in zip(cuts,S):
        P=[drop(x.coords[cyc],i) for x in cells for cyc in polygons(x)]
        chi=euler(cells)[1] if cells else 0; res.append((len(P),chi))
        if P:
            Z=np.vstack(P)[:,cv]
            ims.append(crop(render(P,cols([p[:,cv].mean() for p in P],cmap,Z.min(),Z.max()),elev=elev,azim=azim,size=(900,800),edge_alpha=0.45)))
        else: ims.append(np.ones((10,10,3)))
        tt.append(f'$x_{i+1}={v}$: {len(P)} polygons, $\\chi={chi}$')
    fig,axs=plt.subplots(1,len(ims),figsize=figsize or (4*len(ims),3.9))
    for a,im,t in zip(axs,ims,tt): a.imshow(im); a.axis('off'); a.set_title(t,fontsize=9)
    plt.subplots_adjust(0.01,0.01,0.99,0.9,0.03,0); plt.savefig(FG+out,dpi=200); plt.close()
    print(case,alg,tot,res,flush=True); return tot,res
info={}
w=sys.argv[1:]
if 'balls' in w:
    info['balls']=figure('MOD_R4_TWO_BALLS',[(3,0.0),(3,0.6),(3,0.75),(3,0.86)],'mod_r4_balls.png','GnBu',elev=15,azim=-5,figsize=(14,3.6))
if 'diff' in w:
    info['diff']=figure('MOD_R4_DIFF',[(3,0.0),(3,0.7),(0,0.0),(0,0.6)],'mod_r4_diff.png','OrRd',elev=50,azim=-60,figsize=(14,3.7))
if 'tube' in w:
    info['tube']=figure('MOD_R4_TUBE',[(3,0.0),(3,0.2),(3,0.35),(3,0.62)],'mod_r4_tube.png','viridis',elev=40,azim=-60,figsize=(14,3.9))
if 'ms' in w:   # Euler characteristics of the same sections computed from the simplicial approximations
    for case,cuts in [('MOD_R4_TWO_BALLS',[(3,0.0),(3,0.6),(3,0.75),(3,0.86)]),('MOD_R4_DIFF',[(3,0.0),(3,0.7),(0,0.0),(0,0.6)]),('MOD_R4_TUBE',[(3,0.0),(3,0.3),(3,0.45),(3,0.6)])]:
        if not os.path.exists(R+case+'/ms.pol'): continue
        tot,S=sections(R+case+'/ms.pol',cuts); info[case+'_ms']=(tot,[(sum(len(polygons(x)) for x in s),euler(s)[1] if s else 0) for s in S]); print(case,'ms',info[case+'_ms'],flush=True)
pickle.dump(info,open(R+'r4_'+'_'.join(w)+'.pkl','wb'))
