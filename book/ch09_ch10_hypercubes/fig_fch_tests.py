import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np
from slicer import iter_pol, slice_all, axis
from polfile import polygons
# S4 in R5 (20^5), S5 in R7 (10^7), T4 in R8 (5^8): sections down to surfaces
jobs={
 'S4': (REPO + '/book/runs/fch/HYPERSHPERE_5TO1/fch.pol', [(3,0.3),(4,0.2)]),
 'S5': (REPO + '/book/runs/fch/HYPERSHPERE_N_TO_K/fch.pol', [(3,0.2),(4,0.3),(5,0.1)]),
 'T4': (REPO + '/book/runs/fch/S1_K_TIMES/fch.pol', [(5,0.6),(7,0.6)]),
}
out={}
for name,(f,cuts) in jobs.items():
    keep=[]; total=0
    for c in iter_pol(f):
        total+=1
        X=c.coords
        if all(X[:,i].min()<=v<=X[:,i].max() for i,v in cuts): keep.append(c)
    n=keep[0].coords.shape[1] if keep else 0
    S=keep
    for i,v in cuts: S=slice_all(S,axis(n,i),v)
    P=[p for x in S for p in polygons(x)]
    Q=[x.coords[cyc] for x in S for cyc in polygons(x)]
    out[name]=dict(total=total,kept=len(keep),polys=Q)
    print(name,total,len(keep),len(Q),flush=True)
pickle.dump(out,open(REPO + '/book/ch09_ch10_hypercubes/fch_tests.pkl','wb'))
