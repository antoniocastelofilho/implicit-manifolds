import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np
exec(open('fig_s2s1.py').read().split("# sections of the fine approximation")[0].split("# global view")[0])
from polfile import read_pol
from render import slice3
n,k,D,cells=read_pol(REPO + '/book/runs/simplicial/s2s1_cs.pol')
X0=[c.coords.copy() for c in cells]
def chi(P):
    V={};E=set()
    ids=[[V.setdefault(tuple(np.round(q,7)),len(V)) for q in p] for p in P]
    for I in ids:
        for a,b in zip(I,I[1:]+I[:1]):
            if a!=b: E.add(frozenset((a,b)))
    return len(V)-len(E)+len(P)
for name,ax_,c,w in [('z=0',2,0.0,None),('u=0.2',3,0.2,None),('z+u=1.2',2,1.2,True)]:
    sub=[]
    for cl,X in zip(cells,X0):
        h=X[:,2]+X[:,3] if w else X[:,ax_]
        if h.min()<=c<=h.max():
            if w: Y=X.copy(); Y[:,2]=h; cl.coords=Y
            sub.append(cl)
    P=slice3(sub,ax_,c)
    for cl,X in zip(cells,X0): cl.coords=X
    print(name,len(P),chi(P))
