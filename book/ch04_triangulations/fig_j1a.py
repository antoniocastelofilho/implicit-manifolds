import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds/j1a')
from j1a import *
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'mathtext.fontset':'stix','font.family':'STIXGeneral'})
from matplotlib.collections import PolyCollection, LineCollection
F=lambda p: (p[0]-4.05)**2/6.25+(p[1]-3.95)**2/4.0-1
Gf=lambda p: p[1]-4.6-0.15*(p[0]-4)
def simps(G,b):
    out=[]
    for P in simplices(G,b):
        out.append([tuple(float(x) for x in p) for p in P])
    return out
def is_transition(G,b): return len(refined_faces(G,b))>0
def crossN(P):
    A=np.array([[1,1,1],[F(p) for p in P],[Gf(p) for p in P]]); 
    try: l=np.linalg.solve(A,[1,0,0])
    except np.linalg.LinAlgError: return False
    return (l>=0).all()
def seg(P):
    v=[F(p) for p in P]; pts=[]
    for a,b in ((0,1),(1,2),(0,2)):
        if v[a]*v[b]<0:
            t=v[a]/(v[a]-v[b]); pts.append((P[a][0]+t*(P[b][0]-P[a][0]),P[a][1]+t*(P[b][1]-P[a][1])))
    return pts if len(pts)==2 else None
G=Grid(2,4)
for it in range(3):
    mark=[b for b in G.leaves if b[0]==it and any(crossN(P) for P in simps(G,b))]
    for b in mark:
        if b in G.leaves: G.refine(b)
print(sorted(Counter(b[0] for b in G.leaves).items()))
fig=plt.figure(figsize=(11,4.4))
# panels a-c: single blocks
def block_demo(ax,Phi,title):
    tris=[]
    for s in itertools.product((-1,1),repeat=2):
        for pi in itertools.permutations(range(2)):
            steps=[(pi[j],s[pi[j]]) for j in range(2)]
            h=None
            for hh in (1,2):
                if any(phi<=frozenset(steps[:hh]) for phi in Phi): h=hh;break
            if h is None or h==2:
                v=[0,0];P=[tuple(v)]
                for j in range(2): v[pi[j]]+=s[pi[j]];P.append(tuple(v))
                tris.append((P,False));continue
            for t in (-1,1):
                v=[0.0,0.0];P=[tuple(v)];v[pi[0]]+=s[pi[0]];v[pi[1]]+=s[pi[1]]/2;P.append(tuple(v));v[pi[1]]+=t/2;P.append(tuple(v));tris.append((P,True))
    ax.add_collection(PolyCollection([P for P,_ in tris],facecolors=['0.82' if tr else 'white' for _,tr in tris],edgecolors='k',linewidths=0.8))
    ax.plot([0],[0],'ko',ms=3)
    ax.set_xlim(-1.15,1.15);ax.set_ylim(-1.15,1.15);ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=9)
gs=fig.add_gridspec(2,4)
block_demo(fig.add_subplot(gs[0,0]),[],'(a) basic block')
block_demo(fig.add_subplot(gs[0,1]),[frozenset({(0,1)})],'(b) refined face $\\{1\\}$')
block_demo(fig.add_subplot(gs[1,0]),[frozenset({(0,1)}),frozenset({(1,1)})],'(c) refined faces $\\{1\\},\\{2\\}$')
block_demo(fig.add_subplot(gs[1,1]),[frozenset({(0,1)}),frozenset({(1,-1)}),frozenset({(0,-1)})],'(d) three refined faces')
ax=fig.add_subplot(gs[:,2:])
polys=[];fc=[];segs=[]
for b in G.leaves:
    tr=is_transition(G,b)
    for P in simps(G,b):
        polys.append(P); fc.append('0.85' if tr else 'white')
        s=seg(P)
        if s: segs.append(s)
ax.add_collection(PolyCollection(polys,facecolors=fc,edgecolors='0.35',linewidths=0.3))
ax.add_collection(LineCollection(segs,colors='C3',linewidths=1.6))
xs=np.linspace(0,8,4001); 
# points of N (exact)
from scipy.optimize import brentq
for x0 in (1.0,7.0):
    try:
        xr=brentq(lambda x: F((x,4.6+0.15*(x-4))),x0,4.05); ax.plot(xr,4.6+0.15*(xr-4),'ko',ms=4)
    except Exception as e: print(e)
for b in G.leaves:
    c,hs=G.box(b); c=[float(v) for v in c]; hs=float(hs)
    ax.add_patch(plt.Rectangle((c[0]-hs,c[1]-hs),2*hs,2*hs,fill=False,lw=0.9,ec='k'))
ax.set_xlim(0,8);ax.set_ylim(0,8);ax.set_aspect('equal');ax.axis('off')
ax.set_title('(e)',fontsize=11)
plt.tight_layout()
plt.savefig(REPO + '/book/figures/j1a_2d.pdf'); print(len(polys),len(segs))
