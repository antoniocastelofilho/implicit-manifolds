import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
F = lambda x,y: np.sin(2.2*x+0.3)*np.sin(2.2*y+0.2) + 0.06 - 0.25*(x*x+y*y)
N=8; xs=np.linspace(-1.6,1.6,N+1)+0.013; ys=np.linspace(-1.2,1.2,7)+0.007
def zero(p,q,fp,fq): t=fp/(fp-fq); return p+t*(q-p)
fig,axs=plt.subplots(1,2,figsize=(10,3.3))
for ax in axs:
    for x in xs: ax.plot([x,x],[ys[0],ys[-1]],color='0.75',lw=0.6)
    for y in ys: ax.plot([xs[0],xs[-1]],[y,y],color='0.75',lw=0.6)
    ax.set_aspect('equal'); ax.axis('off')
amb=[]
for i in range(N):
    for j in range(len(ys)-1):
        P=[np.array([xs[i],ys[j]]),np.array([xs[i+1],ys[j]]),np.array([xs[i+1],ys[j+1]]),np.array([xs[i],ys[j+1]])]
        V=[F(*p) for p in P]
        # K1 triangles: v0,v1,v2 and v0,v3,v2 (diagonal v0-v2)
        axs[0].plot([P[0][0],P[2][0]],[P[0][1],P[2][1]],color='0.85',lw=0.5)
        for T in ([0,1,2],[0,3,2]):
            pts=[]
            for u,v in ((T[0],T[1]),(T[1],T[2]),(T[0],T[2])):
                if V[u]*V[v]<0: pts.append(zero(P[u],P[v],V[u],V[v]))
            if len(pts)==2:
                axs[0].plot(*np.array(pts).T,color='C3',lw=1.6); axs[0].plot(*np.array(pts).T,'o',color='C3',ms=2.5)
        # hypercube: crossings on the 4 edges, connected via traced arcs (pairs)
        cr=[]
        for e,(u,v) in enumerate(((0,1),(1,2),(3,2),(0,3))):
            if V[u]*V[v]<0: cr.append((e,zero(P[u],P[v],V[u],V[v])))
        for e,p in cr: axs[1].plot(*p,'o',color='C0',ms=3.5)
        if len(cr)==2: axs[1].plot([cr[0][1][0],cr[1][1][0]],[cr[0][1][1],cr[1][1][1]],color='C0',lw=1.6)
        elif len(cr)==4:
            amb.append((i,j)); axs[1].fill([p[0] for p in P],[p[1] for p in P],color='C1',alpha=0.25)
            # K1 resolution: diagonal v0-v2 separates; pair crossings adjacent to the vertex whose sign differs from the diagonal's
            s=np.sign(V)
            # the sign at v1 and v3 (off-diagonal) are equal; if equal to F_T at diagonal crossing... use PL zero set on triangles
            segs=[]
            for T in ([0,1,2],[0,3,2]):
                pts=[]
                for u,v in ((T[0],T[1]),(T[1],T[2]),(T[0],T[2])):
                    if V[u]*V[v]<0: pts.append(zero(P[u],P[v],V[u],V[v]))
                if len(pts)==2: segs.append(pts)
            for sg in segs: axs[1].plot(*np.array(sg).T,color='C0',lw=1.6)
axs[0].set_title('simplicial: zero set of $F_\\mathcal{T}$ on the $K_1$ triangles',fontsize=9)
axs[1].set_title('hypercube: vertices on grid edges, one arc per square',fontsize=9)
plt.tight_layout(); plt.savefig(REPO + '/book/figures/gamma_2d.pdf'); print('ambiguous',amb)
