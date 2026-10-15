import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import itertools, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
OUT=REPO + '/book/figures/'
plt.rcParams.update({'font.size':9,'mathtext.fontset':'stix','font.family':'STIXGeneral'})
C=['#4C72B0','#DD8452','#55A868','#C44E52','#8172B3','#937860']
def K1(g,pi):
    v=list(g);V=[tuple(v)]
    for p in pi: v[p]+=1;V.append(tuple(v))
    return V
def J1(g,pi):
    s=[(-1)**gj for gj in g];v=[gj+(gj%2) for gj in g];V=[tuple(v)]
    for p in pi: v[p]+=s[p];V.append(tuple(v))
    return V
def grid(ax,nx,ny,kind):
    tris=[]
    for g in itertools.product(range(nx),range(ny)):
        for pi in itertools.permutations(range(2)):
            tris.append(K1(g,pi) if kind=='K' else J1(g,pi))
    ax.add_collection(PolyCollection(tris,facecolors='none',edgecolors='0.55',linewidths=0.6))
    for x in range(nx+1): ax.plot([x,x],[0,ny],color='0.25',lw=0.9)
    for y in range(ny+1): ax.plot([0,nx],[y,y],color='0.25',lw=0.9)
    ax.set_xlim(-0.15,nx+0.15);ax.set_ylim(-0.15,ny+0.15);ax.set_aspect('equal');ax.axis('off')
def lab(v,g): return sum((v[j]-g[j])<<j for j in range(len(v)))
# ---------- Fig 4.1: CFK of I^2 and I^3
fig=plt.figure(figsize=(9,3.8))
ax=fig.add_axes([0,0,0.38,0.9])
for k,pi in enumerate([(0,1),(1,0)]):
    V=K1((0,0),pi); ax.add_patch(plt.Polygon(V,fc=C[k],alpha=0.35,ec='k',lw=1.2))
    cen=np.mean(V,0); ax.text(*cen,r'$\sigma_{(%d,%d)}$'%(pi[0]+1,pi[1]+1),ha='center',va='center',fontsize=11)
for v in itertools.product((0,1),repeat=2):
    ax.plot(*v,'ko',ms=5); ax.annotate(f'{lab(v,(0,0))}',v,xytext=(v[0]*1.12-0.06,v[1]*1.12-0.06),ha='center',va='center',fontsize=10,color='0.2')
ax.plot([0,1],[0,1],color='k',lw=2.2)
ax.set_xlim(-0.2,1.2);ax.set_ylim(-0.2,1.2);ax.set_aspect('equal');ax.axis('off');ax.set_title('(a)',fontsize=11)
E=[(a,b) for a in range(8) for b in range(a+1,8) if bin(a^b).count('1')==1]
P=lambda b: np.array([(b>>i)&1 for i in range(3)],float)
gs=fig.add_gridspec(2,5,left=0,right=1,wspace=0.0,hspace=0.12)
ax.set_position(gs[:,0:2].get_position(fig))
for k,pi in enumerate(itertools.permutations(range(3))):
    a3=fig.add_subplot(gs[k//3,2+k%3],projection='3d')
    V=np.array(K1((0,0,0),pi),float)
    faces=[V[list(f)] for f in itertools.combinations(range(4),3)]
    a3.add_collection3d(Poly3DCollection(faces,facecolors=C[k],alpha=0.6,edgecolors='k',linewidths=0.7))
    for x,y in E: a3.plot(*zip(P(x),P(y)),color='0.45',lw=0.6)
    a3.plot(*zip(P(0),P(7)),color='k',lw=1.6)
    labs=[sum(int(V[i][j])<<j for j in range(3)) for i in range(4)]
    a3.set_box_aspect((1,1,1));a3.view_init(elev=20,azim=-60);a3.axis('off')
    a3.set_xlim(0,1);a3.set_ylim(0,1);a3.set_zlim(0,1)
    a3.set_title('$\\pi=(%d,%d,%d)$: %s'%(pi[0]+1,pi[1]+1,pi[2]+1,','.join(map(str,labs))),fontsize=8,pad=-4)
fig.text(0.7,0.95,'(b)',ha='center',fontsize=11)
plt.savefig(OUT+'tri_cfk.pdf'); plt.close()
# ---------- Fig 4.2: K1 and J1 of R^2
fig,axs=plt.subplots(1,2,figsize=(8.4,3.2))
grid(axs[0],6,4,'K'); grid(axs[1],6,4,'J')
for x in range(0,7):
    for y in range(0,5):
        axs[0].plot(x,y,'o',color='k',ms=2.5)
        if x%2==0 and y%2==0: axs[1].plot(x,y,'o',color=C[3],ms=5)
for g,k in [((1,1),0),((4,2),2)]:
    for kk,pi in enumerate(itertools.permutations(range(2))):
        axs[0].add_patch(plt.Polygon(K1(g,pi),fc=C[k+kk],alpha=0.45,ec='k',lw=1))
        axs[1].add_patch(plt.Polygon(J1(g,pi),fc=C[k+kk],alpha=0.45,ec='k',lw=1))
axs[0].set_title('(a) $K_1$',fontsize=11)
axs[1].set_title('(b) $J_1$',fontsize=11)
plt.tight_layout(); plt.savefig(OUT+'tri_k1j1.pdf'); plt.close()
# ---------- Fig 4.3 / 4.4: a simplex and its three pivots
def pivfig(kind,g,pi,piv,fn,title):
    fig,ax=plt.subplots(figsize=(5.2,3.6)); grid(ax,6,4,kind)
    S=K1 if kind=='K' else J1
    V=S(g,pi)
    ax.add_patch(plt.Polygon(V,fc='0.35',alpha=0.55,ec='k',lw=1.6))
    for i,(gp,pp) in enumerate(piv):
        W=S(gp,pp); ax.add_patch(plt.Polygon(W,fc=C[i],alpha=0.45,ec=C[i],lw=1.4))
        c=np.mean(W,0); ax.text(*c,f'$\\sigma_{i}$',ha='center',va='center',fontsize=11)
    c=np.mean(V,0); ax.text(*c,'$\\sigma$',ha='center',va='center',fontsize=12,color='w')
    for i,v in enumerate(V):
        d=np.array(v,float)-c; d=d/np.linalg.norm(d)*0.22
        ax.text(v[0]+d[0],v[1]+d[1],f'$v_{i}$',ha='center',va='center',fontsize=11)
        ax.plot(*v,'ko',ms=4)
    plt.tight_layout(); plt.savefig(OUT+fn); plt.close()
g=(2,1);pi=(0,1)
pivfig('K',g,pi,[((3,1),(1,0)),((2,1),(1,0)),((2,0),(1,0))],'tri_pivk1.pdf',
       '$\\sigma=K_1((2,1),(1,2))$ and its pivots $\\sigma_i$ across the facets opposite $v_i$')
g=(1,1);pi=(0,1)   # s=(-1,-1), v0=(2,2)
pivfig('J',g,pi,[((0,1),(0,1)),((1,1),(1,0)),((1,2),(0,1))],'tri_pivj1.pdf',
       '$\\sigma$ in the cell $(1,1)$ of $J_1$, $\\pi=(1,2)$, and its pivots $\\sigma_i$')
print('done')
