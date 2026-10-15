import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds'); sys.path.insert(0,REPO + '/python/implicit_manifolds/homology')
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from render import new_ax, _equal
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from dual_homology import dual_top_cubes

# ---------- 1. 2D illustration: grid, Gamma, Khalimsky dual X ----------
F = lambda x,y: (x/0.85)**2 + (y/0.6)**2 + 0.35*np.sin(3*x) - 1
N=7; lo,hi=-1.2,1.2; h=(hi-lo)/N
xs=lo+h*np.arange(N+1)+0.011; ys=lo+h*np.arange(N+1)+0.007
fig,axs=plt.subplots(1,2,figsize=(9.6,4.6))
for ax in axs:
    for x in xs: ax.plot([x,x],[ys[0],ys[-1]],color='0.8',lw=0.6)
    for y in ys: ax.plot([xs[0],xs[-1]],[y,y],color='0.8',lw=0.6)
    ax.set_aspect('equal'); ax.axis('off')
V=F(*np.meshgrid(xs,ys,indexing='ij'))
def z(p,q,a,b): t=a/(a-b); return p+t*(q-p)
pts={}
for i in range(N+1):
    for j in range(N+1):
        if i<N and V[i,j]*V[i+1,j]<0: pts[('x',i,j)]=z(np.array([xs[i],ys[j]]),np.array([xs[i+1],ys[j]]),V[i,j],V[i+1,j])
        if j<N and V[i,j]*V[i,j+1]<0: pts[('y',i,j)]=z(np.array([xs[i],ys[j]]),np.array([xs[i],ys[j+1]]),V[i,j],V[i,j+1])
# arcs per square (no ambiguity here)
for i in range(N):
    for j in range(N):
        e=[k for k in (('x',i,j),('x',i,j+1),('y',i,j),('y',i+1,j)) if k in pts]
        if len(e)==2:
            P=np.array([pts[e[0]],pts[e[1]]]); axs[0].plot(*P.T,color='C0',lw=1.8)
        for kk in e: pass
for p in pts.values(): axs[0].plot(*p,'o',color='C3',ms=4)
# shade a primal top cell (square) and a 1-cell
axs[0].set_title(r'primal: $0$-cells on edges, $1$-cells in squares',fontsize=9)
# dual: for each crossed edge J, D(J) joins the centres of the two squares sharing J
cen=lambda i,j: np.array([xs[i]+h/2, ys[j]+h/2])
for (t,i,j),p in pts.items():
    if t=='x':   a,b=cen(i,j-1),cen(i,j)   # horizontal edge, fixed y -> vertical dual segment
    else:        a,b=cen(i-1,j),cen(i,j)
    axs[1].plot([a[0],b[0]],[a[1],b[1]],color='C2',lw=2.2)
    axs[1].plot(*p,'o',color='C3',ms=3,alpha=0.6)
for i in range(N):
    for j in range(N):
        if any(k in pts for k in (('x',i,j),('x',i,j+1),('y',i,j),('y',i+1,j))):
            axs[1].plot(*cen(i,j),'s',color='k',ms=3.5)
axs[1].set_title(r'dual: cubical set $X$ through the centres of the crossed squares',fontsize=9)
plt.tight_layout(); plt.savefig(REPO + '/book/figures/hom_dual2d.pdf')

# ---------- 2. 3D: dual cubical set X of a torus ----------
torus = lambda x,y,z: (np.sqrt(x**2+y**2)-0.6)**2 + z**2 - 0.25**2
N=16; lo,hi=-1,1; h=(hi-lo)/N
tops,_=dual_top_cubes(torus,3,1,lo,hi,N)
sq=[]
for c,mask in tops:
    fixed=[i for i in range(3) if mask>>i&1]
    base=np.array([lo+(ci+0.5)*h for ci in c])
    a,b=fixed
    ea=np.zeros(3); ea[a]=h; eb=np.zeros(3); eb[b]=h
    sq.append(np.array([base,base+ea,base+ea+eb,base+eb]))
fig=plt.figure(figsize=(6,4.2)); ax=new_ax(fig,elev=35,azim=-60)
L=np.array([0.4,0.3,0.85]); L/=np.linalg.norm(L)
cols=[]
for Q in sq:
    nrm=np.cross(Q[1]-Q[0],Q[3]-Q[0]); s=0.45+0.55*abs(nrm@L)/np.linalg.norm(nrm)
    cols.append(tuple(np.array([0.85,0.55,0.25])*s))
ax.add_collection3d(Poly3DCollection(sq,facecolors=cols,edgecolors=(0.2,0.2,0.2,0.5),linewidths=0.25))
_equal(ax,np.vstack(sq),zoom=1.5)
plt.subplots_adjust(0,0,1,1); plt.savefig(REPO + '/book/figures/hom_dual_torus.png',dpi=220)
print('squares',len(sq))

# ---------- 3. Hausdorff rate ----------
Ns=np.array([8,16,32,64,128]); d=np.array([3.347e-2,7.938e-3,1.961e-3,4.888e-4,1.221e-4]); hh=2/Ns
fig,ax=plt.subplots(figsize=(4.6,3.2))
ax.loglog(hh,d,'o-',color='C0',label=r'$\max_\tau d(p_\tau,M)$')
ax.loglog(hh,0.5*hh**2,'--',color='0.4',label=r'$0.5\,h^2$')
ax.set_xlabel('$h$'); ax.legend(frameon=False,fontsize=8); ax.grid(True,which='both',alpha=0.25)
plt.tight_layout(); plt.savefig(REPO + '/book/figures/hom_rate.pdf')

# ---------- 4. ambiguity under refinement ----------
fig,ax=plt.subplots(figsize=(5.2,3.3))
ax.plot([8,16,32,64,128],[0.58,0.25,0.17,0.08,0.0],'o-',label=r'$T^2\subset\mathbb{R}^3$ ($m=2$)')
ax.plot([8,16,32,48,64],[2.7,1.7,1.5,2.7,2.6],'s-',label=r'$S^2\times S^1\subset\mathbb{R}^4$ ($m=3$)')
ax.plot([8,12,16,20,24],[1.8,6.0,6.0,8.7,10.0],'^-',label=r'$S^3\times S^1\subset\mathbb{R}^5$ ($m=4$)')
ax.set_xscale('log',base=2); ax.set_xlabel('cells per direction $N$'); ax.set_ylabel('ambiguous 2-faces per grid')
ax.legend(frameon=False,fontsize=8); ax.grid(True,alpha=0.25)
plt.tight_layout(); plt.savefig(REPO + '/book/figures/hom_ambiguity.pdf')
