import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, itertools; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from polfile import read_pol, polygons
from render import slice3
plt.rcParams.update({'font.size':9,'mathtext.fontset':'stix','font.family':'STIXGeneral'})
OUT=REPO + '/book/figures/'; R=REPO + '/book/runs/ch5/'
C=['#4C72B0','#DD8452','#55A868','#C44E52','#8172B3','#937860']
Lgt=np.array([0.35,-0.45,0.82]);Lgt/=np.linalg.norm(Lgt)
def shade(Q,base):
    n=np.cross(Q[1]-Q[0],Q[2]-Q[0]); s=0.5+0.5*abs(n@Lgt)/(np.linalg.norm(n)+1e-15)
    return tuple(np.clip(np.array(base)*s,0,1))
def cmapc(v,vmin,vmax,cm='viridis'):
    return np.array(plt.get_cmap(cm)((v-vmin)/(vmax-vmin+1e-12)))[:3]
def ax3(fig,pos,elev=22,azim=-58):
    a=fig.add_subplot(pos,projection='3d'); a.view_init(elev=elev,azim=azim); a.axis('off'); return a
def equal(a,P,z=1.0):
    P=np.vstack(P); c=(P.max(0)+P.min(0))/2; r=(P.max(0)-P.min(0)).max()/2/z
    a.set_xlim(c[0]-r,c[0]+r);a.set_ylim(c[1]-r,c[1]+r);a.set_zlim(c[2]-r,c[2]+r);a.set_box_aspect((1,1,1))
def drawP(a,P,cols,lw=0.15,ec=(0,0,0,0.35)):
    a.add_collection3d(Poly3DCollection(P,facecolors=cols,edgecolors=ec,linewidths=lw))

# ---------- 1. consistent perturbation
fig,ax=plt.subplots(figsize=(6.4,3.6))
rng=np.random.default_rng(3); xi={c:rng.uniform(-0.16,0.16,2) for c in itertools.product((0,1),repeat=2)}
nx,ny=7,4
P=lambda i,j: np.array([i,j])+xi[(i%2,j%2)]
for i in range(nx+1):
    ax.plot([P(i,j)[0] for j in range(ny+1)],[P(i,j)[1] for j in range(ny+1)],color='0.3',lw=0.9)
    ax.plot([i]*2,[0,ny],color='0.85',lw=0.6,ls='--',zorder=0)
for j in range(ny+1):
    ax.plot([P(i,j)[0] for i in range(nx+1)],[P(i,j)[1] for i in range(nx+1)],color='0.3',lw=0.9)
    ax.plot([0,nx],[j]*2,color='0.85',lw=0.6,ls='--',zorder=0)
for g in itertools.product(range(nx),range(ny)):
    ax.text(g[0]+0.5,g[1]+0.5,f'({g[0]},{g[1]})',ha='center',va='center',fontsize=7,color='0.55')
for i in range(nx+1):
    for j in range(ny+1):
        c=(i%2,j%2); q=P(i,j); k=c[0]+2*c[1]
        ax.annotate('',xy=q,xytext=(i,j),arrowprops=dict(arrowstyle='->',color=C[k],lw=1.2))
        ax.plot(*q,'o',color=C[k],ms=5)
for k,c in enumerate([(0,0),(1,0),(0,1),(1,1)]):
    ax.plot([],[],'o',color=C[k],label=f'$\\xi_{{{k}}}$')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.02),ncol=4,frameon=False)
ax.set_aspect('equal');ax.axis('off');plt.tight_layout();plt.savefig(OUT+'ms_pert.pdf');plt.close()

# ---------- 2. cells in tetrahedra (skeleton)
def tetra_panel(a,signs,k,title):
    V=np.array([[0,0,0],[0,0,1],[0,1,1],[1,1,1]],float)   # CFK simplex, pi=(3,2,1): labels 0,4,6,7
    labs=[0,4,6,7]
    for i,j in itertools.combinations(range(4),2): a.plot(*zip(V[i],V[j]),color='0.2',lw=1.0)
    if k==1:
        f=np.array(signs,float); pts=[]
        for i,j in itertools.combinations(range(4),2):
            if f[i]*f[j]<0: t=f[i]/(f[i]-f[j]); pts.append(V[i]+t*(V[j]-V[i]))
        pts=np.array(pts); c=pts.mean(0)
        # order around centroid
        nrm=np.cross(pts[1]-pts[0],pts[2]-pts[0]); u=pts[0]-c; u/=np.linalg.norm(u); w=np.cross(nrm/np.linalg.norm(nrm),u)
        ang=[np.arctan2((p-c)@w,(p-c)@u) for p in pts]; pts=pts[np.argsort(ang)]
        drawP(a,[pts],[ (*matplotlib.colors.to_rgb(C[3]),0.75)],lw=1.2,ec='k')
        a.scatter(*pts.T,color='k',s=14,depthshade=False)
        for i in range(4): a.scatter(*V[i],s=40,color='k' if f[i]>0 else 'w',edgecolors='k',depthshade=False)
    else:
        # k=2: F=(f1,f2) linear; segment = zero set; draw transversal triangles' zeros
        F=np.array(signs,float)
        pts=[]
        for tri in itertools.combinations(range(4),3):
            A=np.vstack([np.ones(3),F[list(tri)].T]); 
            try: l=np.linalg.solve(A,[1,0,0])
            except: continue
            if (l>=0).all(): pts.append(l@V[list(tri)])
        pts=np.array(pts); a.plot(*pts.T,color=C[3],lw=3); a.scatter(*pts.T,color='k',s=16,depthshade=False)
        for i in range(4): a.scatter(*V[i],s=20,color='0.4',depthshade=False)
    for i in range(4): a.text(*(V[i]+np.array([0.06,-0.08,0.05])),str(labs[i]),fontsize=10)
    a.set_xlim(0,1);a.set_ylim(0,1);a.set_zlim(0,1);a.set_box_aspect((1,1,1));a.set_title(title,fontsize=10)
fig=plt.figure(figsize=(9,3.3))
tetra_panel(ax3(fig,131,18,-62),[-1,-1,-1,1],1,'(a) $k=1$, one positive vertex')
tetra_panel(ax3(fig,132,18,-62),[-1,-1,1,1],1,'(b) $k=1$, two positive vertices')
tetra_panel(ax3(fig,133,18,-62),[[-0.3,-0.4],[0.5,-0.2],[0.2,0.7],[-0.6,0.5]],2,'(c) $k=2$')
plt.subplots_adjust(0,0,1,0.9,0,0); plt.savefig(OUT+'ms_cells.pdf'); plt.close()

# ---------- 3. S^3: cells of a few grid cells, slices with K1 and J1
def cells_in(f,cubes):
    n,k,D,cells=read_pol(R+f); return [c for c in cells if tuple(c.grid) in cubes],cells
fig=plt.figure(figsize=(11,3.9))
few={(7,4,4,5),(7,4,5,5),(6,4,4,5)}
sel,cellsK=cells_in('BOOK_SPHERE4_ms.pol',few)
a=ax3(fig,131,20,-50); Ps=[];cols=[]
for i,c in enumerate(sel):
    for cyc in polygons(c):
        Q=c.coords[cyc][:,:3]; Ps.append(Q); cols.append((*shade(Q,matplotlib.colors.to_rgb(C[i%6])),0.55))
drawP(a,Ps,cols,lw=0.5,ec=(0,0,0,0.7)); equal(a,Ps,1.0)
a.set_title(f'(a) $K_1$: {len(sel)} cells in 3 grid cells',fontsize=10)
for j,(f,t) in enumerate([('BOOK_SPHERE4_ms.pol','(b) $K_1$'),('BOOK_SPHERE4_ms_j1.pol','(c) $J_1$')]):
    n,k,D,cells=read_pol(R+f)
    S=slice3([c for c in cells if c.coords[:,3].min()<=0.3<=c.coords[:,3].max()],3,0.3)
    Q=[s[:,:3] for s in S]; a=ax3(fig,132+j,20,-50)
    drawP(a,Q,[shade(q,cmapc(q[:,2].mean(),-1,1,'Blues')*0.8+0.15) for q in Q],lw=0.2)
    equal(a,Q,1.2); a.set_title(f'{t}: slice $x_4=0.3$, {len(Q)} polygons',fontsize=10)
plt.subplots_adjust(0,0,1,0.92,0,0); plt.savefig(OUT+'ms_s3.png',dpi=220); plt.close()

# ---------- 4. torus in R^4: K1 vs J1 (stereographic)
R0=np.sqrt(5)/2
st=lambda P: np.array([p[:3]*R0/(R0-p[3]) for p in P])
fig=plt.figure(figsize=(10,4))
for j,(f,t) in enumerate([('BOOK_TORUS4_ms.pol','(a) $K_1$'),('BOOK_TORUS4_ms_j1.pol','(b) $J_1$')]):
    n,k,D,cells=read_pol(R+f); Q=[];cv=[]
    for c in cells:
        for cyc in polygons(c):
            P=c.coords[cyc]; Q.append(st(P)); cv.append(P[:,2].mean())
    a=ax3(fig,121+j,28,-60)
    drawP(a,Q,[shade(q,cmapc(v,-0.5,0.5)) for q,v in zip(Q,cv)],lw=0.08)
    equal(a,Q,1.45); a.set_title(f'{t}: {len(Q)} polygons',fontsize=10)
plt.subplots_adjust(0,0,1,0.93,0,0); plt.savefig(OUT+'ms_torus.png',dpi=220); plt.close()

# ---------- 5. sphere S^2: K1 vs J1
fig=plt.figure(figsize=(8,3.8))
for j,(f,t) in enumerate([('BOOK_SPHERE3_ms.pol','(a) $K_1$'),('BOOK_SPHERE3_ms_j1.pol','(b) $J_1$')]):
    n,k,D,cells=read_pol(R+f); Q=[]
    for c in cells:
        for cyc in polygons(c): Q.append(c.coords[cyc])
    a=ax3(fig,121+j,25,-45)
    drawP(a,Q,[shade(q,(0.55,0.68,0.88)) for q in Q],lw=0.2,ec=(0,0,0,0.6))
    equal(a,Q,1.35); a.set_title(f'{t}: {len(Q)} polygons',fontsize=10)
plt.subplots_adjust(0,0,1,0.93,0,0); plt.savefig(OUT+'ms_sphere_k1j1.png',dpi=220); plt.close()

# ---------- 6. S^3 on J1a: slices near and far from N
n,k,D,cells=read_pol(REPO + '/book/runs/j1a/s3_L2.pol')
lc={0:(0.80,0.82,0.86),1:(0.55,0.68,0.85),2:(0.85,0.40,0.35)}
fig=plt.figure(figsize=(10,4))
for j,(ax_,c4,keep) in enumerate([(0,0.05,[1,2,3]),(3,-0.4,[0,1,2])]):
    sub=[c for c in cells if c.coords[:,ax_].min()<=c4<=c.coords[:,ax_].max()]
    Q=[];cols=[]
    for c in sub:
        for q in slice3([c],ax_,c4): Q.append(q[:,keep]); cols.append(shade(q[:,keep],lc[c.labels[0][0]]))
    a=ax3(fig,121+j,22,-50); drawP(a,Q,cols,lw=0.1,ec=(0,0,0,0.3)); equal(a,Q,1.25)
    a.set_title(f'({"ab"[j]}) slice $x_{ax_+1}={c4}$, {len(Q)} polygons',fontsize=10)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=lc[l],label=f'{l}-blocks') for l in range(3)],loc='lower center',ncol=3,frameon=False)
plt.subplots_adjust(0,0.07,1,0.93,0,0); plt.savefig(OUT+'ms_j1a_s3.png',dpi=220); plt.close()
print('ok')
