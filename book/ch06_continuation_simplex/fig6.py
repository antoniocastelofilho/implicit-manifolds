import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, itertools; sys.path.insert(0,REPO + '/python/implicit_manifolds'); sys.path.insert(0,REPO + '/book/ch05_marching_simplex')
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from polfile import read_pol, polygons
plt.rcParams.update({'font.size':9,'mathtext.fontset':'stix','font.family':'STIXGeneral'})
OUT=REPO + '/book/figures/'; R=REPO + '/book/runs/ch5/'
C=['#4C72B0','#DD8452','#55A868','#C44E52','#8172B3','#937860']
Lgt=np.array([0.35,-0.45,0.82]);Lgt/=np.linalg.norm(Lgt)
def shade(Q,base):
    n=np.cross(Q[1]-Q[0],Q[2]-Q[0]); s=0.5+0.5*abs(n@Lgt)/(np.linalg.norm(n)+1e-15)
    return tuple(np.clip(np.array(base)*s,0,1))
def equal(a,P,z=1.0):
    P=np.vstack(P); c=(P.max(0)+P.min(0))/2; r=(P.max(0)-P.min(0)).max()/2/z
    a.set_xlim(c[0]-r,c[0]+r);a.set_ylim(c[1]-r,c[1]+r);a.set_zlim(c[2]-r,c[2]+r);a.set_box_aspect((1,1,1))
# ---------- 1. barycentric coordinates
fig,axs=plt.subplots(1,3,figsize=(9,2.9))
V=np.array([[0,0],[4,0.4],[1.4,3.0]]); x=np.array([3.4,2.3])
def foot(p,a,b):
    d=b-a; t=(p-a)@d/(d@d); return a+t*d
for i,ax in enumerate(axs):
    a_,b_=V[(i+1)%3],V[(i+2)%3]; col=C[i]
    ax.add_patch(plt.Polygon(V,fc='#DCE6F2',ec='k',lw=1.0))
    d=b_-a_; ax.plot(*zip(a_-d*0.35,b_+d*0.35),color=col,lw=1.2)
    fv=foot(V[i],a_,b_); fx=foot(x,a_,b_)
    ax.plot(*zip(V[i],fv),color=col,lw=1.6,ls='--'); ax.plot(*zip(x,fx),color='k',lw=1.6)
    ax.text(*((V[i]+fv)/2+[0.12,0]),f'$d_{{v{i+1}}}$',color=col,fontsize=11)
    ax.text(*((x+fx)/2+[0.12,0.05]),f'$d_{{x{i+1}}}$',fontsize=11)
    for j,off in zip(range(3),([-0.35,-0.35],[0.1,-0.35],[-0.1,0.15])): ax.plot(*V[j],'ko',ms=3.5); ax.text(*(V[j]+off),f'$v_{j+1}$',fontsize=11)
    ax.plot(*x,'o',color=C[3],ms=5); ax.text(*(x+[0.1,0.1]),'$x$',fontsize=11,color=C[3])
    sign='' if i>0 else '-'
    ax.set_title(f'$\\lambda_{i+1}={sign}d_{{x{i+1}}}/d_{{v{i+1}}}$',fontsize=11)
    ax.set_xlim(-0.6,4.8);ax.set_ylim(-0.6,3.5);ax.set_aspect('equal');ax.axis('off')
plt.tight_layout();plt.savefig(OUT+'cs_baric.pdf');plt.close()
# ---------- 2. walks in K1 and J1
def verts(kind,g,pi):
    if kind=='K':
        v=np.array(g,float);V=[v.copy()]
        for p in pi: v[p]+=1;V.append(v.copy())
    else:
        s=[(-1)**gi for gi in g]; v=np.array([gi+(gi%2) for gi in g],float);V=[v.copy()]
        for p in pi: v[p]+=s[p];V.append(v.copy())
    return np.array(V)
def pivot(kind,g,pi,i):
    g=list(g);pi=list(pi);n=2
    if kind=='K':
        if i==0: g[pi[0]]+=1; pi=pi[1:]+pi[:1]
        elif i==n: g[pi[-1]]-=1; pi=pi[-1:]+pi[:-1]
        else: pi[i-1],pi[i]=pi[i],pi[i-1]
    else:
        s=lambda d:(-1)**g[d]
        if i==0: g[pi[0]]+=s(pi[0])
        elif i==n: g[pi[-1]]-=s(pi[-1])
        else: pi[i-1],pi[i]=pi[i],pi[i-1]
    return tuple(g),tuple(pi)
def walk(kind,g,pi,x):
    seq=[(g,pi)]
    for _ in range(60):
        V=verts(kind,g,pi); A=np.vstack([np.ones(3),V.T]); l=np.linalg.solve(A,[1,*x])
        if (l>=-1e-12).all(): return seq
        g,pi=pivot(kind,g,pi,int(np.argmin(l))); seq.append((g,pi))
    return seq
fig,axs=plt.subplots(1,2,figsize=(9,3.6))
x=np.array([5.62,3.31])
for ax,kind,t in zip(axs,'KJ',['(a) $K_1$','(b) $J_1$']):
    tris=[verts(kind,g,pi) for g in itertools.product(range(7),range(5)) for pi in itertools.permutations(range(2))]
    ax.add_collection(PolyCollection(tris,facecolors='none',edgecolors='0.7',linewidths=0.5))
    seq=walk(kind,(0,0),(0,1),x)
    cm=plt.get_cmap('Blues')
    for j,(g,pi) in enumerate(seq):
        V=verts(kind,g,pi); ax.add_patch(plt.Polygon(V,fc=cm(0.25+0.6*j/max(1,len(seq)-1)),ec='k',lw=0.8))
        c=V.mean(0); ax.text(*c,str(j),ha='center',va='center',fontsize=7.5,color='k' if j<len(seq)*0.6 else 'w')
    ax.plot(*x,'o',color=C[3],ms=6,zorder=5); ax.text(x[0]+0.12,x[1]+0.12,'$x$',color=C[3],fontsize=12)
    ax.set_xlim(-0.1,7.1);ax.set_ylim(-0.1,5.1);ax.set_aspect('equal');ax.axis('off');ax.set_title(f'{t}: {len(seq)} simplices',fontsize=10)
plt.tight_layout();plt.savefig(OUT+'cs_walk.pdf');plt.close()
print('walk lengths ok')
# ---------- 3. complex curves (K1)
fig=plt.figure(figsize=(11,3.9))
for j,(f,t,el) in enumerate([('BOOK_ZCOSW_cs.pol','(a) $z=\\cos w$',20),('BOOK_Z2COSW_cs.pol','(b) $z^2=\\cos w$',20),('BOOK_Z2W5_cs.pol','(c) $z^2=w^5+1+i$',25)]):
    n,k,D,cells=read_pol(R+f); Q=[];cv=[]
    for c in cells:
        for cyc in polygons(c):
            P=c.coords[cyc]; Q.append(P[:,[0,1,2]]); cv.append(P[:,3].mean())
    cv=np.array(cv); lo,hi=np.percentile(cv,[2,98])
    a=fig.add_subplot(131+j,projection='3d'); a.view_init(elev=el,azim=-55); a.axis('off')
    cm=plt.get_cmap('viridis')
    a.add_collection3d(Poly3DCollection(Q,facecolors=[shade(q,np.array(cm(np.clip((v-lo)/(hi-lo),0,1)))[:3]) for q,v in zip(Q,cv)],edgecolors=(0,0,0,0.25),linewidths=0.08))
    equal(a,Q,1.15); a.set_title(f'{t}: {len(cells)} cells',fontsize=10)
plt.subplots_adjust(0,0,1,0.92,0,0);plt.savefig(OUT+'complex_cs.png',dpi=220);plt.close()
# ---------- 4. z = cos w on J1a
lc={0:(0.80,0.82,0.86),1:(0.55,0.68,0.85),2:(0.85,0.40,0.35)}
fig=plt.figure(figsize=(10,4.2))
for j,L in enumerate([0,2]):
    n,k,D,cells=read_pol(R+f'zcosw_j1a_L{L}.pol'); Q=[];cols=[]
    for c in cells:
        for cyc in polygons(c):
            P=c.coords[cyc][:,[0,1,2]]; Q.append(P); cols.append(shade(P,lc[c.labels[0][0]]))
    a=fig.add_subplot(121+j,projection='3d'); a.view_init(elev=22,azim=-60); a.axis('off')
    a.add_collection3d(Poly3DCollection(Q,facecolors=cols,edgecolors=(0,0,0,0.3),linewidths=0.08)); equal(a,Q,1.2)
    a.set_title(f'({"ab"[j]}) {L} levels, {len(cells)} cells',fontsize=10)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color=lc[l],label=f'{l}-blocks') for l in range(3)],loc='lower center',ncol=3,frameon=False)
plt.subplots_adjust(0,0.07,1,0.93,0,0);plt.savefig(OUT+'cs_j1a_zcosw.png',dpi=220);plt.close()
print('ok')
