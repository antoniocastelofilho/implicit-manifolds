import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'mathtext.fontset':'cm'})
f={(0,0):1.0,(1,0):-1.0,(0,1):-0.8,(1,1):0.6}
def lerp(a,b):
    fa,fb=f[a],f[b]; t=fa/(fa-fb); return np.array(a)+t*(np.array(b)-np.array(a))
B=lerp((0,0),(1,0)); R=lerp((1,0),(1,1)); T=lerp((0,1),(1,1)); L=lerp((0,0),(0,1))
x=np.linspace(0,1,300); X,Y=np.meshgrid(x,x)
Z=f[0,0]*(1-X)*(1-Y)+f[1,0]*X*(1-Y)+f[0,1]*(1-X)*Y+f[1,1]*X*Y
fig,axs=plt.subplots(1,3,figsize=(10,3.6))
def base(ax,title):
    ax.plot([0,1,1,0,0],[0,0,1,1,0],'k',lw=1.2)
    for (i,j),v in f.items():
        ax.plot(i,j,'o',ms=11,color='#c0392b' if v>0 else '#2e86c1',zorder=5)
        ax.text(i+(-.09 if i==0 else .09),j+(-.09 if j==0 else .09),('+' if v>0 else '$-$')+f'{abs(v):g}',ha='center',va='center',fontsize=9)
    for p in (B,R,T,L): ax.plot(*p,'s',color='k',ms=5,zorder=6)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_xlim(-.2,1.2); ax.set_ylim(-.2,1.2); ax.set_title(title,fontsize=10)
# (a) separate positive corners
ax=axs[0]; base(ax,'(a) positive corners separated')
ax.contour(X,Y,Z,[0],colors='0.7',linewidths=1,linestyles='--')
ax.plot(*np.c_[B,L],color='k',lw=2.2); ax.plot(*np.c_[R,T],color='k',lw=2.2)
ax.fill([0,B[0],L[0]],[0,B[1],L[1]],color='#c0392b',alpha=.15); ax.fill([1,R[0],T[0]],[1,R[1],T[1]],color='#c0392b',alpha=.15)
# (b) positive corners joined (asymptotic decider here: f00 f11 - f10 f01 = 0.6-0.8 <0 -> separated) show other pairing
ax=axs[1]; base(ax,'(b) positive corners joined')
ax.contour(X,Y,Z,[0],colors='0.7',linewidths=1,linestyles='--')
ax.plot(*np.c_[B,R],color='k',lw=2.2); ax.plot(*np.c_[L,T],color='k',lw=2.2)
ax.fill([1,B[0],R[0]],[0,B[1],R[1]],color='#2e86c1',alpha=.15); ax.fill([0,L[0],T[0]],[1,L[1],T[1]],color='#2e86c1',alpha=.15)
# (c) Freudenthal-Kuhn diagonal and the PL interpolant
ax=axs[2]; base(ax,'(c) Freudenthal–Kuhn diagonal')
ax.plot([0,1],[0,1],color='0.4',lw=1)
# PL interpolant zero set: lower triangle (00,10,11), upper (00,01,11); diagonal both positive
ax.plot(*np.c_[B,R],color='k',lw=2.2); ax.plot(*np.c_[L,T],color='k',lw=2.2)
ax.fill([1,B[0],R[0]],[0,B[1],R[1]],color='#2e86c1',alpha=.15); ax.fill([0,L[0],T[0]],[1,L[1],T[1]],color='#2e86c1',alpha=.15)
plt.subplots_adjust(0.01,0.01,0.99,0.9,0.05,0)
plt.savefig(REPO + '/book/figures/cmh_ambig.pdf')
print(f[0,0]*f[1,1]-f[1,0]*f[0,1])
