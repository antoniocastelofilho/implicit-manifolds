import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'mathtext.fontset':'cm'})
OUT=REPO + '/book/figures/'
x=np.linspace(-1.6,2.6,700); y=np.linspace(-1.6,1.6,540); X,Y=np.meshgrid(x,y)
f1=X**2+Y**2; f2=(X-1)**2+Y**2
cin='#5dade2'; cout='#f5f5f5'
def panel(ax,f,title,lab=None):
    ax.contourf(X,Y,f,[-1e9,1,1e9],colors=[cin,cout])
    ax.contour(X,Y,f,[1],colors='k',linewidths=1.6)
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title,fontsize=10)
    for s in ax.spines.values(): s.set_color('0.7')
    if lab:
        for (px,py,t) in lab: ax.text(px,py,t,ha='center',va='center',fontsize=11)
fig,axs=plt.subplots(1,4,figsize=(12,2.7))
panel(axs[0],f1,'$\\mathcal{M}_1$',[(0,0,'$\\mathcal{M}_1^-$'),(1.9,1.1,'$\\mathcal{M}_1^+$')])
panel(axs[1],f2,'$\\mathcal{M}_2$',[(1,0,'$\\mathcal{M}_2^-$'),(-0.9,1.1,'$\\mathcal{M}_2^+$')])
panel(axs[2],np.maximum(f1,f2),'$\\mathcal{M}_\\cap=f_{\\max}^{-1}(1)$')
panel(axs[3],np.minimum(f1,f2),'$\\mathcal{M}_\\cup=f_{\\min}^{-1}(1)$')
for a in axs[2:]:
    a.contour(X,Y,f1,[1],colors='0.45',linewidths=0.7,linestyles='--'); a.contour(X,Y,f2,[1],colors='0.45',linewidths=0.7,linestyles='--')
plt.subplots_adjust(0.01,0.02,0.99,0.88,0.06,0); plt.savefig(OUT+'mod_2d.pdf'); plt.close()
# p-norm families
fig,axs=plt.subplots(1,2,figsize=(10,3.6))
ps=[1,2,4,8,32]; cm=plt.get_cmap('viridis')
for ax,(nm,fun,inf) in zip(axs,[('union',lambda p:(f1**-p+f2**-p)**(-1/p),np.minimum(f1,f2)),('intersection',lambda p:(f1**p+f2**p)**(1/p),np.maximum(f1,f2))]):
    for i,p in enumerate(ps):
        with np.errstate(all='ignore'): ax.contour(X,Y,fun(p),[1],colors=[cm(i/(len(ps)))],linewidths=1.3)
        ax.plot([],[],color=cm(i/len(ps)),label=f'$p={p}$')
    ax.contour(X,Y,inf,[1],colors='k',linewidths=1.6,linestyles='--'); ax.plot([],[],'k--',label='$p=\\infty$')
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f'smooth {nm} $\\mathcal{{M}}_{{{"\\cup" if nm=="union" else "\\cap"},p}}$')
    ax.legend(fontsize=8,loc='lower right',frameon=False)
plt.subplots_adjust(0.01,0.02,0.99,0.9,0.04,0); plt.savefig(OUT+'mod_pnorm.pdf'); plt.close()
