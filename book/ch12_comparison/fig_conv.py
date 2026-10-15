import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, os, pickle; sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from slicer import iter_pol
from polfile import polygons
plt.rcParams.update({'font.size':9,'mathtext.fontset':'cm','axes.spines.top':False,'axes.spines.right':False})
R=REPO + '/book/runs/bench/conv/'; FG=REPO + '/book/figures/'
def err(case,X):
    if case=='BOOK_CONV_S2': return np.abs(np.linalg.norm(X,axis=1)-1)
    return np.hypot(np.linalg.norm(X[:,:2],axis=1)-1, np.linalg.norm(X[:,2:4],axis=1)-0.5)
data={}
import os
if os.path.exists(REPO + '/book/ch12_comparison/conv.pkl'): data=pickle.load(open(REPO + '/book/ch12_comparison/conv.pkl','rb'))
for case in ('BOOK_CONV_S2','BOOK_CONV_T2'):
    for N in (8,16,32,64):
        d=R+f'{case}_{N}/'
        if not os.path.exists(d+'fch.pol') or (case,N,'fch') in data: continue
        for alg in ('ms','fch'):
            E=[]; sizes=[]; cells=0
            for c in iter_pol(d+alg+'.pol'):
                cells+=1; E.append(err(case,c.coords))
                for cyc in polygons(c): sizes.append(len(cyc))
            E=np.concatenate(E)
            data[(case,N,alg)]=dict(h=2.25/N,max=E.max(),rms=np.sqrt((E**2).mean()),cells=cells,sizes=np.bincount(sizes))
            print(case,N,alg,cells,E.max(),flush=True)
pickle.dump(data,open(REPO + '/book/ch12_comparison/conv.pkl','wb'))
col={'ms':'#2a78d6','fch':'#eb6834'}; lab={'ms':'Marching Simplex','fch':'FCH'}
fig,axs=plt.subplots(1,2,figsize=(8.4,3.3))
for ax,(case,title) in zip(axs,[('BOOK_CONV_S2','$S^2\\subset\\mathbb{R}^3$'),('BOOK_CONV_T2','$T^2\\subset\\mathbb{R}^4$')]):
    for alg in ('ms','fch'):
        pts=sorted((v['h'],v['max'],v['rms']) for (c,N,a),v in data.items() if c==case and a==alg)
        if not pts: continue
        h,mx,rm=map(np.array,zip(*pts))
        ax.loglog(h,mx,'-o',color=col[alg],lw=2,ms=4,label=lab[alg]+', max')
        ax.loglog(h,rm,'--s',color=col[alg],lw=1.5,ms=3.5,label=lab[alg]+', rms')
    h0=np.array([h.min(),h.max()]); ax.loglog(h0,0.25*h0**2,color='0.55',lw=1); hm=np.sqrt(h0[0]*h0[1]); ax.text(hm,0.25*hm**2*0.5,'slope 2',color='0.4',fontsize=8,ha='left',va='top')
    ax.set_xlabel('grid spacing $h$'); ax.minorticks_off(); ax.set_xticks(h); ax.set_xticklabels([f'{x:.2f}' for x in h]); ax.set_title(title,fontsize=10); ax.grid(color='0.9',lw=0.5,which='both')
axs[0].set_ylabel('distance of vertices to $\\mathcal{M}$'); axs[0].legend(frameon=False,fontsize=7.5)
plt.tight_layout(); plt.savefig(FG+'cmp_conv.pdf')
# vertices per 2-cell at N=32
fig,axs=plt.subplots(1,2,figsize=(8.4,2.9))
for ax,(case,title) in zip(axs,[('BOOK_CONV_S2','$S^2$, $32^3$ cubes'),('BOOK_CONV_T2','$T^2\\subset\\mathbb{R}^4$, $32^2\\cdot16^2$ cubes')]):
    w=0.38
    for j,alg in enumerate(('ms','fch')):
        v=data.get((case,32,alg))
        if v is None: continue
        s=v['sizes']; k=np.arange(len(s)); m=s>0
        ax.bar(k[m]+(j-0.5)*w,s[m]/s.sum(),width=w,color=col[alg],label=f"{lab[alg]} ({v['cells']} cells)")
    ax.set_xlabel('vertices per polygon'); ax.set_xticks(range(3,7)); ax.set_title(title,fontsize=10); ax.grid(axis='y',color='0.9',lw=0.5)
axs[0].set_ylabel('fraction of polygons'); [a.legend(frameon=False,fontsize=7.5) for a in axs]
plt.tight_layout(); plt.savefig(FG+'cmp_cells.pdf')
