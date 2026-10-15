import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import re, os, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8.5, 'mathtext.fontset': 'cm'})
S=REPO + '/book/runs/bench/scale/'
def g(f,pat,cast=float):
    try:
        m=re.search(pat,open(f).read()); return cast(m.group(1)) if m else None
    except FileNotFoundError: return None
cases=[('BOOK_CONV_S2','sphere $S^2\\subset\\mathbb{R}^3$',[16,32,64,128],3,2),('BOOK_CONV_T2','torus $T^2\\subset\\mathbb{R}^4$',[8,16,24,32,48],4,2)]
fig,axs=plt.subplots(1,3,figsize=(10,3.0)); res={}
mk={'MS':'o-','CS':'s-','PTA':'^-','FCH':'D-'}
for ci,(c,title,Ns,n,m) in enumerate(cases):
    T={k:[] for k in mk}; E={'MS':[],'CS':[]}; NN=[]
    for N in Ns:
        d=S+f'{c}_{N}/'
        if not os.path.exists(d+'fch.log'): continue
        t=dict(MS=g(d+'ms.log',r'Time \(MS\): ([\d.e+-]+)'),CS=g(d+'cs.log',r'Time \(CS\): ([\d.e+-]+)'),
               PTA=g(d+'pta.log',r'Time to generate bin \(PTA\): ([\d.e+-]+)'),FCH=g(d+'fch.log',r'Time to generate bin \(FCH\): ([\d.e+-]+)'))
        if None in t.values(): continue
        NN.append(N)
        for k in mk: T[k].append(t[k])
        E['MS'].append(g(d+'ms.log',r'F evaluations: (\d+)',int)); E['CS'].append(g(d+'cs.log',r'F evaluations: (\d+)',int))
    NN=np.array(NN); a=axs[ci]; res[c]={}
    for k in mk:
        y=np.array(T[k]); sl=np.polyfit(np.log(NN[1:]),np.log(y[1:]),1)[0] if len(NN)>2 else np.nan
        res[c][k]=sl; a.loglog(NN,y,mk[k],ms=3.5,label=f'{k} (slope {sl:.1f})')
    a.set_title(title+': time',fontsize=8.5); a.set_xlabel('cubes per direction $N$'); a.legend(fontsize=7); a.set_xticks(NN); a.set_xticklabels(NN); a.minorticks_off()
    if ci==0: a.set_ylabel('seconds')
    for k,st in (('MS','o-'),('CS','s--')):
        y=np.array(E[k],float); sl=np.polyfit(np.log(NN),np.log(y),1)[0]; res[c]['E'+k]=sl
        axs[2].loglog(NN,y,st,ms=3.5,label=f'{k}, $\\mathbb{{R}}^{n}$ (slope {sl:.1f})')
axs[2].set_title('evaluations of $F$',fontsize=8.5); axs[2].set_xlabel('cubes per direction $N$'); axs[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(REPO + '/book/figures/cmp_scale.pdf'); plt.close()
print(res)
