import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import pickle, numpy as np
R=REPO + '/book/runs/ml/'
rows=[('curve','closed curve','$(2,1)$'),('surface','bumpy torus','$(3,1)$'),('swiss','Swiss roll','$(3,1)$'),('torus4','torus','$(4,2)$'),('torus6','torus','$(6,4)$')]
T=pickle.load(open(R+'tda.pkl','rb'))
out=[]
for tag,name,nk in rows:
    r=pickle.load(open(R+tag+'.pkl','rb'))['r']
    from itertools import groupby
    grid='$'+r'\times'.join((f'{k}^{{{len(list(g))}}}' if len(l:=list(g2:=[k]*0))==0 else '') for k,g in [(k,list(g)) for k,g in groupby(r['div_pui'])]) +'$'
    gg=[(k,len(list(g))) for k,g in groupby(r['div_pui'])]
    grid='$'+r'\times'.join(f'{k}^{{{c}}}' if c>1 else f'{k}' for k,c in gg)+'$'
    eb=np.sqrt(np.mean(r['err_before']**2)); ea=np.sqrt(np.mean(r['err_after']**2))
    b=T[tag]['betti']
    out.append(f"{name} & {nk} & {r['npts']:,} & {grid} & {r['cells']:,} & ${r['chi'][1]}$ & $({', '.join(map(str,b))})$ & {eb:.4f} & {ea:.4f} & {r['t_pui']:.1f} & {r['t_fch']['fch']:.2f}\\\\".replace(',',r'\,') if False else
               f"{name} & {nk} & {r['npts']:,} & {grid} & {r['cells']:,} & ${r['chi'][1]}$ & $({', '.join(map(str,b))})$ & {eb:.4f} & {ea:.4f} & {r['t_pui']:.1f} & {r['t_fch']['fch']:.2f}\\\\")
fix=lambda s: s.replace(',000','\\,000').replace(',',',')
tex=r'''\begin{table}[t]
\caption{Partition of unity implicits computed from samples and approximated by Fast Continuation Hypercubes.
$(n,k)$: ambient dimension and codimension; grid of the partition of unity; number of cells of the
approximation, its Euler characteristic and its Betti numbers over $\mathbb{Z}_2$; root mean square distance of
the sample to the exact manifold before and after projection onto the cells; times in seconds of the
construction of the partition of unity (parallel version, two threads) and of the continuation.}
\label{tab:ml:runs}
\centering
\small
\setlength{\tabcolsep}{3pt}
\begin{tabular}{@{}lcrlrrcccrr@{}}
\toprule
Sample & $(n,k)$ & points & grid & cells & $\chi$ & Betti & before & after & PUI & FCH\\
\midrule
'''+'\n'.join(out)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
import re
tex=re.sub(r'(\d),(\d\d\d)',r'\1\\,\2',tex)
open(REPO + '/book/tables/ml_runs.tex','w').write(tex); print(tex)
