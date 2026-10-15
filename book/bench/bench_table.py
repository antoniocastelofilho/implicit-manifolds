import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys; sys.path.insert(0,REPO + '/book/runs/bench')
from parse import res, cases
def t(x): 
    if x is None: return '--'
    return f'{x:.3f}' if x<1 else (f'{x:.2f}' if x<100 else f'{x:.0f}')
def c(x): return '--' if x is None else f'{x:,}'.replace(',','\\,')
rows=[]
for key,name,nk,grid in cases:
    if key not in res: continue
    r=res[key]
    rows.append(f"{name} & ${nk}$ & {grid} & {t(r['ms'][0])} & {t(r['ms_j1'][0])} & {t(r['cs'][0])} & {t(r['pta'][0])} & {t(r['gcch'][0])} & {t(r['fch'][0])}\\\\")
rows2=[]
for key,name,nk,grid in cases:
    if key not in res: continue
    r=res[key]
    rows2.append(f"{name} & {c(r['ms'][1])} & {c(r['ms_j1'][1])} & {c(r['pta'][1])} & {c(r['gcch'][1])} & {r['cs'][2]:.1f} & {r['pta'][2]:.1f} & {r['gcch'][2]:.1f} & {r['fch'][2]:.1f}\\\\")
tex=r'''\begin{table}[t]
\caption{Running times in seconds of the C++ implementations on the same grids: Marching Simplex on $K_1$ and
on $J_1$ (MS), Continuation Simplex on $K_1$ (CS), permutahedral tracing (PTA), Generalized Combinatorial
Continuation Hypercubes (GCCH) and Fast Continuation Hypercubes (FCH). Times are those of the approximation
step; the construction of the skeleton by FCH is reported in the text. One core of an Intel Xeon at 2.8\,GHz.}
\label{tab:cmp:cpp}
\centering
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}lcl rrrrrr@{}}
\toprule
Manifold & $(n,k)$ & grid & MS $K_1$ & MS $J_1$ & CS & PTA & GCCH & FCH\\
\midrule
'''+'\n'.join(rows)+r'''
\bottomrule
\end{tabular}
\end{table}
\begin{table}[t]
\caption{Size of the output and peak memory (MB) of the runs of Table~\ref{tab:cmp:cpp}: crossed $n$-simplices
of $K_1$ and $J_1$ (the cells of Marching and Continuation Simplex), transversal $k$-simplices traced by
PTA, and crossed cubes (the cells of GCCH and FCH).}
\label{tab:cmp:cppcells}
\centering
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}l rrrr rrrr@{}}
\toprule
 & \multicolumn{4}{c}{cells / simplices} & \multicolumn{4}{c}{memory (MB)}\\
\cmidrule(lr){2-5}\cmidrule(lr){6-9}
Manifold & $K_1$ & $J_1$ & PTA & cubes & CS & PTA & GCCH & FCH\\
\midrule
'''+'\n'.join(rows2)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
open(REPO + '/book/tables/cpp_runs.tex','w').write(tex); print(tex)
