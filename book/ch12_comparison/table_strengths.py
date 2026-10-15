import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import re, os
B=REPO + '/book/runs/bench/'
def g(f,pat,cast=float):
    try:
        m=re.search(pat,open(f).read()); return cast(m.group(1)) if m else None
    except FileNotFoundError: return None
comp={'BOOK_BLOBS3':(11,1),'BOOK_BLOBS4':(11,1),'BOOK_BRATU3':(1,1),'BOOK_BRATU5':(1,1)}
rows=[('BOOK_BLOBS3','blobs, $\\R^3$','$(3,1)$','$40^3$'),('BOOK_BLOBS4','blobs, $\\R^4$','$(4,1)$','$24^4$'),
      ('BOOK_BRATU3','Bratu, $N=3$','$(4,3)$','$24\\cdot20^2\\cdot30$'),('BOOK_BRATU5','Bratu, $N=5$','$(6,5)$','$12^2\\cdot10^4$')]
def t(x):
    if x is None: return '--'
    return '%.2g'%x if x<0.01 else ('%.3f'%x if x<1 else ('%.1f'%x if x<100 else '{:,}'.format(round(x)).replace(',','\\,')))
def num(x): return '--' if x is None else '{:,}'.format(x).replace(',','\\,')
out=[]
for c,name,nk,grid in rows:
    d=B+c+'/'
    ms=(g(d+'ms.log',r'Time \(MS\): ([\d.e+-]+)'),g(d+'ms.log',r'F evaluations: (\d+)',int),g(d+'ms.log',r'crossed: (\d+)',int))
    cs=(g(d+'cs.log',r'Time \(CS\): ([\d.e+-]+)'),g(d+'cs.log',r'F evaluations: (\d+)',int),g(d+'cs.log',r'crossed n-simplices: (\d+)',int))
    pta=g(d+'pta.log',r'Time to generate bin \(PTA\): ([\d.e+-]+)'); gc=g(d+'gcch.log',r'Time to generate bin \(GCCH\): ([\d.e+-]+)'); fc=g(d+'fch.log',r'Time to generate bin \(FCH\): ([\d.e+-]+)')
    cm,cc=comp[c]
    out.append(f"{name} & {nk} & {t(ms[0])} & {num(ms[1])} & {cm} & {t(cs[0])} & {num(cs[1])} & {cc} & {t(pta)} & {t(gc)} & {t(fc)}\\\\")
T=r"""\begin{table}[t]
\caption{Marching and continuation on problems chosen to show their strengths. Marching Simplex (MS) and Continuation
Simplex (CS): time in seconds, evaluations of $F$ and number of components found; the continuation methods start from one
point. Last columns: time of the permutahedral tracing algorithm, of the generalized continuation hypercubes and of
Fast Continuation Hypercubes, which find the same component as CS. Grids: $40^3$ and $24^4$ cubes for the blobs,
$24\cdot20^2\cdot30$ and $12^2\cdot10^4$ for the Bratu problems.}
\label{tab:cmp:strengths}
\centering\footnotesize\setlength{\tabcolsep}{2.2pt}
\begin{tabular}{@{}llrrcrrcrrr@{}}
\toprule
& & \multicolumn{3}{c}{MS} & \multicolumn{3}{c}{CS} & PTA & GCCH & FCH\\
\cmidrule(lr){3-5}\cmidrule(lr){6-8}\cmidrule(lr){9-9}\cmidrule(lr){10-10}\cmidrule(l){11-11}
Problem & $(n,k)$ & time & eval. & comp. & time & eval. & comp. & time & time & time\\
\midrule
"""+'\n'.join(out)+r"""
\bottomrule
\end{tabular}
\end{table}
"""
open(REPO + '/book/tables/cmp_strengths.tex','w').write(T); print(T)
