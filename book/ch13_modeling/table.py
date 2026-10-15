import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import re, pickle, os
R=REPO + '/book/runs/modeling/'
info=pickle.load(open(R+'info.pkl','rb'))
rows=[('Two spheres, sharp union','MOD_TWO_SPHERES',3,'$40\\times28^2$'),('Two spheres, $p=8$','MOD_TWO_SPHERES_P8',3,'$40\\times28^2$'),
('Eight tori, $p=10$','MOD_EIGHT_TORI',3,'$80^2\\times20$'),('Knot of spheres, $p=30$','MOD_KNOT',3,'$72^2\\times42$'),
('Three spheres, $p=10$','MOD_INTERSECT3',3,'$64^2\\times44$'),('Pipe/box, sharp','MOD_PIPEBOX',3,'$56^2\\times112$'),('Pipe/box, $p=8$','MOD_PIPEBOX_P8',3,'$56^2\\times112$'),
('Two balls, $p=8$','MOD_R4_TWO_BALLS',4,'$36\\times22^3$'),('Drilled ball, $p=6$','MOD_R4_DIFF',4,'$32^4$'),('Tube around a curve','MOD_R4_TUBE',4,'$44^2\\times22^2$')]
def num(f,pat):
    m=re.search(pat,open(f).read()); return float(m.group(1)) if m else None
def cnt(path):
    n=0
    with open(path) as f:
        pass
    return n
out=[]
for name,case,n,grid in rows:
    d=R+case+'/'
    tms=num(d+'ms.log',r'Time \(MS\): ([0-9.e+-]+)'); cr=int(num(d+'ms.log',r'crossed: (\d+)'))
    tf=num(d+'fch.log',r'Time to generate bin \(FCH\): ([0-9.e+-]+)')
    ts=num(d+'skeleton.log',r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): 0:([0-9.]+)') or 0
    tc=num(d+'fch_to_hc.log',r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): 0:([0-9.]+)') or 0
    fc=info.get((case,'fch'),(None,None))[0]
    if fc is None:
        import sys; sys.path.insert(0,REPO + '/python/implicit_manifolds'); from slicer import iter_pol
        fc=sum(1 for _ in iter_pol(d+'fch.pol'))
    chi=info.get((case,'fch'),(None,None))[1]
    chi_s='--' if n==4 else f'${chi}$'
    th=lambda v: f'{v:,}'.replace(',','\\,')
    out.append(f"{name} & {n} & {grid} & {th(cr)} & {tms:.2f} & {th(fc)} & {tf:.2f} & {tf+tc+ts:.2f} & {chi_s}\\\\")
tex=r'''\begin{table}[t]
\caption{Implicit models approximated by Marching Simplex on $K_1$ (MS) and by Fast Continuation Hypercubes
(FCH) on the same grids: number of cells, time of the approximation (seconds) and, for FCH, total time
including the construction of the skeleton. For the three-dimensional models both methods give the same
Euler characteristic.}
\label{tab:mod:runs}
\centering
\small
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{@{}lclrrrrrr@{}}
\toprule
Model & $n$ & grid & MS cells & MS (s) & FCH cells & FCH (s) & total (s) & $\chi$\\
\midrule
'''+'\n'.join(out)+r'''
\bottomrule
\end{tabular}
\end{table}
'''
open(REPO + '/book/tables/mod_runs.tex','w').write(tex); print(tex)
