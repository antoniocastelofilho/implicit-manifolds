import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import re, os, sys
sys.path.insert(0,REPO + '/python/implicit_manifolds')
R=REPO + '/book/runs/bench/'
cases=[('BOOK_SPHERE3','$S^2\\subset\\R^3$','(3,1)','$12^3$'),('BOOK_SPHERE4','$S^3\\subset\\R^4$','(4,1)','$8^4$'),
('BOOK_TORUS4','$T^2\\subset\\R^4$','(4,2)','$10^4$'),('BOOK_ZCOSW','$z=\\cos w$','(4,2)','$10^4$'),
('BOOK_S2XS1','$S^2\\times S^1\\subset\\R^5$','(5,2)','$12^3\\cdot8^2$'),('BOOK_T3_R6_B','$T^3\\subset\\R^6$','(6,3)','$8^6$')]
def g(f,pat,cast=float):
    try:
        m=re.search(pat,open(f).read()); return cast(m.group(1)) if m else None
    except FileNotFoundError: return None
def rss(f): v=g(f,r'Maximum resident set size \(kbytes\): (\d+)',int); return v/1024 if v else None
def el(f):
    m=re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\d+):([\d.]+)',open(f).read()) if os.path.exists(f) else None
    return int(m.group(1))*60+float(m.group(2)) if m else None
res={}
for c,*_ in cases:
    d=R+c+'/'
    if not os.path.exists(d+'skeleton.log'): continue
    r={}
    r['ms']=(g(d+'ms.log',r'Time \(MS\): ([\d.e+-]+)'),g(d+'ms.log',r'crossed: (\d+)',int),rss(d+'ms.log'),g(d+'ms.log',r'F evaluations: (\d+)',int))
    r['ms_j1']=(g(d+'ms_j1.log',r'Time \(MS\): ([\d.e+-]+)'),g(d+'ms_j1.log',r'crossed: (\d+)',int),rss(d+'ms_j1.log'),None)
    r['cs']=(g(d+'cs.log',r'Time \(CS\): ([\d.e+-]+)'),g(d+'cs.log',r'crossed n-simplices: (\d+)',int),rss(d+'cs.log'),g(d+'cs.log',r'F evaluations: (\d+)',int))
    r['pta']=(g(d+'pta.log',r'Time to generate bin \(PTA\): ([\d.e+-]+)'),g(d+'pta.log',r'(\d+) simplices generated',int),rss(d+'pta.log'),None)
    r['gcch']=(g(d+'gcch.log',r'Time to generate bin \(GCCH\): ([\d.e+-]+)'),g(d+'skeleton.log',r'n_saved_hypercubes: (\d+)',int),rss(d+'gcch.log'),None)
    post=sum(x or 0 for x in (el(d+'fch_hc.log'),el(d+'reorder.log'),el(d+'skeleton.log')))
    r['fch']=(g(d+'fch.log',r'Time to generate bin \(FCH\): ([\d.e+-]+)'),g(d+'fch.log',r'(\d+) cells generated',int),rss(d+'fch.log'),post)
    res[c]=r
if __name__=='__main__':
    for c in res: print(c,res[c])
