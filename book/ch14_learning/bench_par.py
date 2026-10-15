import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, os, pickle; sys.path.insert(0,'.')
import numpy as np
from pui_fast import PUIFast
rng=np.random.default_rng(2026)
def torus6(N):
    u=rng.uniform(0,2*np.pi,N); v=rng.uniform(0,2*np.pi,N)
    return np.c_[np.cos(u),np.sin(u),np.cos(v),np.sin(v),0.5*np.cos(u+v),0.5*np.sin(u+v)]+0.01*rng.normal(size=(N,6))
def torus4(N):
    u=rng.uniform(0,2*np.pi,N); v=rng.uniform(0,2*np.pi,N)
    return np.c_[np.cos(u),np.sin(u),0.6*np.cos(v),0.6*np.sin(v)]+0.01*rng.normal(size=(N,4))
cases={'torus4':(torus4(20000),[-1.27,-1.26,-0.83,-0.82],[1.25,1.26,0.81,0.82],[12,12,8,8]),
       'torus6':(torus6(40000),[-1.17,-1.16,-1.15,-1.14,-0.63,-0.62],[1.13,1.14,1.15,1.16,0.61,0.62],[12,12,12,12,6,6])}
out={}
for name,(X,a,b,d) in cases.items():
    M=PUIFast(X,a,b,d,2,fill=False)   # warm up / compile
    M=PUIFast(X,a,b,d,2,fill=False); out[name]=dict(M.timings, cells=M.ndata, threads=os.environ.get('NUMBA_NUM_THREADS','2'))
    print(name, out[name], flush=True)
pickle.dump(out,open(REPO + '/book/runs/ml/bench_par_%s.pkl'%os.environ.get('NUMBA_NUM_THREADS','2'),'wb'))
