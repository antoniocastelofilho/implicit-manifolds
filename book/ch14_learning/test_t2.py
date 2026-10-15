import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time; sys.path.insert(0,'.'); sys.path.insert(0,REPO + '/python/fch_lib')
import numpy as np
from pui import PUI
rng=np.random.default_rng(1)
Np=20000
u=rng.uniform(0,2*np.pi,Np); v=rng.uniform(0,2*np.pi,Np)
X=np.c_[np.cos(u),np.sin(u),0.6*np.cos(v),0.6*np.sin(v)]+0.01*rng.normal(size=(Np,4))
dmin=np.array([-1.27,-1.26,-0.83,-0.82]); dmax=np.array([1.25,1.26,0.81,0.82]); div=np.array([12,12,8,8])
t=time.time(); M=PUI(X,dmin,dmax,div,m=2); print('cells',len(M.cells),time.time()-t)
shim,ev=M.cfunc()
y=np.zeros(2)
U=rng.uniform(0,2*np.pi,(200,2)); E=[]
for a,b in U:
    x=np.array([np.cos(a),np.sin(a),0.6*np.cos(b),0.6*np.sin(b)]); ev(x,y); E.append(np.abs(y).max())
print('max |F| on M',max(E))
x=np.array([0.,0.,0.,0.]); ev(x,y); print('F(0)',y)
x=np.array([1.2,0,0.6,0]); ev(x,y); print('F(off)',y)
print('eig ratios',np.median(M.eig,0))
from runfch import run_fch
x0=X[0]
pol,tm=run_fch(shim,4,2,dmin,dmax,div*2,x0,'t2_test'); print(pol,tm)
sys.path.insert(0,REPO + '/python/implicit_manifolds'); from polfile import read_pol, euler
n,k,D,cells=read_pol(pol); print(len(cells), euler(cells))
