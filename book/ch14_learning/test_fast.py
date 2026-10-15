import sys, time; sys.path.insert(0,'.')
import numpy as np
from pui import PUI
from pui_fast import PUIFast
rng=np.random.default_rng(2026)
N=20000; u=rng.uniform(0,2*np.pi,N); v=rng.uniform(0,2*np.pi,N)
X=np.c_[np.cos(u),np.sin(u),0.6*np.cos(v),0.6*np.sin(v)]+0.01*rng.normal(size=(N,4))
args=([-1.27,-1.26,-0.83,-0.82],[1.25,1.26,0.81,0.82],[12,12,8,8],2)
t=time.time(); A=PUI(X,*args); ta=time.time()-t
B=PUIFast(X,*args,verbose=True); B=PUIFast(X,*args,verbose=True)
print('slow',ta,'cells',A.ndata,B.ndata)
ea=A.njit_eval(); eb=B.njit_eval(); y1=np.zeros(2); y2=np.zeros(2); d=0
for x in X[:2000]+0.05*rng.normal(size=(2000,4)):
    ea(x,y1); eb(x,y2); d=max(d,abs(np.linalg.norm(y1)-np.linalg.norm(y2)))
print('max diff',d)
