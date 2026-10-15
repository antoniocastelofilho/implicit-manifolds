import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle; sys.path.insert(0,'.'); sys.path.insert(0,REPO + '/python/implicit_manifolds')
import numpy as np
from numba import njit
from pui_fast import PUIFast
from polfile import read_pol, euler
from project import project
from scipy.spatial import cKDTree
rng = np.random.default_rng(7)
N = 40000; u = rng.uniform(0, 2*np.pi, N); v = rng.uniform(0, 2*np.pi, N)
def emb(u, v): return np.c_[np.cos(u), np.sin(u), np.cos(v), np.sin(v), 0.5*np.cos(u+v), 0.5*np.sin(u+v)]
X = emb(u, v) + 0.01*rng.normal(size=(N, 6))      # same sample as neural.py (same seed and order)
dmin=[-1.17,-1.16,-1.15,-1.14,-0.63,-0.62]; dmax=[1.13,1.14,1.15,1.16,0.61,0.62]; div=[12,12,12,12,6,6]
M = PUIFast(X, dmin, dmax, div, 2); ev = M.njit_eval()
@njit
def newton(V, ev, its):
    n = V.shape[1]; y = np.zeros(4); y2 = np.zeros(4); J = np.zeros((4, n)); e = 1e-6
    for i in range(V.shape[0]):
        x = V[i].copy()
        for it in range(its):
            ev(x, y)
            for j in range(n):
                x2 = x.copy(); x2[j] += e; ev(x2, y2)
                J[:, j] = (y2 - y) / e
            x = x - np.linalg.pinv(J) @ y
        V[i] = x
n_, k_, D, cells = read_pol(REPO + '/book/runs/ml/torus6_mlp/results.pol')
for c in cells:
    newton(c.coords, ev, 2)
tt = np.linspace(0, 2*np.pi, 300); U, Vv = np.meshgrid(tt, tt); TG = cKDTree(emb(U.ravel(), Vv.ravel()))
h = (np.asarray(dmax)-np.asarray(dmin))/np.asarray(div)
P, d = project(X, cells, 2, np.asarray(dmin, float), h)
info = pickle.load(open(REPO + '/book/runs/ml/torus6_mlp_info.pkl','rb'))
info['err_after_newton'] = np.sqrt(np.mean(TG.query(P)[0]**2)); info['chi_newton'] = euler(cells)
allv = np.vstack([c.coords for c in cells]); info['vert_err_newton'] = np.sqrt(np.mean(TG.query(allv)[0]**2))
print(info['err_after_newton'], info['chi_newton'], info['vert_err_newton'])
pickle.dump(info, open(REPO + '/book/runs/ml/torus6_mlp_info.pkl','wb'))
