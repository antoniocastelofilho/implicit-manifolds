"""Learning the implicit function: a small neural network is trained to reproduce the partition of unity
implicit near the data (distillation), and its forward pass, compiled with numba, replaces F in the
continuation.  Its cost does not depend on the 3^n neighbouring cells visited by the PUI evaluation."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, time, pickle; sys.path.insert(0, REPO + '/python/implicit_manifolds/learning'); sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np
from numba import njit, cfunc, types, carray
from sklearn.neural_network import MLPRegressor
from pui_fast import PUIFast
from runfch import run_fch
from polfile import read_pol, euler
from project import project
rng = np.random.default_rng(7)
R = REPO + '/book/runs/ml/'


def mlp_njit(model, xm, xs, ys):
    Ws = [np.ascontiguousarray(w) for w in model.coefs_]; bs = [np.ascontiguousarray(b) for b in model.intercepts_]
    W0, W1, W2 = Ws; b0, b1, b2 = bs
    n = W0.shape[0]; k = W2.shape[1]

    @njit(nogil=True, fastmath=True)
    def ev(x, y):
        z = (x - xm) / xs
        h1 = np.tanh(z @ W0 + b0)
        h2 = np.tanh(h1 @ W1 + b1)
        o = h2 @ W2 + b2
        for j in range(k): y[j] = o[j] * ys[j]
    return ev


def cf(ev, n, k):
    sig = types.void(types.CPointer(types.double), types.CPointer(types.double))
    @cfunc(sig, nopython=True, fastmath=True)
    def shim(xp, yp):
        ev(carray(xp, n), carray(yp, k))
    return shim


def timing(ev, X, reps=1):
    y = np.zeros(64); ev(X[0], y[:1] if False else y)  # noqa
    return None


@njit
def bench(ev, Xs, k):
    y = np.zeros(k); s = 0.0
    for i in range(Xs.shape[0]):
        ev(Xs[i], y); s += y[0]
    return s


if __name__ == '__main__':
    N = 40000; u = rng.uniform(0, 2 * np.pi, N); v = rng.uniform(0, 2 * np.pi, N)
    def emb(u, v): return np.c_[np.cos(u), np.sin(u), np.cos(v), np.sin(v), 0.5 * np.cos(u + v), 0.5 * np.sin(u + v)]
    X = emb(u, v) + 0.01 * rng.normal(size=(N, 6))
    dmin = [-1.17, -1.16, -1.15, -1.14, -0.63, -0.62]; dmax = [1.13, 1.14, 1.15, 1.16, 0.61, 0.62]; div = [12, 12, 12, 12, 6, 6]
    M = PUIFast(X, dmin, dmax, div, 2)
    evp = M.njit_eval()
    # training set: samples around the data, targets from the PUI
    S = np.concatenate([X + s * rng.normal(size=X.shape) for s in (0.0, 0.03, 0.08, 0.2)])
    Y = np.zeros((len(S), 4)); y = np.zeros(4)
    for i in range(len(S)):
        evp(S[i], y); Y[i] = y
    xm = S.mean(0); xs = S.std(0); ys = np.abs(Y).std(0) + 1e-12
    t = time.time()
    model = MLPRegressor(hidden_layer_sizes=(96, 96), activation='tanh', max_iter=120, batch_size=512,
                         learning_rate_init=2e-3, tol=1e-7, n_iter_no_change=20, random_state=0)
    model.fit((S - xm) / xs, Y / ys)
    ttrain = time.time() - t
    evn = mlp_njit(model, xm, xs, ys)
    # accuracy near the manifold
    tt = np.linspace(0, 2 * np.pi, 300); U, V = np.meshgrid(tt, tt); G = emb(U.ravel(), V.ravel())
    yp = np.zeros(4); yn = np.zeros(4); ep = []; en = []
    for x in G[::7]:
        evp(x, yp); evn(x, yn); ep.append(np.linalg.norm(yp)); en.append(np.linalg.norm(yn))
    # cost of one evaluation
    Q = X[:20000] + 0.05 * rng.normal(size=(20000, 6))
    bench(evp, Q[:10], 4); bench(evn, Q[:10], 4)
    t = time.time(); bench(evp, Q, 4); tp = (time.time() - t) / len(Q)
    t = time.time(); bench(evn, Q, 4); tn = (time.time() - t) / len(Q)
    print('train %.1fs  |F| on M: pui %.4f mlp %.4f  eval pui %.2e s  mlp %.2e s' % (ttrain, np.mean(ep), np.mean(en), tp, tn), flush=True)
    # continuation with the learned function
    shim = cf(evn, 6, 4)
    vals = []
    for x in X[:300]:
        evn(x, yn); vals.append(np.linalg.norm(yn))
    x0 = X[int(np.argmin(vals))]
    pol, tm = run_fch(shim, 6, 4, dmin, dmax, div, x0, 'torus6_mlp')
    n_, k_, D, cells = read_pol(pol); chi = euler(cells)
    from scipy.spatial import cKDTree
    TG = cKDTree(emb(U.ravel(), V.ravel()))
    h = (np.asarray(dmax) - np.asarray(dmin)) / np.asarray(div)
    P, d = project(X, cells, 2, np.asarray(dmin, float), h)
    eb = TG.query(X)[0]; ea = TG.query(P)[0]
    out = dict(train=ttrain, loss=model.loss_, iters=model.n_iter_, Fm_pui=np.mean(ep), Fm_mlp=np.mean(en), eval_pui=tp, eval_mlp=tn,
               fch=tm, cells=len(cells), chi=chi, err_before=np.sqrt(np.mean(eb ** 2)), err_after=np.sqrt(np.mean(ea ** 2)))
    print(out)
    pickle.dump(out, open(R + 'torus6_mlp_info.pkl', 'wb'))
