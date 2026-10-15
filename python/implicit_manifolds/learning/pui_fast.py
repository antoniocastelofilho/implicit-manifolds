"""Parallel construction of the partition of unity implicit (same model as pui.PUI).

All the local approximations are independent: the points are binned, every point is assigned to the
3^n cells whose supports contain it, the (cell, point) pairs are sorted by cell, and the local PCA and
least-squares fits are computed for all cells in parallel (numba prange).  The alignment of the normal
frames is a breadth-first traversal on the sorted cell keys.  The same batched structure (many small
covariance matrices, eigendecompositions and least-squares systems) maps directly to GPU batched
linear algebra."""
import numpy as np, itertools, time
from numba import njit, prange
from pui import PUI


@njit(cache=True)
def _bsp(t):
    t = abs(t)
    if t <= 0.5: return 0.75 - t * t
    if t < 1.5: return 0.5 * (1.5 - t) ** 2
    return 0.0


@njit(parallel=True, cache=True)
def _local_fits(X, pidx, starts, ckeys, div, strides, dmin, h, m, nc, minpts, flat_max, tau_min, lam):
    C = starts.shape[0] - 1; n = X.shape[1]; k = n - m
    P = np.zeros((C, n)); T = np.zeros((C, n, m)); N = np.zeros((C, n, k)); coef = np.zeros((C, k, nc))
    eig = np.zeros((C, n)); ok = np.zeros(C, np.bool_)
    for c in prange(C):
        a = starts[c]; b = starts[c + 1]; cnt = b - a
        if cnt < minpts: continue
        # centre of the cell
        ctr = np.empty(n); key = ckeys[c]
        for i in range(n - 1, -1, -1):
            q = key // strides[i]; key -= q * strides[i]; ctr[i] = dmin[i] + (q + 0.5) * h[i]
        w = np.empty(cnt); W = 0.0
        p = np.zeros(n)
        for j in range(cnt):
            x = X[pidx[a + j]]; ww = 1.0
            for i in range(n): ww *= _bsp((x[i] - ctr[i]) / h[i])
            w[j] = ww; W += ww
            for i in range(n): p[i] += ww * x[i]
        if W <= 0: continue
        p /= W
        Cv = np.zeros((n, n))
        for j in range(cnt):
            d = X[pidx[a + j]] - p
            for r in range(n):
                for s in range(n): Cv[r, s] += w[j] * d[r] * d[s]
        Cv /= W
        lam_, V = np.linalg.eigh(Cv)
        if lam_[k - 1] > flat_max * lam_[k]: continue
        Tc = V[:, k:][:, ::-1].copy(); Nc = V[:, :k].copy()
        # local coordinates and coverage
        U = np.empty((cnt, m)); Vn = np.empty((cnt, k))
        for j in range(cnt):
            d = X[pidx[a + j]] - p
            for t in range(m):
                s = 0.0
                for i in range(n): s += Tc[i, t] * d[i]
                U[j, t] = s
            for t in range(k):
                s = 0.0
                for i in range(n): s += Nc[i, t] * d[i]
                Vn[j, t] = s
        area = 1.0; hs = 0.0
        for i in range(n): hs += h[i] * h[i]
        sup = 3.0 * np.sqrt(hs / n)
        for t in range(m):
            lo = 1e300; hi = -1e300
            for j in range(cnt):
                if w[j] > 0:
                    lo = min(lo, U[j, t]); hi = max(hi, U[j, t])
            area *= max(hi - lo, 0.0) / sup
        deg2 = area >= tau_min and nc > 1 + m
        ncu = nc if deg2 else 1 + m
        A = np.empty((cnt, ncu))
        for j in range(cnt):
            A[j, 0] = 1.0
            for t in range(m): A[j, 1 + t] = U[j, t]
            if deg2:
                q = 1 + m
                for t1 in range(m):
                    for t2 in range(t1, m):
                        A[j, q] = U[j, t1] * U[j, t2]; q += 1
        G = np.zeros((ncu, ncu)); R = np.zeros((ncu, k))
        for j in range(cnt):
            for r in range(ncu):
                for s in range(ncu): G[r, s] += w[j] * A[j, r] * A[j, s]
                for t in range(k): R[r, t] += w[j] * A[j, r] * Vn[j, t]
        for r in range(ncu): G[r, r] += lam
        sol = np.linalg.solve(G, R)
        for t in range(k):
            for r in range(ncu): coef[c, t, r] = sol[r, t]
        P[c] = p; T[c] = Tc; N[c] = Nc; ok[c] = True
        for i in range(n): eig[c, i] = lam_[n - 1 - i]
    return P, T, N, coef, eig, ok


@njit(cache=True)
def _align(keys, N, coef, offl, valid_off, div, strides):
    """BFS over the cells (sorted keys); valid_off[c, o] tells whether neighbour o exists inside the grid"""
    C = keys.shape[0]; k = N.shape[2]
    done = np.zeros(C, np.bool_); queue = np.empty(C, np.int64)
    for s0 in range(C):
        if done[s0]: continue
        done[s0] = True; head = 0; tail = 0; queue[tail] = s0; tail += 1
        while head < tail:
            a = queue[head]; head += 1
            for o in range(offl.shape[0]):
                if not valid_off[a, o]: continue
                kb = keys[a] + offl[o]
                b = np.searchsorted(keys, kb)
                if b >= C or keys[b] != kb or done[b]: continue
                M = np.zeros((k, k))
                for o2 in range(offl.shape[0]):
                    if not valid_off[b, o2]: continue
                    kc = keys[b] + offl[o2]
                    c = np.searchsorted(keys, kc)
                    if c < C and keys[c] == kc and done[c]:
                        M += N[b].T @ N[c]
                U, S, Vt = np.linalg.svd(M)
                Rr = U @ Vt
                N[b] = N[b] @ Rr
                coef[b] = Rr.T @ coef[b]          # the normal coordinates rotate with the frame
                done[b] = True; queue[tail] = b; tail += 1


class PUIFast(PUI):
    def __init__(self, X, dmin, dmax, div, m, minpts=None, deg=2, lam=1e-8, fill=True, tau_min=0.4, flat_max=0.3, verbose=False):
        t0 = time.time()
        X = np.ascontiguousarray(X, float); self.X = X
        self.n = n = X.shape[1]; self.m = m; self.k = k = n - m
        self.dmin = np.asarray(dmin, float); self.dmax = np.asarray(dmax, float)
        self.div = np.asarray(div, np.int64); self.h = (self.dmax - self.dmin) / self.div
        self.deg = deg; nc = 1 + m + (m * (m + 1) // 2 if deg == 2 else 0); self.ncoef = nc
        minpts = minpts or max(3 * nc, m + 2)
        strides = np.cumprod(np.r_[1, self.div[:-1]]).astype(np.int64); self.strides = strides
        I = np.floor((X - self.dmin) / self.h).astype(np.int64)
        keep = np.all((I >= 0) & (I < self.div), 1); X = X[keep]; I = I[keep]
        offs = np.array(list(itertools.product((-1, 0, 1), repeat=n)), np.int64); self.offs = offs
        # (cell, point) pairs, in chunks of offsets to bound memory
        keys = []; pts = []
        for o in offs:
            J = I + o; v = np.all((J >= 0) & (J < self.div), 1)
            keys.append(J[v] @ strides); pts.append(np.nonzero(v)[0])
        keys = np.concatenate(keys); pts = np.concatenate(pts)
        order = np.argsort(keys, kind='stable'); keys = keys[order]; pts = pts[order]
        ck, starts = np.unique(keys, return_index=True); starts = np.r_[starts, len(keys)]
        t1 = time.time()
        P, T, N, coef, eig, ok = _local_fits(X, pts, starts, ck, self.div, strides, self.dmin, self.h, m, nc, minpts, flat_max, tau_min, lam)
        t2 = time.time()
        ck = ck[ok]; P = P[ok]; T = T[ok]; N = N[ok]; coef = coef[ok]; eig = eig[ok]
        # alignment
        multi = np.stack([(ck // strides[i]) % self.div[i] for i in range(n)], 1)
        valid = np.ones((len(ck), len(offs)), np.bool_)
        for oi, o in enumerate(offs):
            J = multi + o; valid[:, oi] = np.all((J >= 0) & (J < self.div), 1)
        _align(ck, N, coef, offs @ strides, valid, self.div, strides)
        t3 = time.time()
        self.cells = [tuple(r) for r in multi]; self.P = P; self.T = T; self.N = N; self.coef = coef; self.eig = eig
        self.ndata = len(ck); self.isfill = [False] * len(ck)
        self.index = {c: i for i, c in enumerate(self.cells)}
        if fill: self._fill()
        self._build_arrays()
        t4 = time.time()
        self.timings = dict(pairs=t1 - t0, fits=t2 - t1, align=t3 - t2, fill=t4 - t3, total=t4 - t0)
        if verbose: print(self.timings)

    def _fill(self):
        from scipy.spatial import cKDTree
        n = self.n; strides = self.strides
        multi = np.array(self.cells, np.int64); ck = multi @ strides
        if np.prod(self.div) <= 400000:
            allk = np.arange(int(np.prod(self.div)), dtype=np.int64)
        else:
            cand = []
            for o in self.offs:
                J = multi + o; v = np.all((J >= 0) & (J < self.div), 1); cand.append(J[v] @ strides)
            allk = np.unique(np.concatenate(cand))
        allk = np.setdiff1d(allk, ck)
        if len(allk) == 0: return
        am = np.stack([(allk // strides[i]) % self.div[i] for i in range(n)], 1)
        _, nb = cKDTree(multi.astype(float)).query(am.astype(float))
        cf = self.coef[nb].copy(); cf[:, :, 1 + self.m:] = 0.0
        self.cells = self.cells + [tuple(r) for r in am]
        self.P = np.concatenate([self.P, self.P[nb]]); self.T = np.concatenate([self.T, self.T[nb]])
        self.N = np.concatenate([self.N, self.N[nb]]); self.coef = np.concatenate([self.coef, cf])
        self.isfill = self.isfill + [True] * len(allk)
        self.index = None

    def _build_arrays(self):
        self.grid = -np.ones(int(np.prod(self.div)), np.int64)
        keys = np.array(self.cells, np.int64) @ self.strides
        self.grid[keys] = np.arange(len(keys))
