"""Partition of unity implicits (PUI) for point clouds sampled near an m-manifold of R^n.

Local approximations are computed in each cell of a grid from the points of the cell and of its 3^n
neighbours: local PCA gives a tangent basis T (n x m) and a normal basis N (n x k), and the k normal
coordinates are fitted by quadratic polynomials of the m tangent coordinates.  The local implicit
function is f_c(x) = N^T (x - p) - q(T^T (x - p)) in R^k.  Normal frames are aligned between
neighbouring cells by orthogonal Procrustes, and the global implicit function is the partition of unity
F(x) = sum_c w_c(x) f_c(x) / sum_c w_c(x), with tensor-product quadratic B-spline weights."""
import numpy as np, itertools, ctypes
from collections import deque
from numba import njit, cfunc, types, carray


def bspline2(t):
    t = np.abs(t)
    return np.where(t <= 0.5, 0.75 - t * t, np.where(t < 1.5, 0.5 * (1.5 - t) ** 2, 0.0))


def quad_terms(U):
    """[1, u_i, u_i u_j (i<=j)] for rows of U (N x m)"""
    m = U.shape[1]
    cols = [np.ones(len(U))] + [U[:, i] for i in range(m)]
    cols += [U[:, i] * U[:, j] for i in range(m) for j in range(i, m)]
    return np.stack(cols, 1)


class PUI:
    def __init__(self, X, dmin, dmax, div, m, minpts=None, deg=2, lam=1e-8, fill=True, knn=20, tau_min=0.4, flat_max=0.3):
        X = np.asarray(X, float); self.X = X
        self.n = n = X.shape[1]; self.m = m; self.k = k = n - m
        self.dmin = np.asarray(dmin, float); self.dmax = np.asarray(dmax, float)
        self.div = np.asarray(div, int); self.h = (self.dmax - self.dmin) / self.div
        self.deg = deg
        nc = 1 + m + (m * (m + 1) // 2 if deg == 2 else 0)
        self.ncoef = nc
        minpts = minpts or max(3 * nc, m + 2)
        # bin points
        idx = np.floor((X - self.dmin) / self.h).astype(int)
        ok = np.all((idx >= 0) & (idx < self.div), 1)
        X, idx = X[ok], idx[ok]
        bins = {}
        for i, key in enumerate(map(tuple, idx)):
            bins.setdefault(key, []).append(i)
        offs = np.array(list(itertools.product((-1, 0, 1), repeat=n)))
        self.offs = offs
        # candidate cells: every cell with points in its 3^n block
        cand = set()
        for key in bins:
            for o in offs:
                c = tuple(np.add(key, o))
                if all(0 <= c[i] < self.div[i] for i in range(n)): cand.add(c)
        cells = []; P = []; T = []; Nm = []; ev = []; sup = []
        for c in sorted(cand):
            pts = [j for o in offs for j in bins.get(tuple(np.add(c, o)), [])]
            if len(pts) < minpts: continue
            Y = X[pts]; ctr = self.dmin + (np.array(c) + 0.5) * self.h
            w = np.prod(bspline2((Y - ctr) / self.h), 1)
            if w.sum() <= 0: continue
            p = (w[:, None] * Y).sum(0) / w.sum()
            Cv = ((Y - p) * w[:, None]).T @ (Y - p) / w.sum()
            lam_, V = np.linalg.eigh(Cv)                 # ascending
            # flatness test: the support must look like a piece of an m-dimensional manifold
            # (largest normal variance small with respect to the smallest tangent variance)
            if lam_[k - 1] > flat_max * lam_[k]: continue
            cells.append(c); P.append(p); T.append(V[:, k:][:, ::-1]); Nm.append(V[:, :k]); ev.append(lam_[::-1]); sup.append(pts)
        self.ndata = len(cells)
        self.isfill = [False] * len(cells)
        self.cells = cells; self.P = np.array(P); self.T = np.array(T); self.N = np.array(Nm)
        self.eig = np.array(ev[:self.ndata]); self.support = sup
        self.index = {c: i for i, c in enumerate(cells)}
        self._align()
        self.tau_min = tau_min
        self._fit(X, lam)
        if fill: self._fill()
        self._build_arrays()

    def _align(self):
        """align normal frames by Procrustes along a breadth-first traversal of the cells"""
        C = len(self.cells); done = np.zeros(C, bool)
        for s in range(C):
            if done[s]: continue
            done[s] = True; Q = deque([s])
            while Q:
                a = Q.popleft()
                for o in self.offs:
                    b = self.index.get(tuple(np.add(self.cells[a], o)))
                    if b is None or done[b]: continue
                    # rotation R (k x k) minimizing |N_b R - (frames of aligned neighbours of b)|
                    M = np.zeros((self.k, self.k))
                    for o2 in self.offs:
                        c = self.index.get(tuple(np.add(self.cells[b], o2)))
                        if c is not None and done[c]: M += self.N[b].T @ self.N[c]
                    U, S, Vt = np.linalg.svd(M)
                    self.N[b] = self.N[b] @ (U @ Vt)
                    done[b] = True; Q.append(b)


    def _fill(self):
        """every other cell of the grid receives the plane (linear part) of the nearest data cell, so that
        F is defined in the whole box and keeps the orientation of the data cells away from the sample"""
        from scipy.spatial import cKDTree
        have = set(self.cells)
        dc = np.array(self.cells, float)
        tree = cKDTree(dc)
        if np.prod(self.div) <= 400000:
            allc = [c for c in itertools.product(*[range(d) for d in self.div]) if c not in have]
        else:                                   # high dimension: fill a band of one cell around the data cells
            band = set()
            for c in self.cells:
                for o in self.offs:
                    b = tuple(np.add(c, o))
                    if b not in have and all(0 <= b[i] < self.div[i] for i in range(self.n)): band.add(b)
            allc = sorted(band)
        if not allc: return
        _, nb = tree.query(np.array(allc, float))
        m = self.m
        P = [self.P]; T = [self.T]; N = [self.N]; coef = [self.coef]
        cf = self.coef[nb].copy(); cf[:, :, 1 + m:] = 0.0
        self.cells = self.cells + allc
        self.P = np.concatenate([self.P, self.P[nb]]); self.T = np.concatenate([self.T, self.T[nb]])
        self.N = np.concatenate([self.N, self.N[nb]]); self.coef = np.concatenate([self.coef, cf])
        self.isfill = self.isfill + [True] * len(allc)
        self.index = {c: i for i, c in enumerate(self.cells)}

    def _fit(self, X, lam):
        C = len(self.cells); k = self.k
        self.coef = np.zeros((C, k, self.ncoef))
        self.deg_used = np.zeros(C, int)
        for i, c in enumerate(self.cells):
            Y = X[self.support[i]]; ctr = self.dmin + (np.array(c) + 0.5) * self.h
            w = np.prod(bspline2((Y - ctr) / self.h), 1)
            D = Y - self.P[i]; U = D @ self.T[i]; Vn = D @ self.N[i]
            # coverage criterion: area of the box of the projected points relative to the projected support
            ext = np.ptp(U[w > 0], 0) if np.any(w > 0) else np.zeros(self.m)
            tau = np.prod(ext) / np.prod(np.full(self.m, 3.0 * np.linalg.norm(self.h) / np.sqrt(self.n)))
            deg = self.deg if tau >= self.tau_min else 1
            self.deg_used[i] = deg
            A = quad_terms(U) if deg == 2 else np.c_[np.ones(len(U)), U]
            Wh = np.sqrt(w)[:, None]
            G = (A * Wh).T @ (A * Wh) + lam * np.eye(A.shape[1])
            cf = np.linalg.solve(G, (A * Wh).T @ (Vn * Wh)).T
            self.coef[i, :, :cf.shape[1]] = cf

    def _build_arrays(self):
        n = self.n
        self.grid = -np.ones(int(np.prod(self.div)), np.int64)
        strides = np.cumprod(np.r_[1, self.div[:-1]]).astype(np.int64); self.strides = strides
        for i, c in enumerate(self.cells): self.grid[int(np.dot(c, strides))] = i

    def njit_eval(self):
        n, m, k, nc = self.n, self.m, self.k, self.ncoef
        dmin, h, div, strides = self.dmin.copy(), self.h.copy(), self.div.copy(), self.strides.copy()
        grid, P, T, N, coef, offs = self.grid, self.P, self.T, self.N, self.coef, self.offs.astype(np.int64)
        deg = self.deg

        @njit(nogil=True, fastmath=True)
        def ev(x, y):
            base = np.empty(n, np.int64)
            for i in range(n): base[i] = int(np.floor((x[i] - dmin[i]) / h[i]))
            W = 0.0
            acc = np.zeros(k)
            u = np.empty(m); d = np.empty(n); terms = np.empty(nc)
            for o in range(offs.shape[0]):
                key = 0; inside = True
                for i in range(n):
                    ci = base[i] + offs[o, i]
                    if ci < 0 or ci >= div[i]:
                        inside = False; break
                    key += ci * strides[i]
                if not inside: continue
                c = grid[key]
                if c < 0: continue
                w = 1.0
                for i in range(n):
                    t = abs((x[i] - dmin[i]) / h[i] - (base[i] + offs[o, i] + 0.5))
                    if t <= 0.5: w *= 0.75 - t * t
                    elif t < 1.5: w *= 0.5 * (1.5 - t) ** 2
                    else:
                        w = 0.0; break
                if w == 0.0: continue
                for i in range(n): d[i] = x[i] - P[c, i]
                for a in range(m):
                    s = 0.0
                    for i in range(n): s += T[c, i, a] * d[i]
                    u[a] = s
                terms[0] = 1.0
                for a in range(m): terms[1 + a] = u[a]
                if deg == 2:
                    q = 1 + m
                    for a in range(m):
                        for b in range(a, m):
                            terms[q] = u[a] * u[b]; q += 1
                for j in range(k):
                    v = 0.0
                    for i in range(n): v += N[c, i, j] * d[i]
                    qv = 0.0
                    for t in range(nc): qv += coef[c, j, t] * terms[t]
                    acc[j] += w * (v - qv)
                W += w
            if W < 1e-12:
                for j in range(k): y[j] = 1.0
            else:
                for j in range(k): y[j] = acc[j] / W
        return ev

    def cfunc(self):
        ev = self.njit_eval(); n, k = self.n, self.k
        sig = types.void(types.CPointer(types.double), types.CPointer(types.double))

        @cfunc(sig, nopython=True, fastmath=True)
        def shim(xp, yp):
            ev(carray(xp, n), carray(yp, k))
        return shim, ev
