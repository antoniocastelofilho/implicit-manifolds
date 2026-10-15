"""Radial basis function (RBF) and Hermite RBF (HRBF) local approximations inside the partition of unity.

The frames (p, T, N), the alignment, the robustness criteria and the fill cells are those of PUIFast.
In each data cell the normal coordinates v = N^T (x - p) are approximated, as functions of the tangent
coordinates u = T^T (x - p), by

    s(u) = q(u) + sum_j a_j phi(u - c_j)                                   (kind='rbf')
    s(u) = q(u) + sum_j [ a_j phi(u - c_j) - sum_l b_jl d_l phi(u - c_j) ]   (kind='hrbf')

with phi(u) = |u|^3 (polyharmonic, conditionally positive definite of order 2), q the weighted least-squares
quadric of PUIFast and centres c_j chosen by farthest-point sampling among the tangent coordinates of the
points of the block.  The RBF part fits the residual v - q(u) on all the points of the block (weighted least
squares) and, for the Hermite version, also the residual slopes J_i - Dq(u_i), where J_i = dv/du is
estimated at each sample point from the tangent space of its k nearest neighbours (no normals are given).
The coefficients satisfy the side conditions P^T gamma = 0 (orthogonality to linear polynomials), and the
fit is regularized by the native-space semi-norm gamma^T K gamma, K the (Hermite) interpolation matrix of the
centres: a smoothing (Hermite) RBF in the sense of Wahba/Macedo-Gois-Velho."""
import time, numpy as np
from numba import njit, cfunc, types, carray
from scipy.spatial import cKDTree
from scipy.linalg import null_space
from pui_fast import PUIFast, _bsp


def phi_basis(U, C, hermite):
    """columns: phi(U - c_j) and, if hermite, -d_l phi(U - c_j); phi = r^3, grad = 3 r d"""
    D = U[:, None, :] - C[None, :, :]; r = np.linalg.norm(D, axis=2)
    cols = [r ** 3]
    if hermite:
        cols += [-3.0 * r * D[:, :, l] for l in range(C.shape[1])]
    return np.concatenate(cols, 1)


def phi_basis_grad(U, C, hermite):
    """derivative d/du_a of each column, returned as (len(U), m, ncols)"""
    D = U[:, None, :] - C[None, :, :]; r = np.linalg.norm(D, axis=2); m = C.shape[1]
    out = []
    for a in range(m):
        cols = [3.0 * r * D[:, :, a]]
        if hermite:
            rs = np.where(r > 0, r, 1.0)
            for l in range(m):
                cols.append(-3.0 * (D[:, :, l] * D[:, :, a] / rs + (r if l == a else 0.0)))
        out.append(np.concatenate(cols, 1))
    return np.stack(out, 1)


def gram(C, hermite):
    """(Hermite) interpolation matrix K of the centres and side-condition matrix P"""
    M, m = C.shape
    if not hermite:
        K = phi_basis(C, C, False)
        P = np.c_[np.ones(M), C]
        return K, P
    V = phi_basis(C, C, True)                         # value rows
    G = phi_basis_grad(C, C, True)                    # gradient rows (M, m, cols)
    K = np.concatenate([V] + [G[:, a, :] for a in range(m)], 0)
    P = np.concatenate([np.c_[np.ones(M), C]] + [np.c_[np.zeros(M), np.tile(np.eye(m)[a], (M, 1))] for a in range(m)], 0)
    # column ordering of K: [phi_j | -d_1 phi_j | ... ], row ordering the same blocks: symmetric
    return 0.5 * (K + K.T), P


def fps(U, M, rng):
    i0 = int(rng.integers(len(U))); sel = [i0]; d = np.linalg.norm(U - U[i0], axis=1)
    for _ in range(M - 1):
        i = int(np.argmax(d)); sel.append(i); d = np.minimum(d, np.linalg.norm(U - U[i], axis=1))
    return np.array(sel)


def point_tangents(X, m, knn=20):
    """tangent basis of every sample point from the PCA of its knn nearest neighbours"""
    _, nb = cKDTree(X).query(X, knn)
    Y = X[nb] - X[nb].mean(1, keepdims=True)
    C = np.einsum('pji,pjk->pik', Y, Y)
    _, V = np.linalg.eigh(C)
    return V[:, :, -m:]


class PUIRBF(PUIFast):
    def __init__(self, X, dmin, dmax, div, m, kind='rbf', ncent=16, reg=1e-4, eta=1.0, knn=20, seed=0, tau_rbf=0.8, **kw):
        super().__init__(X, dmin, dmax, div, m, **kw)
        self.kind = kind; hermite = kind == 'hrbf'
        t0 = time.time()
        X = self.X; n, k = self.n, self.k
        nb = ncent * (m + 1 if hermite else 1); self.nbasis = nb; self.ncent = ncent
        C_all = np.zeros((len(self.cells), ncent, m)); G_all = np.zeros((len(self.cells), k, nb))
        has = np.zeros(len(self.cells), np.int64)
        TP = point_tangents(X, m, knn) if hermite else None
        I = np.floor((X - self.dmin) / self.h).astype(np.int64)
        keep = np.all((I >= 0) & (I < self.div), 1)
        bins = {}
        for i in np.nonzero(keep)[0]:
            bins.setdefault(int(I[i] @ self.strides), []).append(i)
        offk = self.offs @ self.strides
        rng = np.random.default_rng(seed)
        self.rbf_cells = 0
        for c in range(self.ndata):
            cm = np.array(self.cells[c]); key = int(cm @ self.strides)
            ids = []
            for o, ok_ in zip(self.offs, offk):
                J = cm + o
                if np.all((J >= 0) & (J < self.div)): ids += bins.get(key + int(ok_), [])
            ids = np.array(ids)
            Y = X[ids]; ctr = self.dmin + (cm + 0.5) * self.h
            w = np.prod(_bsp_np((Y - ctr) / self.h), 1); sel = w > 0
            Y, w, ids = Y[sel], w[sel], ids[sel]
            if len(Y) < 3 * nb: continue
            D = Y - self.P[c]; U = D @ self.T[c]; Vn = D @ self.N[c]
            # coverage: radial functions extrapolate badly where the block is not covered by the sample
            tau = np.prod(np.ptp(U, 0)) / (3.0 * np.linalg.norm(self.h) / np.sqrt(n)) ** m
            if tau < tau_rbf: continue
            q, dq = self._poly(c, U)
            Rv = Vn - q                                       # residual values (N, k)
            cen = U[fps(U, ncent, rng)]
            B = phi_basis(U, cen, hermite)
            K, Pm = gram(cen, hermite)
            Z = null_space(Pm.T)                               # gamma = Z beta
            W = np.sqrt(w)[:, None]
            A = (B @ Z) * W; rhs = Rv * W
            if hermite:
                Tc = self.T[c].T @ TP[ids]; Nc = self.N[c].T @ TP[ids]          # (N, m, m), (N, k, m)
                good = np.abs(np.linalg.det(Tc)) > 0.5
                Jv = np.einsum('pkm,pmj->pkj', Nc[good], np.linalg.inv(Tc[good]))   # dv/du (N, k, m)
                Rd = Jv - dq[good]                                               # residual slopes
                Gb = phi_basis_grad(U[good], cen, True)                          # (N, m, nb)
                s = eta * np.mean(self.h) * W[good][:, :, None]
                A = np.concatenate([A] + [(Gb[:, a, :] @ Z) * s[:, :, 0] for a in range(m)], 0)
                rhs = np.concatenate([rhs] + [Rd[:, :, a] * s[:, :, 0] for a in range(m)], 0)
            H = A.T @ A; Kz = Z.T @ K @ Z
            lam = reg * np.trace(H) / max(np.trace(Kz), 1e-300)
            beta = np.linalg.solve(H + lam * Kz + 1e-12 * np.trace(H) * np.eye(len(H)), A.T @ rhs)
            C_all[c] = cen; G_all[c] = (Z @ beta).T; has[c] = 1; self.rbf_cells += 1
        self.cen = C_all; self.gam = G_all; self.has = has
        self.timings['rbf'] = time.time() - t0

    def _poly(self, c, U):
        m = self.m; cf = self.coef[c]                      # (k, nc)
        A = [np.ones(len(U))] + [U[:, a] for a in range(m)]
        dA = [np.zeros((len(U), m))] + [np.tile(np.eye(m)[a], (len(U), 1)) for a in range(m)]
        if self.deg == 2:
            for a in range(m):
                for b in range(a, m):
                    A.append(U[:, a] * U[:, b]); g = np.zeros((len(U), m)); g[:, a] += U[:, b]; g[:, b] += U[:, a]; dA.append(g)
        A = np.stack(A, 1); dA = np.stack(dA, 2)            # (N, nc), (N, m, nc)
        return A @ cf.T, np.einsum('pmt,kt->pkm', dA, cf)

    def njit_eval(self):
        n, m, k, nc = self.n, self.m, self.k, self.ncoef
        dmin, h, div, strides = self.dmin.copy(), self.h.copy(), self.div.copy(), self.strides.copy()
        grid, P, T, N, coef, offs = self.grid, self.P, self.T, self.N, self.coef, self.offs.astype(np.int64)
        deg = self.deg; cen, gam, has = self.cen, self.gam, self.has
        hermite = self.kind == 'hrbf'; M = self.ncent

        @njit(nogil=True, fastmath=True)
        def ev(x, y):
            base = np.empty(n, np.int64)
            for i in range(n): base[i] = int(np.floor((x[i] - dmin[i]) / h[i]))
            W = 0.0; acc = np.zeros(k)
            u = np.empty(m); d = np.empty(n); terms = np.empty(nc); dd = np.empty(m)
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
                if c < has.shape[0] and has[c] == 1:
                    for jc in range(M):
                        r2 = 0.0
                        for a in range(m):
                            dd[a] = u[a] - cen[c, jc, a]; r2 += dd[a] * dd[a]
                        r = np.sqrt(r2); p3 = r2 * r
                        for j in range(k):
                            s = gam[c, j, jc] * p3
                            if hermite:
                                for a in range(m): s -= gam[c, j, M * (1 + a) + jc] * 3.0 * r * dd[a]
                            acc[j] -= w * s
                W += w
            if W < 1e-12:
                for j in range(k): y[j] = 1.0
            else:
                for j in range(k): y[j] = acc[j] / W
        return ev


def _bsp_np(t):
    t = np.abs(t)
    return np.where(t <= 0.5, 0.75 - t * t, np.where(t < 1.5, 0.5 * (1.5 - t) ** 2, 0.0))
