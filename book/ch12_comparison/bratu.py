"""Discrete Bratu problem: (u_{i-1}-2u_i+u_{i+1})/h^2 + lam*exp(u_i) = 0, u_0=u_{N+1}=0.  Unknowns x=(lam,u_1..u_N)."""
import numpy as np
def F(x, N):
    lam = x[0]; u = np.r_[0, x[1:], 0]; h = 1.0 / (N + 1)
    return (u[:-2] - 2 * u[1:-1] + u[2:]) / h ** 2 + lam * np.exp(u[1:-1])
def J(x, N):
    lam = x[0]; u = x[1:]; h = 1.0 / (N + 1)
    A = (np.diag(-2 * np.ones(N)) + np.diag(np.ones(N - 1), 1) + np.diag(np.ones(N - 1), -1)) / h ** 2 + lam * np.diag(np.exp(u))
    return np.c_[np.exp(u), A]
def natural(N, lams):
    """natural-parameter continuation: Newton in u for fixed lam, previous solution as guess"""
    u = np.zeros(N); out = []
    for lam in lams:
        ok = False
        for it in range(30):
            x = np.r_[lam, u]; r = F(x, N)
            if np.linalg.norm(r) < 1e-10: ok = True; break
            u = u - np.linalg.solve(J(x, N)[:, 1:], r)
            if not np.all(np.isfinite(u)) or np.abs(u).max() > 50: break
        if not ok: return np.array(out), lam
        out.append(np.r_[lam, u])
    return np.array(out), None
def arclength(N, ds=0.02, smax=40, umax=8.0):
    """pseudo-arclength continuation for the reference curve"""
    x = np.zeros(N + 1); t = np.zeros(N + 1); t[0] = 1; P = [x.copy()]
    for _ in range(int(smax / ds)):
        xp = x + ds * t
        for it in range(20):
            r = np.r_[F(xp, N), t @ (xp - x) - ds]
            Jx = np.r_[J(xp, N), t[None, :]]
            d = np.linalg.solve(Jx, r); xp = xp - d
            if np.linalg.norm(d) < 1e-12: break
        Jx = J(xp, N); _, _, Vt = np.linalg.svd(Jx); tn = Vt[-1]
        if tn @ t < 0: tn = -tn
        x, t = xp, tn; P.append(x.copy())
        if x[1:].max() > umax or x[0] < 0: break
    return np.array(P)
if __name__ == '__main__':
    for N in (3, 5):
        P = arclength(N); i = np.argmax(P[:, 0]); print(N, 'fold lam', P[i, 0], 'u', P[i, 1:].round(3), 'end', P[-1].round(2), len(P))
        nat, fail = natural(N, np.arange(0, 4.0, 0.05)); print(' natural fails at', fail, 'last', nat[-1].round(3))
        print(' u ranges', P[:, 1:].min(0).round(3), P[:, 1:].max(0).round(3))
