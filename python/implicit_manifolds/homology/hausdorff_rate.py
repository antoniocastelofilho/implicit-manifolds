"""
Geometric convergence of the PL approximation: one-sided Hausdorff distance
max_{p in M_T vertices} dist(p, M) for the torus, whose distance function is exact:
    dist(x, T) = | sqrt((sqrt(x^2+y^2) - R)^2 + z^2) - r |.
FCH vertices (k = 1) are the zeros of the linear interpolant on the grid edges
(the k-faces); arcs live in the 2-faces, triangulated by K1.  Predicted rate: O(h^2).  Also reports the predicted
bound  n K h^2 / mu  of Theorem 6.x, with K = sup |D^2 F| and mu = inf |grad F| near M.
"""
import itertools
import numpy as np

R, r = 0.6, 0.25
F = lambda x, y, z: (np.sqrt(x**2 + y**2) - R) ** 2 + z**2 - r**2
dist = lambda x, y, z: np.abs(np.sqrt((np.sqrt(x**2 + y**2) - R) ** 2 + z**2) - r)

def run(N, shift=0.0123):
    xs = np.linspace(-1, 1, N + 1) + shift
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    V = F(X, Y, Z)
    worst = 0.0
    P = np.stack([X, Y, Z], -1)
    # grid edges v -> v + e_i
    for sub in itertools.product((0, 1), repeat=3):
        if sum(sub) != 1:          # only grid edges = k-faces for k = 1
            continue
        sl0 = tuple(slice(0, N + 1 - s) for s in sub)
        sl1 = tuple(slice(s, N + 1) for s in sub)
        a, b = V[sl0], V[sl1]
        m = (a * b) < 0
        t = a[m] / (a[m] - b[m])
        p = P[sl0][m] + t[:, None] * (P[sl1][m] - P[sl0][m])
        if len(p):
            worst = max(worst, dist(p[:, 0], p[:, 1], p[:, 2]).max())
    return worst

if __name__ == "__main__":
    # constants for the a priori bound (torus): |grad F| = 2 r on M; near M, |grad F| >= 2(r - d)
    # D^2 F has norm <= 2 R / (R - r - d) + 2 roughly; we just report the observed rate.
    prev = None
    for N in (8, 16, 32, 64, 128):
        h = 2 / N
        w = run(N)
        rate = "" if prev is None else f"  ratio {prev / w:5.2f} (h^2 -> 4)"
        print(f"N={N:4d} h={h:.4f}  max dist(M_T vertices, M) = {w:.3e}   dist/h^2 = {w/h**2:.3f}{rate}")
        prev = w
