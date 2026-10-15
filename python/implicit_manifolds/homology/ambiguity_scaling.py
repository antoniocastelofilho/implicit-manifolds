"""
Does grid refinement remove ambiguity?  (k = 1)

Heuristic: a face with j fixed directions is ambiguous only if M is nearly tangent
(within O(h^2)) to the coordinate flat {x_fixed = c} spanned by that face.  The
lattice of such flats has spacing h in R^j and the near-critical values form a
hypersurface in R^j, so the expected number of ambiguous faces with j fixed
directions behaves like h^(2-j):
    j = 1 :  O(h)  -> 0
    j = 2 :  O(1)
    j >= 3:  grows like h^(2-j)
Ambiguous (k+1)-faces have j = n-k-1 = m-1 fixed directions.

For k=1 a 2-face is ambiguous iff more than two of its edges carry a transversal
segment of M_T (the PL zero set in a triangulated square is then two arcs).
We count ambiguous 2-faces over random rotations of M and random grid offsets.
"""
import itertools, sys
import numpy as np

def crossings(F, n, lo, hi, N, shift):
    xs = [np.linspace(lo, hi, N + 1) + shift[i] for i in range(n)]
    g = np.meshgrid(*xs, indexing="ij")
    s = F(*g) > 0
    return s

def ambiguous_2faces(s, n, N):
    count = 0
    for a, b in itertools.combinations(range(n), 2):
        def sl(da, db):
            return tuple(slice(da, N + da) if i == a else slice(db, N + db) if i == b
                         else slice(0, N + 1) for i in range(n))
        s00, s10, s01, s11 = s[sl(0, 0)], s[sl(1, 0)], s[sl(0, 1)], s[sl(1, 1)]
        ncross = ((s00 != s10).astype(int) + (s01 != s11) + (s00 != s01) + (s10 != s11))
        count += int((ncross > 2).sum())
    return count

def rot(n, rng):
    q, r = np.linalg.qr(rng.standard_normal((n, n)))
    return q * np.sign(np.diag(r))

def make(kind, Q):
    if kind == "T2":      # m = 2 in R^3
        def F(x, y, z):
            X = [Q[i, 0] * x + Q[i, 1] * y + Q[i, 2] * z for i in range(3)]
            return (np.sqrt(X[0]**2 + X[1]**2) - 0.6) ** 2 + X[2]**2 - 0.25**2
        return F, 3
    if kind == "S2xS1":   # m = 3 in R^4
        def F(a, b, c, d):
            X = [Q[i, 0] * a + Q[i, 1] * b + Q[i, 2] * c + Q[i, 3] * d for i in range(4)]
            return (np.sqrt(X[0]**2 + X[1]**2) - 0.6) ** 2 + X[2]**2 + X[3]**2 - 0.25**2
        return F, 4
    if kind == "S3xS1":   # m = 4 in R^5
        def F(*x):
            X = [sum(Q[i, j] * x[j] for j in range(5)) for i in range(5)]
            return (np.sqrt(X[0]**2 + X[1]**2) - 0.6) ** 2 + X[2]**2 + X[3]**2 + X[4]**2 - 0.25**2
        return F, 5

if __name__ == "__main__":
    kind = sys.argv[1]
    Ns = [int(v) for v in sys.argv[2].split(",")]
    trials = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    rng = np.random.default_rng(1)
    rots = [rot({"T2": 3, "S2xS1": 4, "S3xS1": 5}[kind], rng) for _ in range(trials)]
    shifts = [rng.uniform(0, 1, 5) for _ in range(trials)]
    for N in Ns:
        h = 2.0 / N
        vals = []
        for Q, sh in zip(rots, shifts):
            F, n = make(kind, Q)
            s = crossings(F, n, -1.0, 1.0, N, sh[:n] * h)
            vals.append(ambiguous_2faces(s, n, N))
        frac = np.mean([v > 0 for v in vals])
        print(f"{kind:6s} N={N:3d} h={h:.4f} | ambiguous 2-faces per grid: mean {np.mean(vals):7.2f} "
              f"| values {vals} | P(any) = {frac:.2f}", flush=True)
