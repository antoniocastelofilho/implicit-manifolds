"""
Test of the hypothesis: "if a hypercube has no ambiguity, its polytope approximates
M cap H topologically".

k = 1.  Gamma is built exactly as FCH does: vertices on crossed grid edges, arcs in the
2-faces (K1 diagonal from the lower to the upper corner; with four crossings the arcs
cut off the two corners off the diagonal).

For every active n-cube H we record:
  * fully unambiguous: Gamma_L connected for every face L of H with dim L >= 2;
  * top cell regular: the boundary of the cell H (in the Khalimsky model, which is the
    intrinsic skeleton inside a fully unambiguous cube) has the Z2-homology of S^{m-1}.
A fully unambiguous cube with an irregular top cell is a counterexample to the
hypothesis (connected boundary graph, but M cap H is not a ball).
"""
import itertools, sys, time
import numpy as np
from dual_homology import dual_top_cubes, primal_cells
from convergence import boundary_betti_z2

class UF:
    def __init__(self): self.p = {}
    def f(self, a):
        self.p.setdefault(a, a)
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]; a = self.p[a]
        return a
    def u(self, a, b): self.p[self.f(a)] = self.f(b)

def analyse(F, n, lo, hi, N, label, seed_shift=None):
    t0 = time.time()
    shift = seed_shift if seed_shift is not None else np.array([1e-3 * np.sqrt(2 + i) for i in range(n)])
    xs = [np.linspace(lo, hi, N + 1) + shift[i] for i in range(n)]
    S = F(*np.meshgrid(*xs, indexing="ij")) > 0
    def cross(u, i):
        v = list(u); v[i] += 1
        return S[tuple(u)] != S[tuple(v)]
    tops, _ = dual_top_cubes(F, n, 1, lo, hi, N)      # same shift convention (1e-3*sqrt(2+i))
    cells = primal_cells(tops, n)
    full = (1 << n) - 1
    cubes = [c for c in cells if c[1] == full]
    n_unamb = n_unamb_bad = n_amb = n_amb_bad = 0
    examples = []
    for corner, _ in cubes:
        # arcs inside each 2-face of the cube, union-find per face
        unamb = True
        for dimL in range(2, n + 1):
            for A in itertools.combinations(range(n), dimL):
                fixed = [i for i in range(n) if i not in A]
                for offs in itertools.product((0, 1), repeat=len(fixed)):
                    p = list(corner)
                    for i, o in zip(fixed, offs): p[i] += o
                    uf = UF(); nodes = []
                    for T in itertools.product((0, 1), repeat=dimL):
                        u = list(p)
                        for i, t in zip(A, T): u[i] += t
                        for i, t in zip(A, T):
                            if t == 0 and cross(u, i):
                                nodes.append((tuple(u), i)); uf.f((tuple(u), i))
                    for a, b in itertools.combinations(A, 2):
                        others = [i for i in A if i not in (a, b)]
                        for T in itertools.product((0, 1), repeat=len(others)):
                            u = list(p)
                            for i, t in zip(others, T): u[i] += t
                            ua = list(u); ua[a] += 1
                            ub = list(u); ub[b] += 1
                            e1, e2 = (tuple(u), a), (tuple(u), b)
                            e3, e4 = (tuple(ua), b), (tuple(ub), a)
                            cr = [e for e in (e1, e2, e3, e4) if cross(list(e[0]), e[1])]
                            if len(cr) == 2:
                                uf.u(cr[0], cr[1])
                            elif len(cr) == 4:
                                uf.u(e1, e3); uf.u(e2, e4)
                    if len({uf.f(x) for x in nodes}) > 1:
                        unamb = False
        regular = boundary_betti_z2((corner, full), cells, n, 1)
        if unamb:
            n_unamb += 1
            if not regular:
                n_unamb_bad += 1
                if len(examples) < 3: examples.append(corner)
        else:
            n_amb += 1
            n_amb_bad += (not regular)
    print(f"{label:16s} N={N:3d} | cubes {len(cubes):6d} | fully unambiguous {n_unamb:6d} "
          f"(irregular top cell: {n_unamb_bad}) | ambiguous {n_amb:4d} (irregular: {n_amb_bad}) "
          f"| {time.time()-t0:5.1f}s", flush=True)
    return n_unamb_bad, examples

if __name__ == "__main__":
    which = sys.argv[1]
    rng = np.random.default_rng(3)
    def rot(n):
        q, r = np.linalg.qr(rng.standard_normal((n, n))); return q * np.sign(np.diag(r))
    if which == "T2":
        Q = rot(3)
        def F(*x):
            X = [sum(Q[i, j] * x[j] for j in range(3)) for i in range(3)]
            return (np.sqrt(X[0]**2 + X[1]**2) - 0.6) ** 2 + X[2]**2 - 0.25**2
        for N in (8, 16, 32): analyse(F, 3, -1, 1, N, "rot T2")
    if which == "S2xS1":
        Q = rot(4)
        def F(*x):
            X = [sum(Q[i, j] * x[j] for j in range(4)) for i in range(4)]
            return (np.sqrt(X[0]**2 + X[1]**2) - 0.6) ** 2 + X[2]**2 + X[3]**2 - 0.25**2
        for N in (8, 12, 16, 20): analyse(F, 4, -1, 1, N, "rot S2xS1")
    if which == "S3xS1":
        Q = rot(5)
        def F(*x):
            X = [sum(Q[i, j] * x[j] for j in range(5)) for i in range(5)]
            return (np.sqrt(X[0]**2 + X[1]**2) - 0.6) ** 2 + X[2]**2 + X[3]**2 + X[4]**2 - 0.3**2
        for N in (8, 10, 12): analyse(F, 5, -1, 1, N, "rot S3xS1")
