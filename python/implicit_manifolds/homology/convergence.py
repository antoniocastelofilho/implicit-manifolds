"""
Refinement study for the convergence conjecture.

For a sequence of grids, compute
  * Betti numbers of the dual cubical set X (Z2 and Q),
  * Betti numbers of the primal skeleton S (Z2), from the face incidences,
  * the number of cells of S that are NOT homology cells, i.e. whose boundary
    cl(x) minus {x} does not have the Z2 homology of the sphere S^{rk x - 1}.
Uses the Khalimsky model (cells = grid faces of dim >= k meeting M_T).
"""
import itertools, sys, time
import numpy as np
from dual_homology import (dual_top_cubes, build_complex, coreduce, betti, rank_mod,
                           primal_cells, primal_complex, P)

def boundary_betti_z2(cell, cells, n, k):
    corner, m = cell
    dirs = [i for i in range(n) if m >> i & 1]
    sub = []
    for choice in itertools.product((0, 1, 2), repeat=len(dirs)):
        if all(ch == 2 for ch in choice):
            continue                                  # drop the cell itself
        q = list(corner); mm = m
        for i, ch in zip(dirs, choice):
            if ch == 1: q[i] += 1
            if ch != 2: mm &= ~(1 << i)
        f = (tuple(q), mm)
        if f in cells: sub.append(f)
    r = bin(m).count("1") - k
    if not sub:
        return r == 0
    order, dim, bd, _ = primal_complex(set(sub), n, k)
    D = max(dim)
    cnt = [dim.count(d) for d in range(D + 1)]
    rk = {d: rank_mod([bd[i] for i in range(len(order)) if dim[i] == d], 2) for d in range(1, D + 1)}
    b = [cnt[d] - rk.get(d, 0) - rk.get(d + 1, 0) for d in range(D + 1)]
    sphere = [0] * (r)                                # S^{r-1}: reduced homology Z2 in degree r-1
    if r == 1:
        return b == [2]
    target = [1] + [0] * (r - 2) + [1]
    return b == target

def run(F, n, k, lo, hi, N, label):
    t = time.time()
    tops, _ = dual_top_cubes(F, n, k, lo, hi, N)
    dim, bd, cbd = build_complex(tops, n)
    alive = coreduce(dim, bd, cbd)
    b2 = betti(dim, bd, alive, 2); bq = betti(dim, bd, alive, P)
    pcells = primal_cells(tops, n)
    _, pdim, pbd, pcbd = primal_complex(pcells, n, k)
    palive = coreduce(pdim, pbd, pcbd)
    pb2 = betti(pdim, pbd, palive, 2)
    bad = sum(1 for c in pcells if not boundary_betti_z2(c, pcells, n, k))
    print(f"{label:12s} N={N:3d} h={2/N:.4f} | cells {len(pcells):7d} | non-homology cells {bad:4d} "
          f"| b(S;Z2) {pb2} | b(X;Z2) {b2} | b(X;Q) {bq} | {time.time()-t:5.1f}s", flush=True)
    return dict(N=N, cells=len(pcells), bad=bad, prim=pb2, dual2=b2, dualq=bq)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "torus"
    if which == "torus":
        thin = lambda x, y, z: (np.sqrt(x**2 + y**2) - 0.6) ** 2 + z**2 - 0.07**2
        for N in (8, 12, 16, 20, 24, 28, 32, 40, 48):
            run(thin, 3, 1, -1, 1, N, "thin T2")
    elif which == "s2s1":
        thin = lambda a, b, c, d: (np.sqrt(a*a + b*b) - 0.6) ** 2 + c*c + d*d - 0.12**2
        for N in (6, 8, 10, 12, 14, 16, 20):
            run(thin, 4, 1, -1, 1, N, "thin S2xS1")
