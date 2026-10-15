"""
Manifolds with boundary: M_T meets the boundary of the domain D.

Extended dual X_d = union of the boxes B(J,H) (Khalimsky lattice), over 0-cells J and
n-cubes H inside D containing J.  This equals the union of the truncated dual cubes
D(J) cap D, whose face poset is S' together with the duals of the cells of dS.
Expected: H(X_d) = H(M) = H(S) (primal), while the naive (transposed) dual gives the
reverse of the primal Betti numbers, i.e. H(M, dM) by Lefschetz duality.
"""
import itertools, sys, time
import numpy as np
from dual_homology import (k1_simplices_of_face, property_P, build_complex, coreduce,
                           betti, P)

def zero_cells(F, n, k, lo, hi, N):
    xs = [np.linspace(lo[i], hi[i], N + 1) + 1e-3 * np.sqrt(2 + i) for i in range(n)]
    grids = np.meshgrid(*xs, indexing="ij")
    Fv = F(*grids)
    Fv = np.stack(Fv, axis=-1) if isinstance(Fv, (list, tuple)) else Fv[..., None]
    out = []
    for A in itertools.combinations(range(n), k):
        base = tuple(slice(0, N) if i in A else slice(0, N + 1) for i in range(n))
        shape = Fv[base].shape[:-1]
        cnt = np.zeros(shape, dtype=np.int8)
        for s in k1_simplices_of_face(A):
            vals = [Fv[tuple(slice(1, N + 1) if i in vset else base[i] for i in range(n))] for vset in s]
            V = np.stack(vals, axis=-2).reshape(-1, k + 1, k)
            cnt += property_P(V).reshape(shape).astype(np.int8)
        for p in np.argwhere(cnt > 0):
            out.append((tuple(int(x) for x in p), frozenset(A)))
    return out

def extended_dual_boxes(zc, n, N):
    boxes = []
    for p, A in zc:
        fixed = [i for i in range(n) if i not in A]
        for sides in itertools.product((-1, 1), repeat=len(fixed)):
            corner = [2 * p[i] + 1 if i in A else 2 * p[i] for i in range(n)]
            mask = 0; ok = True
            for i, s in zip(fixed, sides):
                hidx = p[i] - 1 if s < 0 else p[i]          # index of the n-cube in dir i
                if hidx < 0 or hidx > N - 1:
                    ok = False; break
                if s < 0: corner[i] -= 1
                mask |= 1 << i
            if ok:
                boxes.append((tuple(corner), mask))
    return boxes

def primal(zc, n, k, N):
    full = (1 << n) - 1
    cells = set()
    for p, A in zc:
        fixed = [i for i in range(n) if i not in A]
        freemask = sum(1 << i for i in A)
        for choice in itertools.product((0, 1, 2), repeat=len(fixed)):
            q = list(p); m = freemask; ok = True
            for i, ch in zip(fixed, choice):
                if ch == 1: q[i] -= 1; m |= 1 << i
                elif ch == 2: m |= 1 << i
                if m >> i & 1 and not (0 <= q[i] <= N - 1): ok = False
            if ok: cells.add((tuple(q), m))
    order = list(cells); ids = {c: j for j, c in enumerate(order)}
    dim = [bin(c[1]).count("1") - k for c in order]
    bd = [dict() for _ in order]; cbd = [dict() for _ in order]
    for c, j in ids.items():
        corner, m = c
        for d in range(n):
            if m >> d & 1:
                for sh in (0, 1):
                    q = list(corner); q[d] += sh
                    f = ids.get((tuple(q), m & ~(1 << d)))
                    if f is not None: bd[j][f] = 1; cbd[f][j] = 1
    alive = coreduce(dim, bd, cbd)
    return betti(dim, bd, alive, 2)

def run(F, n, k, lo, hi, N, label, expect, expect_rel):
    t = time.time()
    zc = zero_cells(F, n, k, lo, hi, N)
    boxes = extended_dual_boxes(zc, n, N)
    dim, bd, cbd = build_complex(boxes, n)
    alive = coreduce(dim, bd, cbd)
    b2 = betti(dim, bd, alive, 2); bq = betti(dim, bd, alive, P)
    pb = primal(zc, n, k, N)
    print(f"{label:16s} n={n} k={k} N={N} | 0-cells {len(zc):6d} | ext. dual b(Z2) {b2} b(Q) {bq} "
          f"| primal b(Z2) {pb} | naive dual (reverse) {pb[::-1]} "
          f"| expected H(M) {expect}, H(M,dM) {expect_rel} | {time.time()-t:4.1f}s", flush=True)

if __name__ == "__main__":
    R, r = 0.6, 0.25
    run(lambda x, y, z: x*x + y*y + z*z - 0.64, 3, 1, [-1, -1, 0.1], [1, 1, 1], 16,
        "disc (hemisph.)", "(1,0,0)", "(0,0,1)")
    run(lambda x, y, z: (np.sqrt(x*x + y*y) - R) ** 2 + z*z - r*r, 3, 1, [0.05, -1, -1], [1, 1, 1], 16,
        "cylinder", "(1,1,0)", "(0,1,1)")
    run(lambda a, b, c, d: (np.sqrt(a*a + b*b) - 0.6) ** 2 + c*c + d*d - 0.09, 4, 1,
        [0.05, -1, -1, -1], [1, 1, 1, 1], 14, "S2 x [0,1]", "(1,0,1,0)", "(0,1,0,1)")
    run(lambda a, b, c, d: (a*a + b*b - 0.5, c*c + d*d - 0.5), 4, 2,
        [0.05, -1, -1, -1], [1, 1, 1, 1], 12, "annulus in R4", "(1,1,0)", "(0,1,1)")
