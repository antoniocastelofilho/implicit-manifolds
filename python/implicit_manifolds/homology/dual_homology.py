"""
Homology of the dual cubical complex Q of the FCH/GCMH combinatorial skeleton.

For each k-face J of the grid (lower corner p, free directions A) containing a
transversal Freudenthal--Kuhn (K1) k-simplex, the dual cell is the elementary
(n-k)-cube
    D(J) = prod_{i in A} {p_i}  x  prod_{i not in A} [p_i - 1, p_i]
in the lattice of hypercube indices.  Q = union of the D(J).  By the dual theorem
of the companion notes, Q realizes the dual complex M'.

Homology is computed over GF(2) and GF(P) (P a large prime, i.e. rational Betti
numbers), after Mrozek--Batko coreduction.  Different Betti numbers over GF(2) and
over Q signal 2-torsion.
"""
import itertools, sys, time
from collections import deque
import numpy as np

P = 1_000_003

# ---------------------------------------------------------------------------
# transversality (property P) on K1 simplices of k-faces
# ---------------------------------------------------------------------------
def k1_simplices_of_face(A):
    out = []
    for perm in itertools.permutations(A):
        path = [frozenset()]
        for a in perm:
            path.append(path[-1] | {a})
        out.append(path)
    return out

def property_P(V):
    m, kp1, k = V.shape
    A = np.ones((m, k + 1, k + 1))
    A[:, 1:, :] = np.transpose(V, (0, 2, 1))
    det = np.linalg.det(A)
    ok = np.abs(det) > 1e-13
    res = np.zeros(m, dtype=bool)
    if ok.any():
        b = np.zeros((int(ok.sum()), k + 1, 1)); b[:, 0, 0] = 1
        lam = np.linalg.solve(A[ok], b)[..., 0]
        res[ok] = np.all(lam > 0, axis=1)
    return res

def dual_top_cubes(F, n, k, lo, hi, N):
    xs = [np.linspace(lo, hi, N + 1) + 1e-3 * np.sqrt(2 + i) for i in range(n)]
    grids = np.meshgrid(*xs, indexing="ij")
    Fv = F(*grids)
    Fv = np.stack(Fv, axis=-1) if isinstance(Fv, (list, tuple)) else Fv[..., None]
    tops, multi = [], 0
    for A in itertools.combinations(range(n), k):
        base = tuple(slice(0, N) if i in A else slice(0, N + 1) for i in range(n))
        shape = Fv[base].shape[:-1]
        cnt = np.zeros(shape, dtype=np.int8)
        for s in k1_simplices_of_face(A):
            vals = []
            for vset in s:
                sl = tuple(slice(1, N + 1) if i in vset else base[i] for i in range(n))
                vals.append(Fv[sl])
            V = np.stack(vals, axis=-2).reshape(-1, k + 1, k)
            cnt += property_P(V).reshape(shape).astype(np.int8)
        idx = np.argwhere(cnt > 0)
        multi += int((cnt > 1).sum())
        fixed = [i for i in range(n) if i not in A]
        mask = sum(1 << i for i in fixed)
        for p in idx:
            c = p.copy()
            c[fixed] -= 1
            if (c < 0).any() or (c > N - 1).any():
                raise RuntimeError("manifold touches the domain boundary")
            tops.append((tuple(int(x) for x in c), mask))
    return tops, multi

# ---------------------------------------------------------------------------
# cubical complex, coreduction, ranks
# ---------------------------------------------------------------------------
def faces_of(cell, n):
    corner, mask = cell
    dirs = [i for i in range(n) if mask >> i & 1]
    out = []
    for j, i in enumerate(dirs):
        sgn = 1 if j % 2 == 0 else -1
        m2 = mask & ~(1 << i)
        c_plus = list(corner); c_plus[i] += 1
        out.append(((tuple(c_plus), m2), sgn))
        out.append(((corner, m2), -sgn))
    return out

def build_complex(tops, n):
    cells = {}
    order = []
    stack = list(set(tops))
    for c in stack:
        cells[c] = None
    while stack:
        c = stack.pop()
        order.append(c)
        for f, _ in faces_of(c, n):
            if f not in cells:
                cells[f] = None
                stack.append(f)
    ids = {c: i for i, c in enumerate(order)}
    dim = [bin(c[1]).count("1") for c in order]
    bd = [dict() for _ in order]
    cbd = [dict() for _ in order]
    for c, i in ids.items():
        for f, s in faces_of(c, n):
            j = ids[f]
            bd[i][j] = s
            cbd[j][i] = s
    return dim, bd, cbd

def coreduce(dim, bd, cbd):
    alive = [True] * len(dim)
    v0 = dim.index(0)
    alive[v0] = False
    for c in list(cbd[v0]):
        del bd[c][v0]
    q = deque(cbd[v0].keys())
    removed = 1
    while q:
        a = q.popleft()
        if not alive[a] or len(bd[a]) != 1:
            continue
        (b, _), = bd[a].items()
        for x in (a, b):
            alive[x] = False
            for f in bd[x]:
                cbd[f].pop(x, None)
            for c in cbd[x]:
                bd[c].pop(x, None)
                q.append(c)
        removed += 2
    return alive

def rank_mod(rows, p):
    """rows: list of dict col->coeff. Gaussian elimination mod p."""
    pivots = {}
    r = 0
    for row in rows:
        v = {c: x % p for c, x in row.items() if x % p}
        while v:
            c = max(v)
            if c in pivots:
                pr = pivots[c]
                f = v[c] * pow(pr[c], p - 2, p) % p
                for cc, xx in pr.items():
                    nv = (v.get(cc, 0) - f * xx) % p
                    if nv: v[cc] = nv
                    else: v.pop(cc, None)
            else:
                pivots[c] = v; r += 1; break
    return r

def betti(dim, bd, alive, p):
    D = max(dim)
    cnt = [0] * (D + 1)
    rows = {d: [] for d in range(D + 1)}
    for i, a in enumerate(alive):
        if a:
            cnt[dim[i]] += 1
            rows[dim[i]].append(bd[i])
    rk = {d: rank_mod(rows[d], p) if d > 0 else 0 for d in range(D + 1)}
    b = [cnt[d] - rk[d] - rk.get(d + 1, 0) for d in range(D + 1)]
    b[0] += 1  # the vertex removed at the start of the coreduction
    return b

# ---------------------------------------------------------------------------
# primal complex (Khalimsky model): cells = grid faces L, dim L >= k, containing
# a 0-cell; rank = dim L - k; boundary = facets of L that are cells.
# ---------------------------------------------------------------------------
def primal_cells(tops, n):
    """tops are dual cubes (corner c, mask of fixed dirs); recover the 0-cells J
    (corner p, free dirs A) and generate all faces containing them."""
    full = (1 << n) - 1
    cells = set()
    for c, fixedmask in tops:
        p = list(c)
        fixed = [i for i in range(n) if fixedmask >> i & 1]
        for i in fixed:
            p[i] += 1                      # undo the shift in dual_top_cubes
        freemask = full & ~fixedmask
        # each fixed dir: keep (0), extend left (1), extend right (2)
        for choice in itertools.product((0, 1, 2), repeat=len(fixed)):
            q = list(p); m = freemask
            for i, ch in zip(fixed, choice):
                if ch == 1: q[i] -= 1; m |= 1 << i
                elif ch == 2: m |= 1 << i
            cells.add((tuple(q), m))
    return cells

def primal_complex(cells, n, k):
    order = list(cells)
    ids = {c: i for i, c in enumerate(order)}
    dim = [bin(c[1]).count("1") - k for c in order]
    bd = [dict() for _ in order]; cbd = [dict() for _ in order]
    for c, i in ids.items():
        corner, m = c
        for i_dir in range(n):
            if m >> i_dir & 1:
                for sh in (0, 1):
                    q = list(corner); q[i_dir] += sh
                    f = (tuple(q), m & ~(1 << i_dir))
                    j = ids.get(f)
                    if j is not None:
                        bd[i][j] = 1; cbd[j][i] = 1
    return order, dim, bd, cbd

def closure_betti_z2(top, cells, n, k):
    """Z2 Betti numbers of cl(top) inside the primal complex."""
    corner, m = top
    dirs = [i for i in range(n) if m >> i & 1]
    sub = []
    for choice in itertools.product((0, 1, 2), repeat=len(dirs)):
        q = list(corner); mm = m
        for i, ch in zip(dirs, choice):
            if ch == 0: q[i] += 0; mm &= ~(1 << i)
            elif ch == 1: q[i] += 1; mm &= ~(1 << i)
        f = (tuple(q), mm)
        if f in cells: sub.append(f)
    order, dim, bd, cbd = primal_complex(set(sub), n, k)
    D = max(dim)
    cnt = [dim.count(d) for d in range(D + 1)]
    rk = {d: rank_mod([bd[i] for i in range(len(order)) if dim[i] == d], 2) for d in range(1, D + 1)}
    return [cnt[d] - rk.get(d, 0) - rk.get(d + 1, 0) for d in range(D + 1)]

def run(F, n, k, lo, hi, N, label, expect):
    t = time.time()
    tops, multi = dual_top_cubes(F, n, k, lo, hi, N)
    # dual
    dim, bd, cbd = build_complex(tops, n)
    ncell = len(dim)
    chi = sum((-1) ** d for d in dim)
    alive = coreduce(dim, bd, cbd)
    b2 = betti(dim, bd, alive, 2)
    bq = betti(dim, bd, alive, P)
    # primal (independent code path, Z2)
    pcells = primal_cells(tops, n)
    _, pdim, pbd, pcbd = primal_complex(pcells, n, k)
    palive = coreduce(pdim, pbd, pcbd)
    pb2 = betti(pdim, pbd, palive, 2)
    # regularity of primal top cells
    toplevel = [c for c in pcells if c[1] == (1 << n) - 1]
    bad = 0
    for c in toplevel:
        bb = closure_betti_z2(c, pcells, n, k)
        if bb[0] != 1 or any(bb[1:]):
            bad += 1
    pal = "yes" if b2 == b2[::-1] else "NO"
    print(f"{label:22s} n={n} k={k} N={N:2d} | 0-cells {len(tops):6d} | dual cells {ncell:7d} "
          f"| b(X;Z2) {b2} b(X;Q) {bq} | primal b(M;Z2) {pb2} | palindromic {pal} "
          f"| non-acyclic top cells {bad}/{len(toplevel)} | expected {expect} "
          f"| {time.time()-t:5.1f}s", flush=True)

if __name__ == "__main__":
    which = sys.argv[1:] or ["all"]
    R, r = 0.6, 0.3
    torus = lambda x, y, z: (np.sqrt(x**2 + y**2) - 0.6) ** 2 + z**2 - 0.25**2
    thin = lambda x, y, z: (np.sqrt(x**2 + y**2) - 0.6) ** 2 + z**2 - 0.07**2
    cases = {
        "T2": (torus, 3, 1, -1, 1, 16, "(1,2,1)"),
        "thinT2_N8": (thin, 3, 1, -1, 1, 8, "(1,2,1)"),
        "thinT2_N16": (thin, 3, 1, -1, 1, 16, "(1,2,1)"),
        "thinT2_N40": (thin, 3, 1, -1, 1, 40, "(1,2,1)"),
        "Clifford": (lambda a, b, c, d: (a*a + b*b - 0.5, c*c + d*d - 0.5), 4, 2, -1, 1, 12, "(1,2,1)"),
        "S2xS1": (lambda x, y, z, w: (np.sqrt(x**2 + y**2) - R) ** 2 + z**2 + w**2 - r**2,
                  4, 1, -1, 1, 16, "(1,1,1,1)"),
        "S3xS1": (lambda a, b, c, d, e: (np.sqrt(a**2 + b**2 + c**2 + d**2) - R) ** 2 + e**2 - r**2,
                  5, 1, -1, 1, 10, "(1,1,0,1,1)"),
        "T3": (lambda a, b, c, d, e, f: (a*a + b*b - 1, c*c + d*d - 1, e*e + f*f - 1),
               6, 3, -1.25, 1.25, 7, "(1,3,3,1)"),
        "S2xS2": (lambda a, b, c, d, e, f: (a*a + b*b + c*c - 0.7, d*d + e*e + f*f - 0.7),
                  6, 2, -1.1, 1.1, 7, "(1,0,2,0,1)"),
    }
    for name, (Ff, n, k, lo, hi, N, ex) in cases.items():
        if "all" in which or name in which:
            run(Ff, n, k, lo, hi, N, name, ex)
