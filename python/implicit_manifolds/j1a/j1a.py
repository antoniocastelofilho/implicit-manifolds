"""J1a triangulation: construction (Castelo et al.) and conformity check."""
import itertools, math, numpy as np
from collections import Counter
from fractions import Fraction as Fr

class Grid:
    def __init__(s, m, N, H=Fr(1)):
        s.m, s.N, s.H = m, N, H        # N level-0 blocks per direction, half-size H at level 0
        s.leaves = {(0, g) for g in itertools.product(range(N), repeat=m)}
    def box(s, b):
        r, g = b; hs = s.H / 2**r
        c = tuple(2*hs*gi + hs for gi in g); return c, hs
    def contains(s, b, x):
        c, hs = s.box(b); return all(abs(xi-ci) <= hs for xi, ci in zip(x, c))
    def leaves_at(s, x): return [b for b in s.leaves if s.contains(b, x)]
    def split(s, b):
        r, g = b; s.leaves.remove(b)
        for o in itertools.product((0,1), repeat=s.m):
            s.leaves.add((r+1, tuple(2*gi+oi for gi, oi in zip(g, o))))
    def neighbours(s, b):
        c, hs = s.box(b); out=set()
        for d in itertools.product((-1,0,1), repeat=s.m):
            if any(d):
                x = tuple(ci + di*(hs + hs/8) for ci, di in zip(c, d))
                out.update(s.leaves_at(x))
        return out
    def refine(s, b):
        """split b, enforcing the balance restriction (adjacent levels differ by <= 1)"""
        r = b[0]
        for nb in list(s.neighbours(b)):
            if nb in s.leaves and nb[0] < r: s.refine(nb)
        if b in s.leaves: s.split(b)

def refined_faces(G, b):
    m = G.m; r = b[0]; c, hs = G.box(b); out = []
    for k in range(1, m):                       # number of fixed coordinates
        for axes in itertools.combinations(range(m), k):
            for sides in itertools.product((-1,1), repeat=k):
                lab = frozenset(zip(axes, sides))
                x = list(c)
                for a, sd in lab: x[a] = c[a] + sd*hs
                free = [i for i in range(m) if i not in axes]
                for i in free: x[i] = c[i] + hs/2 - hs/16     # inside one quadrant of the face
                if any(L[0] == r+1 for L in G.leaves_at(tuple(x))):
                    out.append(lab)
    return out

def simplices(G, b):
    m = G.m; c, hs = G.box(b); Phi = refined_faces(G, b); S = []
    for s in itertools.product((-1,1), repeat=m):
        for pi in itertools.permutations(range(m)):
            steps = [(pi[j], s[pi[j]]) for j in range(m)]
            h = None
            for hh in range(1, m+1):
                if any(phi <= frozenset(steps[:hh]) for phi in Phi): h = hh; break
            if h is None or h == m:
                v = [0]*m; P = [tuple(v)]
                for j in range(m): v[pi[j]] += s[pi[j]]; P.append(tuple(v))
                S.append(P); continue
            for t in itertools.product((-1,1), repeat=m-h):
                v = [Fr(0)]*m; P = [tuple(v)]
                for j in range(h-1): v[pi[j]] += s[pi[j]]; P.append(tuple(v))
                v[pi[h-1]] += s[pi[h-1]]
                for j in range(h, m): v[pi[j]] += Fr(s[pi[j]], 2)       # centre of the sub-face
                P.append(tuple(v))
                for j in range(h, m): v[pi[j]] += Fr(t[j-h], 2); P.append(tuple(v))
                S.append(P)
    return [[tuple(ci + hs*vi for ci, vi in zip(c, p)) for p in P] for P in S]

def check(G):
    m = G.m; allS = []
    for b in G.leaves: allS += simplices(G, b)
    vol = sum(abs(Fr(int(round(np.linalg.det(np.array([[float(x) for x in np.subtract(p, P[0])] for p in P[1:]]))*10**9)), 10**9)) for P in allS) / math.factorial(m)
    fc = Counter()
    for P in allS:
        for i in range(m+1): fc[frozenset(P[:i]+P[i+1:])] += 1
    bad = [f for f, n in fc.items() if n > 2 or (n == 1 and not on_boundary(G, f))]
    degenerate = sum(1 for P in allS if len(set(P)) < m+1)
    return len(allS), float(vol), (2*G.N*G.H)**m, len(bad), degenerate

def on_boundary(G, f):
    L = 2*G.N*G.H
    for a in range(G.m):
        if all(p[a] == 0 for p in f) or all(p[a] == L for p in f): return True
    return False

if __name__ == '__main__':
    for m, N, pts in [(2, 4, [(Fr(13,10), Fr(11,10))]), (3, 3, [(Fr(13,10), Fr(11,10), Fr(9,10))])]:
        G = Grid(m, N)
        for lev in range(2):
            for x in pts:
                for b in [b for b in G.leaves_at(x)]:
                    if b in G.leaves: G.refine(b)
        print(m, Counter(b[0] for b in G.leaves), check(G))
