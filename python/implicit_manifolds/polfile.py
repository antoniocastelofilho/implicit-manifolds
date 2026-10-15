"""Reader for the combinatorial-skeleton (.pol) files written by the
Marching/Continuation Simplex and Hypercube codes, plus helpers to build a
global cell complex (vertices identified by coordinates) and compute
Euler characteristics."""
import numpy as np


class Cell:
    """One polytope (connected component) of the approximation inside a
    simplex or hypercube.  levels[0] = vertex coords (m x n);
    labels = mesh-vertex/face labels printed with each vertex;
    faces[j] (j>=1) = list of lists with indices (1-based in file, stored
    0-based here) of (j-1)-faces bounding each j-face."""
    __slots__ = ("g", "grid", "coords", "labels", "faces")


def read_pol(path):
    tok_lines = [ln.split() for ln in open(path)]
    lines = [t for t in tok_lines if t]
    n, k = int(lines[0][0]), int(lines[0][1])
    division = [int(v) for v in lines[1]]
    i = 2
    cells = []
    d = n - k
    while i < len(lines):
        if lines[i][0] == "-1":
            break
        head = [int(v) for v in lines[i]]
        g, grid = head[0], head[1:]
        i += 1
        nc = int(lines[i][0]); i += 1
        for _ in range(nc):
            nv = int(lines[i][0]); i += 1
            coords, labels = [], []
            for _ in range(nv):
                row = lines[i]; i += 1
                coords.append([float(v) for v in row[-n:]])
                labels.append([int(v) for v in row[:-n]])
            faces = [None]
            for j in range(1, d + 1):
                m = int(lines[i][0]); i += 1
                fj = []
                for _ in range(m):
                    fj.append([int(v) - 1 for v in lines[i]]); i += 1
                faces.append(fj)
            c = Cell()
            c.g, c.grid = g, grid
            c.coords = np.array(coords)
            c.labels = labels
            c.faces = faces
            cells.append(c)
    return n, k, division, cells


def cell_vertex_sets(c, j):
    """vertex index sets (local) of every j-face of cell c"""
    if j == 0:
        return [frozenset([i]) for i in range(len(c.coords))]
    prev = cell_vertex_sets(c, j - 1)
    return [frozenset().union(*[prev[b] for b in f]) for f in c.faces[j]]


def polygons(c):
    """ordered vertex cycles (local indices) of the 2-faces of c"""
    if len(c.faces) < 3:
        return []
    edges = [tuple(e) for e in c.faces[1]]
    out = []
    for f in c.faces[2]:
        es = [edges[e] for e in f]
        adj = {}
        for a, b in es:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        start = es[0][0]
        cyc, prev, cur = [start], None, start
        while True:
            nxt = [v for v in adj[cur] if v != prev]
            if not nxt:
                break
            nx = nxt[0]
            if nx == start:
                break
            cyc.append(nx)
            prev, cur = cur, nx
            if len(cyc) > len(adj):
                break
        out.append(cyc)
    return out


def global_complex(cells, ndig=6):
    """identify vertices by rounded coordinates and return, for each
    dimension j, the set of distinct j-faces (as frozensets of global ids)"""
    key2id, X = {}, []
    allfaces = {}
    for c in cells:
        gid = []
        for p in c.coords:
            key = tuple(np.round(p, ndig))
            if key not in key2id:
                key2id[key] = len(X)
                X.append(p)
            gid.append(key2id[key])
        d = len(c.faces) - 1
        for j in range(0, d + 1):
            for s in cell_vertex_sets(c, j):
                allfaces.setdefault(j, set()).add(frozenset(gid[v] for v in s))
    return np.array(X), allfaces


def euler(cells):
    X, F = global_complex(cells)
    counts = [len(F[j]) for j in sorted(F)]
    chi = sum((-1) ** j * c for j, c in enumerate(counts))
    return counts, chi
