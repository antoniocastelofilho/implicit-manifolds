"""Slicing of polytopal cells (combinatorial skeletons read by polfile) by hyperplanes.
A cell of dimension d is cut by {x : a.x = c} into a cell of dimension d-1 whose
j-faces are the cuts of the (j+1)-faces of the original cell."""
import numpy as np
from polfile import Cell

def slice_cell(c, a, cval):
    s = c.coords @ a - cval
    pos = s >= 0
    d = len(c.faces) - 1
    if pos.all() or (~pos).all():
        return None
    # new vertices on crossing edges
    E = c.faces[1]; newv = {}; pts = []
    for ei, (p, q) in enumerate(E):
        if pos[p] != pos[q]:
            t = s[p] / (s[p] - s[q]); newv[ei] = len(pts)
            pts.append(c.coords[p] + t * (c.coords[q] - c.coords[p]))
    if not pts:
        return None
    faces = [None]
    # level j of the old cell (j>=2) gives level j-1 of the new one
    prev_map = {ei: [vi] for ei, vi in newv.items()}      # old (j-1)-face -> new (j-2)-faces
    for j in range(2, d + 1):
        newfaces = []; cur_map = {}
        for fi, f in enumerate(c.faces[j]):
            sub = [x for b in f for x in prev_map.get(b, [])]
            if not sub:
                continue
            if j == 2:
                if len(sub) == 2:
                    groups = [sub]
                else:            # non-convex polygon: pair the crossings in cyclic order
                    sub = order_polygon_crossings(c, f, newv)
                    groups = [sub[i:i + 2] for i in range(0, len(sub) - 1, 2)]
                ids = []
                for g in groups:
                    ids.append(len(newfaces)); newfaces.append(list(g))
                cur_map[fi] = ids
            else:
                # the section of a non-convex face may be disconnected: split it into components,
                # two (j-2)-faces being connected when they share a (j-3)-face
                sub = sorted(set(sub)); lower = faces[-1]
                parent = {x: x for x in sub}
                def find(x):
                    while parent[x] != x:
                        parent[x] = parent[parent[x]]; x = parent[x]
                    return x
                owner = {}
                for x in sub:
                    for b in lower[x]:
                        if b in owner: parent[find(x)] = find(owner[b])
                        else: owner[b] = x
                comps = {}
                for x in sub: comps.setdefault(find(x), []).append(x)
                ids = []
                for comp in comps.values():
                    ids.append(len(newfaces)); newfaces.append(sorted(comp))
                cur_map[fi] = ids
        faces.append(newfaces); prev_map = cur_map
    n = Cell(); n.g = c.g; n.grid = c.grid; n.coords = np.array(pts); n.labels = [c.labels[0]] * len(pts)
    n.faces = faces
    return n

def order_polygon_crossings(c, f, newv):
    E = c.faces[1]; es = [tuple(E[e]) for e in f]
    adj = {}
    for k, (p, q) in enumerate(es):
        adj.setdefault(p, []).append((q, f[k])); adj.setdefault(q, []).append((p, f[k]))
    start = es[0][0]; prev = None; cur = start; out = []; seen = set()
    while True:
        nxt = [(v, e) for v, e in adj[cur] if e not in seen]
        if not nxt: break
        v, e = nxt[0]; seen.add(e)
        if e in newv: out.append(newv[e])
        cur = v
        if cur == start: break
    return out

def slice_all(cells, a, cval):
    a = np.asarray(a, float); out = []
    for c in cells:
        r = slice_cell(c, a, cval)
        if r is not None and len(r.faces) > 1: out.append(r)
    return out

def axis(n, i):
    a = np.zeros(n); a[i] = 1.0; return a

def iter_pol(path):
    """streaming version of polfile.read_pol: yields cells one by one"""
    f = open(path)
    def toks():
        for ln in f:
            t = ln.split()
            if t: yield t
    it = toks()
    n, k = map(int, next(it)[:2]); next(it); d = n - k
    for head in it:
        if head[0] == '-1': break
        h = [int(v) for v in head]; nc = int(next(it)[0])
        for _ in range(nc):
            nv = int(next(it)[0]); co = []; lab = []
            for _ in range(nv):
                r = next(it); co.append([float(v) for v in r[-n:]]); lab.append([int(v) for v in r[:-n]])
            faces = [None]
            for j in range(1, d + 1):
                m = int(next(it)[0]); faces.append([[int(v) - 1 for v in next(it)] for _ in range(m)])
            c = Cell(); c.g = h[0]; c.grid = h[1:]; c.coords = np.array(co); c.labels = lab; c.faces = faces
            yield c
