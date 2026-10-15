"""Topological data analysis of the samples (persistent homology of Rips complexes, gudhi) and Betti
numbers of the reconstructed cell complexes (fan triangulation of the polygons, vertex identification)."""
import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, pickle; sys.path.insert(0, REPO + '/python/implicit_manifolds')
import numpy as np, gudhi
from polfile import read_pol, polygons
R = REPO + '/book/runs/ml/'


def betti_of_cells(cells, m, nd=6):
    key = {}; st = gudhi.SimplexTree()
    def vid(p):
        k = tuple(np.round(p, nd))
        if k not in key: key[k] = len(key)
        return key[k]
    for c in cells:
        ids = [vid(p) for p in c.coords]
        if m == 1:
            for e in c.faces[1]: st.insert([ids[e[0]], ids[e[1]]])
        else:
            for cyc in polygons(c):
                g = len(key); key[('c', id(c), tuple(cyc))] = g     # centre vertex of the fan
                for i in range(len(cyc)):
                    st.insert([g, ids[cyc[i]], ids[cyc[(i + 1) % len(cyc)]]])
    st.compute_persistence(persistence_dim_max=True)
    b = st.betti_numbers()
    return (b + [0] * 3)[:m + 1]


def rips_persistence(X, npts, max_edge, dim, sparse=None, seed=0):
    rng = np.random.default_rng(seed)
    # farthest point subsample, which covers the sample evenly
    idx = [int(rng.integers(len(X)))]; d = np.linalg.norm(X - X[idx[0]], axis=1)
    for _ in range(npts - 1):
        i = int(np.argmax(d)); idx.append(i); d = np.minimum(d, np.linalg.norm(X - X[i], axis=1))
    Y = X[idx]
    rc = gudhi.RipsComplex(points=Y, max_edge_length=max_edge, sparse=sparse)
    st = rc.create_simplex_tree(max_dimension=dim + 1)
    dg = st.persistence()
    return Y, dg, st.num_simplices()


if __name__ == '__main__':
    out = {}
    for tag, m, npts, me, sp in [('curve', 1, 150, 1.2, None), ('surface', 2, 350, 0.9, 0.3), ('torus4', 2, 350, 1.4, 0.3),
                                  ('torus6', 2, 400, 1.6, 0.3), ('swiss', 2, 300, 0.8, 0.3)]:
        D = pickle.load(open(R + tag + '.pkl', 'rb'))
        n_, k_, Dv, cells = read_pol(D['cells_path'])
        if D['trim'] is not None:
            from scipy.spatial import cKDTree
            dd, _ = cKDTree(D['X']).query(np.array([c.coords.mean(0) for c in cells])); cells = [c for c, x in zip(cells, dd) if x < D['trim']]
        b = betti_of_cells(cells, m)
        Y, dg, ns = rips_persistence(D['X'], npts, me, m, sp)
        out[tag] = dict(betti=b, dgm=dg, nsimp=ns, npts=npts, max_edge=me)
        long = {}
        for dim, (bi, de) in dg:
            if de - bi > 0.25 * me or de == float('inf'): long[dim] = long.get(dim, 0) + 1
        print(tag, 'betti of reconstruction', b, 'long bars', long, 'simplices', ns, flush=True)
    pickle.dump(out, open(R + 'tda.pkl', 'wb'))
