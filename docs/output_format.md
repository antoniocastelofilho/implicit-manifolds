# Output format (`.pol`)

All programs write a text file with the following structure (Appendix A of the book). Indices of vertices and cells
are 1-based within each polytope; grid labels are 0-based.

```
n  k
N_1 N_2 ... N_n

<block of grid cell>
<block of grid cell>
...
-1
```

Each block describes the polytopes of one grid cell (one simplex for the simplicial methods, one cube for the
hypercube methods):

```
g G_1 ... G_n       cell label and coordinates
p                   number of polytopes
 (for each polytope)
 V                  number of vertices
 l_0..l_k x_1..x_n  labels of the k-face, then the coordinates
 E_1                number of edges
 a b                vertices of each edge
 E_2                number of 2-cells
 i_1 i_2 ...        edges of each 2-cell
 ...
 E_m                number of m-cells
 i_1 i_2 ...        (m-1)-cells of each m-cell
```

For the simplicial methods `p = 1` and each vertex line starts with the `k+1` labels of the vertices of the
`k`-face that contains it. For the hypercube methods `p` may exceed 1, and the vertex lines start with the labels
that identify the vertex in the local frame of the cube (2 labels in the MATLAB CMH code, $2^k$ in the generalized
methods). Readers should rely on the **last `n` entries** of each vertex line, the coordinates.

A `j`-cell is given by its `(j-1)`-faces, so the vertices of a cell are obtained by descending the lists; the cyclic
order of the vertices of a 2-cell is recovered by walking along its edges. The global complex is obtained by
identifying vertices with equal coordinates across cells. `python/implicit_manifolds/polfile.py` implements all of
this.
