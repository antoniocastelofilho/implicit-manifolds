# MATLAB / Octave codes

The original research implementations, kept as they were published (they run in MATLAB and in GNU Octave ≥ 8).
Each writes the combinatorial skeleton in the `.pol` format of [../docs/output_format.md](../docs/output_format.md).

| Directory | Algorithms | Book | Examples |
|---|---|---|---|
| `marching_simplex/` | Marching Simplex, Continuation Simplex on $K_1$ | Chaps. 5–6 | `Esf.m`, `EsfCont.m` (sphere), `Complex.m` (complex curve) |
| `marching_hypercubes/` | Combinatorial Marching and Continuation Hypercubes ($k=1$), and the simplicial methods | Chap. 8 | `No.m` (knot), `Toro.m` (torus of $\mathbb R^4$) |
| `generalized_marching_hypercubes/` | Generalized Combinatorial Marching and Continuation Hypercubes (any $k$) | Chap. 9 | `Klein_Bottle.m`, `Torus_R10.m`, `S1_x_S1_x_S1.m`, `Complex_Cosine.m`, `Generalized_Hyperbola.m`, `Eight_Torus_R3.m` |

Usage pattern (Octave):

```matlab
F = @(x) x(1)^2 + x(2)^2 + x(3)^2 - 1;                          % sphere of R^3
MarchingSimplex(3, 1, [-1.1 -1.1 -1.1], [1.1 1.1 1.1], [12 12 12], F, 'sphere_ms.pol');
GeneralizedMarchingHyperCube(3, 1, [-1.1 -1.1 -1.1], [1.1 1.1 1.1], [12 12 12], 'sphere_gmhc.pol', F);
```

The runs of the table of interpreted codes in Chapter 12 of the book are in `book/common/octave_runs/` (`jobs*.m`,
called through `runjob.m`, which appends the timings to `timings.txt`).

**Erratum.** In the examples of the complex curve $z^2=w^5+1+i$, the signs of the real and imaginary parts of $w^5$
in the map $F$ are incorrect, and so is the starting point; the corrected map and the starting point $z=1$, $w=-i$
are those of Chapter 6 of the book (C++ test case `BOOK_Z2W5`).
