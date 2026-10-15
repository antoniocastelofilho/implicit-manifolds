# Python interface

## `python/fch_lib` — Fast Continuation Hypercubes from Python

`fch_wrapper.FCH_CPU_Runner` loads the shared library `python/fch_lib/lib/liblib_fch_cpu.so` (Linux x86-64) and
exposes the four steps of a run (Section 10.6 of the book):

1. `find_first_simplex(ndim, kdim, math_func, domain_min, domain_max, domain_div, out_path, PRE_ROOT=..., max_iter_per_exponent=...)`
   finds a transversal $k$-simplex near the point `PRE_ROOT`, which should be near $\mathcal M$ and in general position
   with respect to the grid; `read_find_first_simplex_binary(out_path)` reads the result.
2. `run_threadpool(math_func, domain_min, domain_max, domain_div, first_coord, first_eidx, out_path, n_threads=0)`
   (or `run_batched`) traces the graph $\Gamma$ in parallel threads and writes its edges;
   `read_fch_output_binary(path)` reads them.
3. `skeleton(in_path, out_path)` builds the combinatorial skeleton in the format of [output_format.md](output_format.md).
4. `skeleton_simplicial_decomp(in_path, out_path)` subdivides the cells into simplices for rendering.

The map is a **plain Python function** `f(x, y)` that writes $F(x)$ into `y`; the interface compiles it with Numba
into a C callback, so the C++ threads call it without the Python interpreter. A compiled Numba `cfunc` (with an
`.address`) is also accepted, which is how the partition of unity implicits of Chapter 14 are passed
(`python/implicit_manifolds/learning/pui.py`, method `cfunc`).

Ready-made examples (`python/fch_lib/examples`): springs in 1, 2 and 3 dimensions, solid spheres and the
Cahn–Hilliard Galerkin system in $\mathbb R^{100}$ of Chapter 10:

```sh
cd python/fch_lib
python run_example.py my_math_spring2d spring2d       # results in out/spring2d/
```

The steps are also available as separate scripts `step1_compute_first_simplex.py` … `step5_plot_simplices.py`.

## `python/implicit_manifolds` — tools used throughout the book

| Module | Contents |
|---|---|
| `polfile.py` | reader of the `.pol` format: cells with their faces, cyclic order of 2-faces (`polygons`), identification of vertices (`global_complex`), Euler characteristic (`euler`) |
| `slicer.py` | sections of 3-manifolds by hyperplanes (fans from barycentres, cutting, chaining) |
| `zbuf.py` | small z-buffer renderer for polygon soups (`render`, `crop`), used for most pictures |
| `render.py` | matplotlib helpers for projections of complexes of $\mathbb R^n$ |
| `homology/` | homology of the combinatorial skeleton (Chapter 11): dual complexes, boundary case, convergence and ambiguity experiments |
| `j1a/` | the adaptive triangulation $J_1^a$ in Python, with a checker of the pivoting rules (Chapter 4) |
| `modeling/` | implicit modelling (Chapter 13): $p$-norm union and intersection (`uminp`, `imaxp`) and the skeletons of the knot and of the tube of $\mathbb R^4$; the models of the chapter are the C++ test problems `MOD_*` and the scripts of `book/ch13_modeling` |
| `learning/` | partition of unity implicits (Chapter 14): `pui.py` (reference), `pui_fast.py` (parallel Numba), `pui_rbf.py` (RBF and Hermite RBF local approximations), `project.py` (projection of a sample onto the polytopes), `runfch.py` (FCH on a PUI), `tda.py` (Betti numbers of a complex, Rips persistence) |

Scripts add these directories to `sys.path`; the modules are imported by their file names
(`from polfile import read_pol`, `from pui_fast import PUIFast`).
