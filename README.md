# implicit-manifolds

Companion code of the book

> **Approximation of Implicitly Defined Manifolds — Simplicial and Hypercube Methods Based on Counting and Enumeration**
> Antonio Castelo, Marcio Gameiro, Lucas Martinelli Reia, Manoel Jarra and Juliana Souza Paz. Springer, 2026.

The book studies the piecewise linear approximation of manifolds given implicitly,

$$\mathcal M = F^{-1}(0) = \{x\in\mathbb R^n : F(x)=0\},\qquad F\colon\mathbb R^n\to\mathbb R^k,\ n>k,$$

in any dimension $n$ and codimension $k$, by two families of algorithms: **simplicial methods**, which work on a
triangulation of the domain, and **hypercube methods**, which keep the cubes of a Cartesian grid. Both come in a
**marching** version, which visits every cell of the grid, and a **continuation** version, which starts from a point of
$\mathcal M$ and visits only the cells it crosses. Every object (simplex, cube, face, permutation) is an integer label
computed by **counting and enumeration**, so that the same code works in every dimension.

This repository contains all the implementations, the test problems and the scripts that produced every table and
figure of the book, together with the figures themselves.

| Algorithm | Book | Implementation |
|---|---|---|
| Marching Simplex on $K_1$, $J_1$ and the adaptive $J_1^a$ | Chap. 5 | C++ `cpp/simplicial`, MATLAB `matlab/marching_simplex` |
| Continuation Simplex on $K_1$, $J_1$ and $J_1^a$ | Chap. 6 | C++ `cpp/simplicial`, MATLAB `matlab/marching_simplex` |
| Permutahedral Tracing Algorithm (Boissonnat, Kachanovich, Wintraecken) | Chap. 7 | C++ `cpp/pta` |
| Combinatorial Marching / Continuation Hypercubes ($k=1$) | Chap. 8 | MATLAB `matlab/marching_hypercubes` |
| Generalized Combinatorial Marching / Continuation Hypercubes | Chap. 9 | MATLAB `matlab/generalized_marching_hypercubes`, C++ `cpp/gcch` |
| Fast Continuation Hypercubes (FCH) | Chap. 10 | C++ `cpp/fch`, Python `python/fch_lib` |
| Homology of the combinatorial skeleton | Chap. 11 | Python `python/implicit_manifolds/homology` |
| Implicit modelling (CSG, blends, skeletons, $\mathbb R^n$) | Chap. 13 | Python `python/implicit_manifolds/modeling`, C++ test cases `MOD_*` |
| Partition of unity implicits, manifold learning, RBF / Hermite RBF | Chap. 14 | Python `python/implicit_manifolds/learning` |

## Repository layout

```
implicit-manifolds/
├── cpp/                     C++20 library and programs (CMake)
│   ├── common/              grid, labels, counting, linear algebra, test problems (common/test_cases)
│   ├── simplicial/          Marching and Continuation Simplex on K1, J1 and J1a
│   ├── pta/                 permutahedral tracing algorithm and its skeleton
│   ├── gcch/                generalized combinatorial continuation hypercubes
│   ├── fch/                 Fast Continuation Hypercubes and post-processing
│   └── argparse/            command-line parser (header-only, MIT licence)
├── python/
│   ├── fch_lib/             Python interface of FCH (Numba callbacks, threads)
│   └── implicit_manifolds/  readers, Euler characteristic and homology, rendering, slicing,
│                            J1a, implicit modelling, partition of unity implicits
├── matlab/                  original MATLAB/Octave research codes
├── book/                    one directory per chapter with the scripts of the book,
│   ├── figures/             all the figures of the book (and of the companion article)
│   ├── tables/              all the tables of the book, as LaTeX
│   └── bench/               benchmark scripts of Chapter 12 and of the C++ test problems
└── docs/                    programs, Python interface, output format, test problems, reproduction
```

## Quick start

### C++ programs

Requirements: a C++20 compiler (GCC ≥ 11 or Clang ≥ 14), CMake ≥ 3.20 and Eigen 3 (`apt install libeigen3-dev`,
`brew install eigen`, or a copy in `cpp/eigen`). The problem is chosen at compile time, so that $F$ is inlined:

```sh
cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=BOOK_TORUS4
cmake --build build -j
cd build
./ms                     # Marching Simplex on K1
./cs torus_cs.pol        # Continuation Simplex, writes the skeleton
./fch -o fch.bin         # Fast Continuation Hypercubes (edges of the graph Gamma)
./fch_to_hypercubes -i fch.bin -o hc.bin
./reorder_hypercubes -i hc.bin -o hcr.bin
./skeleton_hypercubes -i hcr.bin -o torus_fch.pol -ff 15.6
```

See [docs/programs.md](docs/programs.md) for all programs and [docs/test_problems.md](docs/test_problems.md) for the
test problems of the book and how to add one.

### Python interface (no compilation of the map)

```sh
pip install -r requirements.txt
```

```python
import sys, numpy as np
sys.path.insert(0, 'python/fch_lib')
from fch_wrapper import FCH_CPU_Runner

def f(x, y):                                  # a torus of R^4
    y[0] = x[0]**2 + x[1]**2 - 1.0
    y[1] = x[2]**2 + x[3]**2 - 0.25

r = FCH_CPU_Runner('python/fch_lib/lib/')
r.find_first_simplex(ndim=4, kdim=2, math_func=f,
    domain_min=np.array([-1.1, -1.1, -0.51, -0.51]), domain_max=np.array([1.1, 1.1, 0.51, 0.51]),
    domain_div=np.array([10, 10, 10, 10], np.uint32), out_path='params.bin',
    PRE_ROOT=np.array([0.9553, 0.2955, 0.3824, 0.3221]))
p = r.read_find_first_simplex_binary('params.bin')
r.run_threadpool(math_func=f, domain_min=p['dmin'], domain_max=p['dmax'], domain_div=p['ddiv'],
    first_coord=p['first_coord'], first_eidx=p['first_eidx'], out_path='torus.bin')
r.skeleton(in_path='torus.bin', out_path='torus.pol')
```

and then, for the topology of the result,

```python
sys.path.insert(0, 'python/implicit_manifolds'); sys.path.insert(0, 'python/implicit_manifolds/learning')
from polfile import read_pol, euler
from tda import betti_of_cells
n, k, div, cells = read_pol('torus.pol')
print(euler(cells)[1], betti_of_cells(cells, 2))     # 0 [1, 2, 1]
```

See [docs/python_interface.md](docs/python_interface.md).

## Reproducing the book

Every figure and table of the book has a script in `book/`; [book/README.md](book/README.md) lists, for each figure,
the chapter, the script and what must be run before it, and [docs/reproducing.md](docs/reproducing.md) explains the
order of the runs. Scripts find the root of the repository by themselves (or through the environment variable
`IMPLICIT_MANIFOLDS_ROOT`), write intermediate results to `book/runs/` (not versioned) and figures to `book/figures/`.
Random samples use fixed seeds; timings depend on the machine (those of the book were measured on a virtual machine
with two cores).

## Output format

All programs write the combinatorial skeleton of the approximation in the same text format (`.pol`), described in
[docs/output_format.md](docs/output_format.md) and in Appendix A of the book; `python/implicit_manifolds/polfile.py`
reads it.

## Original implementations and credits

The algorithms were first implemented in separate repositories, whose code is included here (`matlab/`, `cpp/pta`,
`cpp/gcch`, `cpp/fch`):

* Marching and Continuation Simplex (MATLAB/Octave) — A. Castelo.
* Combinatorial Marching Hypercubes — A. Castelo, L. M. Bueno, M. Gameiro.
* Generalized Combinatorial Marching Hypercubes — A. Castelo, G. K. Nakassima, L. M. Bueno, M. Gameiro.
* Fast Continuation Hypercubes, GCCH and the permutahedral tracing algorithm (C++) — L. M. Reia, M. Gameiro,
  T. B. M. Ribeiro, A. Castelo (<https://github.com/lucasmreia/FCH>).
* Visualization of high-dimensional polytopal complexes: TrueNgine (G. S. Bicalho) and NDRenderer (T. B. M. Ribeiro).

Please cite the book and, for each algorithm, the paper in which it was published:

* A. Castelo, L. M. Bueno, M. Gameiro. A combinatorial marching hypercubes algorithm. *Computers & Graphics* 102 (2022) 67–77.
* A. Castelo, G. Nakassima, L. M. Bueno, M. Gameiro. A generalized combinatorial marching hypercube algorithm.
  *Computational and Applied Mathematics* 43 (2024) 127.
* L. M. Reia, M. Gameiro, T. B. M. Ribeiro, A. Castelo. A fast high-dimensional continuation hypercubes algorithm.
  *Computers & Graphics* 129 (2025) 104237.
* J.-D. Boissonnat, S. Kachanovich, M. Wintraecken. Tracing isomanifolds in $\mathbb R^d$ in time polynomial in $d$
  using Coxeter–Freudenthal–Kuhn triangulations. *SIAM J. Comput.* 52 (2023) 452–486.

## Licence

**To be defined by the authors before the repository is made public** (the original FCH repository asks for
citation of the paper; `cpp/argparse` is MIT-licensed; Eigen, used as an external dependency, is MPL-2.0).
