# C++ programs

All programs are built from `cpp/` with CMake. The test problem is selected at compile time with
`-DTEST_CASE=<NAME>` (see [test_problems.md](test_problems.md)); the map $F$ is then inlined in the inner loops.

```sh
cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=BOOK_S2XS1
cmake --build build -j            # or: --target ms cs fch ...
```

| Program | Algorithm (book chapter) | Usage | Output |
|---|---|---|---|
| `ms` | Marching Simplex on $K_1$ (5) | `ms [out.pol]` | skeleton (`.pol`), counts of cubes, simplices, crossed simplices, evaluations of $F$ |
| `ms_j1` | Marching Simplex on $J_1$ (5) | `ms_j1 [out.pol]` | idem |
| `cs` | Continuation Simplex on $K_1$ (6) | `cs [out.pol\|-] [log2 hash size]` | skeleton of the component through `FIRST_POINT` |
| `cs_j1` | Continuation Simplex on $J_1$ (6) | idem | idem |
| `ms_j1a` | Marching Simplex on the adaptive $J_1^a$ (4, 5) | `ms_j1a <levels> [out.pol\|-] [check]` | skeleton; `check` verifies the pivoting rules |
| `cs_j1a` | Continuation Simplex on $J_1^a$ (6) | `cs_j1a <levels> [out.pol\|-]` | skeleton |
| `pta` | Permutahedral tracing algorithm (7) | `pta -o out.bin` | crossed simplices (binary) |
| `skeleton_pta` | skeleton of the PTA output | `skeleton_pta -i out.bin -o out.pol -ff 15.6` | skeleton |
| `gcch` | Generalized Combinatorial Continuation Hypercubes (9) | `gcch -o out.bin` | crossed cubes (binary) |
| `fch` | Fast Continuation Hypercubes (10) | `fch -o fch.bin` | edges of the graph $\Gamma$ (binary) |
| `fch_to_hypercubes` | edges of $\Gamma$ grouped by cube | `fch_to_hypercubes -i fch.bin -o hc.bin` | binary |
| `reorder_hypercubes` | cubes sorted by label | `reorder_hypercubes -i hc.bin -o hcr.bin` | binary |
| `skeleton_hypercubes` | combinatorial skeleton (8–11) | `skeleton_hypercubes -i hcr.bin -o out.pol -ff 15.6` | skeleton |

The `-ff` option sets the number of significant digits written for the coordinates. Every program prints its running
time and its counts (cells, simplices, evaluations of $F$); `book/bench/bench.sh` collects them for all methods.

## The four strategies

| | marching (whole grid) | continuation (from `FIRST_POINT`) |
|---|---|---|
| **simplices** | `ms`, `ms_j1`, `ms_j1a` | `cs`, `cs_j1`, `cs_j1a`, `pta` |
| **hypercubes** | MATLAB `matlab/marching_hypercubes`, `matlab/generalized_marching_hypercubes` | `gcch`, `fch` (+ MATLAB versions) |

Marching finds every component and needs no starting point; continuation costs in proportion to the output and
follows the component through `FIRST_POINT` (Chapter 12 of the book).

## Benchmarks

```sh
book/bench/bench.sh BOOK_TORUS4          # builds the case and runs every method, results in book/runs/bench/BOOK_TORUS4
book/bench/scale.sh                      # scaling with the grid (Fig. 12.4)
book/bench/conv.sh                       # convergence of the vertices (Fig. 12.6)
book/bench/run_case.sh MOD_KNOT ...      # MS + FCH pipeline for the modelling cases of Chapter 13
python book/bench/parse.py               # parses the logs; bench_table.py writes tables/cpp_runs.tex
```
