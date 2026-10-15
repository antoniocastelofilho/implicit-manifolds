# C++ implementations

C++20 implementations of all the algorithms of the book, sharing the grid, the labels, the test problems and the
output format. The permutahedral tracing, GCCH and FCH codes come from <https://github.com/lucasmreia/FCH>
(L. M. Reia; original README in `README_FCH_original.md`); the simplicial codes (`simplicial/`) were written for the
book in the same framework.

```
common/            definitions, utility, bmi2 (PDEP/PEXT, with the zp7 fallback), test_cases/
simplicial/        simplicial.h (K1, J1, cells of a simplex), j1a.h (adaptive J1a), main_ms.cpp, main_cs.cpp,
                   main_j1a.cpp, main_cs_j1a.cpp
pta/               permutahedral tracing algorithm and skeleton_pta
gcch/              generalized combinatorial continuation hypercubes, reorder and skeleton of hypercubes
fch/               fast continuation hypercubes and fch_to_hypercubes
argparse/          https://github.com/p-ranav/argparse (MIT)
```

## Build

```sh
sudo apt install cmake g++ libeigen3-dev        # or: brew install cmake eigen
cmake -S cpp -B build -DCMAKE_BUILD_TYPE=Release -DTEST_CASE=BOOK_SPHERE3
cmake --build build -j
```

If `cpp/eigen` contains a copy of Eigen (`git submodule add https://gitlab.com/libeigen/eigen cpp/eigen`), it is used
instead of the system one. Release builds use `-O3 -march=native`.

Programs: see [../docs/programs.md](../docs/programs.md). Test problems: [../docs/test_problems.md](../docs/test_problems.md).
