# Test problems

A test problem is a specialization of `TestCaseImpl` in `cpp/common/test_cases/`, selected at compile time with
`-DTEST_CASE=<NAME>`. It gives the dimensions, the box, the grid, a starting point for the continuation methods and
the map:

```cpp
template <>
struct TestCaseImpl<TestCase::BOOK_SPHERE3> {
  static constexpr size_t NDIM = 3, KDIM = 1;
  static constexpr std::array<double, NDIM> DOMAIN_MIN{-1.1, -1.1, -1.1};
  static constexpr std::array<double, NDIM> DOMAIN_MAX{ 1.1,  1.1,  1.1};
  static constexpr std::array<size_t, NDIM> DOMAIN_DIV{12, 12, 12};
  static constexpr std::array<double, NDIM> FIRST_POINT{0.6, 0.48, 0.64};
  static void func(const double *x, double *f) {
    f[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1.0;
  }
};
```

To add a problem: write the specialization in a header of `cpp/common/test_cases/` (included by `all_cases.h`) and
add its name to the `TestCase` enumeration of `test_case_template.h`. The starting point should be near the
manifold and not on a hyperplane of the grid; asymmetric boxes avoid degenerate positions.

## Problems of the book

| Name | n | k | grid | description | header |
|---|---|---|---|---|---|
| `BOOK_SPHERE3` | 3 | 1 | 12 × 12 × 12 |  | `book_cases.h` |
| `BOOK_SPHERE4` | 4 | 1 | 8 × 8 × 8 × 8 |  | `book_cases.h` |
| `BOOK_TORUS4` | 4 | 2 | 10 × 10 × 10 × 10 |  | `book_cases.h` |
| `BOOK_ZCOSW` | 4 | 2 | 10 × 10 × 10 × 10 |  | `book_cases.h` |
| `BOOK_Z2COSW` | 4 | 2 | 10 × 10 × 10 × 10 |  | `book_cases.h` |
| `BOOK_Z2W5` | 4 | 2 | 12 × 12 × 12 × 12 |  | `book_cases.h` |
| `BOOK_S2XS1` | 5 | 2 | 12 × 12 × 12 × 8 × 8 |  | `book_cases.h` |
| `BOOK_AMBIG_SPLIT` | 3 | 1 | 1 × 1 × 1 |  | `book_cases.h` |
| `BOOK_AMBIG_TUBE` | 3 | 1 | 1 × 1 × 1 |  | `book_cases.h` |
| `BOOK_J1A_ELLIPSE` | 2 | 1 | 4 × 4 | ellipse and a line (Fig. 4.5) | `book_cases.h` |
| `BOOK_J1A_TORUS` | 3 | 1 | 8 × 8 × 3 | torus refined around a wavy curve on it | `book_cases.h` |
| `BOOK_J1A_S3` | 4 | 1 | 4 × 4 × 4 × 4 | S^3 in R^4 refined around a 2-sphere | `book_cases.h` |
| `BOOK_J1A_ZCOSW` | 4 | 2 | 5 × 5 × 5 × 5 | z = cos w refined around the real curve Im w = 0.05 | `book_cases.h` |
| `MOD_TWO_SPHERES` | 3 | 1 | 40 × 28 × 28 | sharp union of two spheres | `book_cases.h` |
| `MOD_TWO_SPHERES_P8` | 3 | 1 | 40 × 28 × 28 | smooth union, p = 8 | `book_cases.h` |
| `MOD_EIGHT_TORI` | 3 | 1 | 80 × 80 × 20 | smooth union of eight tori, p = 10 | `book_cases.h` |
| `MOD_KNOT` | 3 | 1 | 72 × 72 × 42 | smooth union of 151 spheres on a knot, p = 30 | `book_cases.h` |
| `MOD_INTERSECT3` | 3 | 1 | 64 × 64 × 44 | smooth intersection of three spheres, p = 10 | `book_cases.h` |
| `MOD_R4_TWO_BALLS` | 4 | 1 | 36 × 22 × 22 × 22 | smooth union of two 3-spheres in R^4, p = 8 | `book_cases.h` |
| `MOD_R4_DIFF` | 4 | 1 | 32 × 32 × 32 × 32 | 4-ball minus the solid cylinder x1^2+x2^2 < 0.45^2, p-norm difference (p = 6) | `book_cases.h` |
| `MOD_R4_TUBE` | 4 | 1 | 44 × 44 × 22 × 22 | smooth union of 200 3-spheres along a closed curve of R^4, p = 30 | `book_cases.h` |
| `MOD_PIPEBOX` | 3 | 1 | 56 × 56 × 112 |  | `book_cases.h` |
| `MOD_PIPEBOX_P8` | 3 | 1 | 56 × 56 × 112 |  | `book_cases.h` |
| `BOOK_HYPERBOLA4` | 4 | 1 | 9 × 9 × 9 × 9 | generalized hyperbola x1 x2 x3 x4 = 1e-4 (GCMH article) | `book_cases.h` |
| `BOOK_T3_R6_B` | 6 | 3 | 8 × 8 × 8 × 8 × 8 × 8 |  | `book_cases.h` |
| `BOOK_CONV_S2` | 3 | 1 | CONV_N × CONV_N × CONV_N |  | `book_cases.h` |
| `BOOK_CONV_T2` | 4 | 2 | CONV_N × CONV_N × CONV_N/2 × CONV_N/2 | torus in R^4: x1^2+x2^2=1, x3^2+x4^2=1/4 | `book_cases.h` |
| `BOOK_BRATU3` | 4 | 3 | 24 × 20 × 20 × 30 | discrete Bratu problem, x=(lambda,u_1..u_3) | `book_cases_ch12.h` |
| `BOOK_BRATU5` | 6 | 5 | 12 × 10 × 10 × 12 × 10 × 10 | discrete Bratu problem, x=(lambda,u_1..u_5) | `book_cases_ch12.h` |
| `BOOK_BLOBS3` | 3 | 1 | 40 × 40 × 40 | 14 Gaussian blobs, level 1/2 | `book_cases_ch12.h` |
| `BOOK_BLOBS4` | 4 | 1 | 24 × 24 × 24 × 24 | 12 Gaussian blobs, level 1/2 | `book_cases_ch12.h` |
| `CIRCLE_9TO8` | 9 | 8 | 11 × 11 × 11 × 11 × 11 × 11 × 11 × 11 × 11 | a circle in R^9 (k = 8) | `circle_9to8.h` |
| `HYPERSHPERE_5TO1` | 5 | 1 | 20 × 20 × 20 × 20 × 20 | sphere S^4 in R^5 | `hypersphere_5to1.h` |
| `HYPERSHPERE_N_TO_K` | 7 | 2 | see header | products of spheres, n and k set in the header | `hypersphere_n_to_k.h` |
| `KLEIN_5TO3` | 5 | 3 | 20 × 20 × 20 × 20 × 20 | graph of a parametrized Klein bottle in R^5 (original FCH repository) | `klein_5to3.h` |
| `S1_K_TIMES` | – | 4 | see header | torus (S^1)^k, k set in the header | `s1_k_times.h` |

`BOOK_CONV_S2` and `BOOK_CONV_T2` take the number of cubes per direction from the macro `CONV_N`
(`-DCMAKE_CXX_FLAGS=-DCONV_N=32`); they are used by `book/bench/conv.sh` and `book/bench/scale.sh`.
`book/ch12_comparison/gen_cases.py` generates `book_cases_ch12.h` (Bratu problems and Gaussian blobs).
