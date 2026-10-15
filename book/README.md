# Scripts, figures and tables of the book

Each directory holds the scripts of one chapter (Chapters 9 and 10 share one). Scripts locate the root of the
repository by themselves (or use `IMPLICIT_MANIFOLDS_ROOT`), write intermediate results to `book/runs/` (not
versioned) and figures to `book/figures/`. See [../docs/reproducing.md](../docs/reproducing.md) for the order of the runs.

| Directory | Contents |
|---|---|
| `ch01_introduction/` | gallery, $S^2\times S^1$ and its sections, growth of the cost, cells of MS vs CMH, ambiguity |
| `ch02_manifolds/` | the graph $\Gamma$ in the plane |
| `ch04_triangulations/` | $K_1$, $J_1$, Coxeter–Freudenthal–Kuhn, pivoting, $J_1^a$ in 2D and on a torus |
| `ch05_marching_simplex/` | perturbation, cells, spheres, $S^3$, torus, $J_1^a$ |
| `ch06_continuation_simplex/` | barycentric walk, door-in/door-out, complex curves, $J_1^a$ |
| `ch08_cmh/` | ambiguous squares, generalized hyperbola |
| `ch09_ch10_hypercubes/` | GCMH/GCCH and FCH examples up to $\mathbb R^{10}$, re-rendering, tests of FCH; data in `.pkl` |
| `ch11_homology/` | duals, rates and ambiguity of the homology study |
| `ch12_comparison/` | strategies, Gaussian blobs, Bratu problems, scaling, convergence, cells, table of strengths |
| `ch13_modeling/` | 2D and $p$-norm figures, models in $\mathbb R^3$ and $\mathbb R^4$, table of runs |
| `ch14_learning/` | partition of unity implicits: experiments, dimensionality reduction, TDA, parallel and learned versions, RBF/Hermite; `paper_simods/` has the scripts of the companion article |
| `common/` | scripts shared by several chapters (`fig_main.py`, `fig_sphere3.py`, `make_table.py`) |
| `bench/` | build-and-run scripts of the C++ programs, parsing of the logs, tables of Chapter 12 |
| `figures/` | every figure of the book; `figures/paper/` those of the article, `figures/extra/` figures not used in the final text |
| `tables/` | every table of the book, as included by LaTeX |

## Figures

| Figure file | Chapter | Script |
|---|---|---|
| `intro_ambiguity.png` | Ch. 01 | `ch01_introduction/fig_ambig.py` |
| `intro_cells.png` | Ch. 01 | `ch01_introduction/fig_cells.py` |
| `intro_gallery.png` | Ch. 01 | `ch01_introduction/fig_gallery.py` |
| `intro_growth.pdf` | Ch. 01 | `ch01_introduction/fig_growth.py` |
| `intro_s2s1_global.png` | Ch. 01 | `ch01_introduction/fig_s2s1.py` |
| `intro_s2s1_sections.png` | Ch. 01 | `ch01_introduction/fig_s2s1.py` |
| `j1a_2d.pdf` | Ch. 01, Ch. 04 | `ch04_triangulations/fig_j1a.py` |
| `tri_k1j1.pdf` | Ch. 01, Ch. 04 | `ch04_triangulations/fig4.py` |
| `gamma_2d.pdf` | Ch. 02 | `ch02_manifolds/gamma2d.py` |
| `j1a_torus.png` | Ch. 04 | `ch04_triangulations/fig_torus.py` |
| `tri_cfk.pdf` | Ch. 04 | `ch04_triangulations/fig4.py` |
| `tri_pivj1.pdf` | Ch. 04 | `ch04_triangulations/fig4.py` |
| `tri_pivk1.pdf` | Ch. 04 | `ch04_triangulations/fig4.py` |
| `ms_cells.pdf` | Ch. 05 | `ch05_marching_simplex/fig5.py` |
| `ms_j1a_s3.png` | Ch. 05 | `ch05_marching_simplex/fig5.py` |
| `ms_pert.pdf` | Ch. 05 | `ch05_marching_simplex/fig5.py` |
| `ms_s3.png` | Ch. 05 | `ch05_marching_simplex/fig5.py` |
| `ms_sphere_k1j1.png` | Ch. 05 | `ch05_marching_simplex/fig5.py` |
| `ms_torus.png` | Ch. 05 | `ch05_marching_simplex/fig5.py` |
| `complex_cs.png` | Ch. 06 | `ch06_continuation_simplex/fig6.py`, `common/fig_main.py` |
| `cs_baric.pdf` | Ch. 06 | `ch06_continuation_simplex/fig6.py` |
| `cs_j1a_zcosw.png` | Ch. 06 | `ch06_continuation_simplex/fig6.py` |
| `cs_walk.pdf` | Ch. 06 | `ch06_continuation_simplex/fig6.py` |
| `cmh_ambig.pdf` | Ch. 08 | `ch08_cmh/fig_ambig.py` |
| `cmh_hyperbola.png` | Ch. 08 | `ch08_cmh/fig_hyp.py` |
| `s3_cmh.png` | Ch. 08 | `ch09_ch10_hypercubes/fig_rerender.py`, `common/fig_main.py` |
| `sphere_ms_cmh.png` | Ch. 08 | `ch09_ch10_hypercubes/fig_rerender.py`, `common/fig_sphere3.py` |
| `complex_gcch.png` | Ch. 09 | `ch09_ch10_hypercubes/fig_rerender.py`, `common/fig_main.py` |
| `hc_s2_r10.png` | Ch. 09 | `ch09_ch10_hypercubes/fig_ch9.py` |
| `hc_s2s1_sections.png` | Ch. 09 | `ch09_ch10_hypercubes/fig_ch9.py` |
| `hc_t3_section.png` | Ch. 09 | `ch09_ch10_hypercubes/fig_ch9.py` |
| `hc_torus_r10.png` | Ch. 09 | `ch01_introduction/fig_gallery.py`, `ch09_ch10_hypercubes/fig_ch9.py` |
| `klein.png` | Ch. 09 | `ch01_introduction/fig_gallery.py`, `ch09_ch10_hypercubes/fig_ch9.py` |
| `torus4_hypercube.png` | Ch. 09 | `ch09_ch10_hypercubes/fig_rerender.py`, `common/fig_main.py` |
| `cahn_hilliard100.png` | Ch. 10 | `ch01_introduction/fig_gallery.py`, `ch09_ch10_hypercubes/fig_fch.py` |
| `fch_tests.png` | Ch. 10 | `ch09_ch10_hypercubes/fig_ch10.py`, `ch09_ch10_hypercubes/fig_fch_tests.py` |
| `springs.png` | Ch. 10 | `common/fig_main.py` |
| `hom_ambiguity.pdf` | Ch. 11 | `ch11_homology/fig_hom.py` |
| `hom_dual2d.pdf` | Ch. 11 | `ch11_homology/fig_hom.py` |
| `hom_dual_torus.png` | Ch. 11 | `ch11_homology/fig_hom.py` |
| `hom_rate.pdf` | Ch. 11 | `ch11_homology/fig_hom.py` |
| `cmp_blobs.png` | Ch. 12 | `ch12_comparison/fig_blobs_bratu.py` |
| `cmp_bratu.pdf` | Ch. 12 | `ch12_comparison/fig_blobs_bratu.py` |
| `cmp_cells.pdf` | Ch. 12 | `ch12_comparison/fig_conv.py` |
| `cmp_conv.pdf` | Ch. 12 | `ch12_comparison/fig_conv.py` |
| `cmp_scale.pdf` | Ch. 12 | `ch12_comparison/fig_scale.py` |
| `cmp_strategies.pdf` | Ch. 12 | `ch12_comparison/fig_strategies.py` |
| `eight_tori.png` | Ch. 13 | `ch13_modeling/fig_3d.py`, `common/fig_main.py` |
| `intersection.png` | Ch. 13 | `ch13_modeling/fig_2d.py`, `ch13_modeling/fig_3d.py` |
| `knot.png` | Ch. 13 | `ch01_introduction/fig_gallery.py`, `ch13_modeling/fig_3d.py` |
| `mod_2d.pdf` | Ch. 13 | `ch13_modeling/fig_2d.py` |
| `mod_pnorm.pdf` | Ch. 13 | `ch13_modeling/fig_2d.py` |
| `mod_r4_balls.png` | Ch. 13 | `ch13_modeling/fig_r4.py` |
| `mod_r4_diff.png` | Ch. 13 | `ch13_modeling/fig_r4.py` |
| `mod_r4_tube.png` | Ch. 13 | `ch01_introduction/fig_gallery.py`, `ch13_modeling/fig_r4.py` |
| `pipebox.png` | Ch. 13 | `ch01_introduction/fig_gallery.py`, `ch13_modeling/fig_3d.py` |
| `two_spheres.png` | Ch. 13 | `ch13_modeling/fig_3d.py`, `common/fig_main.py` |
| `ml_dimred_fail.png` | Ch. 14 | `ch14_learning/dimred_fail.py` |
| `ml_pui2d.pdf` | Ch. 14 | `ch14_learning/figs.py` |
| `ml_rbf.pdf` | Ch. 14 | `ch14_learning/rbf_fig.py` |
| `ml_surface.png` | Ch. 14 | `ch14_learning/figs.py` |
| `ml_swiss.png` | Ch. 14 | `ch14_learning/figs.py` |
| `ml_tda.pdf` | Ch. 14 | `ch14_learning/figs.py` |
| `ml_tori.png` | Ch. 14 | `ch14_learning/figs.py` |

## Tables

| Table file | Chapter | Script |
|---|---|---|
| `octave_runs.tex` | Ch. 12 | runs of the MATLAB/Octave codes (`matlab/`), compiled by `common/make_table.py` |
| `cpp_runs.tex` | Ch. 12 | `bench/bench.sh` on the `BOOK_*` cases, then `bench/bench_table.py` |
| `cmp_strengths.tex` | Ch. 12 | `bench/bench.sh` on `BOOK_BLOBS*` and `BOOK_BRATU*`, then `ch12_comparison/table_strengths.py` |
| `mod_runs.tex` | Ch. 13 | `bench/run_case.sh MOD_...`, then `ch13_modeling/table.py` |
| `ml_runs.tex` | Ch. 14 | `ch14_learning/experiments.py`, then `ch14_learning/table.py` |
| `ml_rbf.tex` | Ch. 14 | `ch14_learning/rbf_compare.py`, then `ch14_learning/rbf_fig.py` |
