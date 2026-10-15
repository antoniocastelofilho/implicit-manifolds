# Reproducing the results of the book

1. **Install**: `pip install -r requirements.txt`; a C++20 compiler, CMake and Eigen for the C++ programs; GNU Octave
   for the MATLAB codes (Table of interpreted codes only).
2. **C++ runs** (Chapters 5–6, 9–13): `book/bench/bench.sh <CASE>` builds a test problem and runs every method,
   `book/bench/run_case.sh <CASE> ...` runs Marching Simplex and the FCH pipeline (modelling cases `MOD_*`),
   `book/bench/conv.sh` and `book/bench/scale.sh` run the convergence and scaling studies. Results go to
   `book/runs/bench/` and `book/runs/modeling/`.
3. **Python runs** (Chapters 10, 14): the scripts of `book/ch09_ch10_hypercubes/` call FCH through `python/fch_lib`;
   `book/ch14_learning/experiments.py curve surface torus4 torus6 swiss swiss2` builds the partition of unity
   implicits and their approximations; `dimred.py`, `dimred_fail.py`, `bench_par.py`, `neural.py`,
   `neural_correct.py` and `rbf_compare.py` run the other experiments of the chapter.
4. **Figures and tables**: run the figure scripts listed in [../book/README.md](../book/README.md). Each figure script
   reads the results of the corresponding runs and writes to `book/figures/`.

The scripts were run with Python 3.13, NumPy 2.5, SciPy 1.18, Numba 0.68, Matplotlib 3.11, scikit-learn 1.9,
GUDHI 3.13 and umap-learn 0.5, GCC 13 and CMake 3.28. Random samples use fixed seeds. Timings depend on the machine;
those of the book were measured on a virtual machine with two cores (Octave runs: one core of an Intel Xeon at
2.8 GHz). The longest runs are Marching Simplex on the Bratu problem in $\mathbb R^6$ (about 20 minutes) and the
partition of unity of the torus of $\mathbb R^6$ (about one minute).
