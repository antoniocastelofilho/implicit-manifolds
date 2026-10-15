# Runs of the MATLAB/Octave codes (Chapter 12 of the book)

`runjob(alg, n, k, First, Last, D, P0, f, fname)` runs one algorithm (`MS`, `CS`, `MHC`, `CHC`, `GMHC`, `GCHC`) on the
map `f`, writes the skeleton `fname` and appends the elapsed time to `timings.txt`. `jobsA.m` … `jobsK.m` are the batches
of the book (spheres, tori, complex curves, knot, $S^2\times S^1$, ...); `wrapF.m`, `knotF.m`, `toro8F.m` define maps used
by them. The algorithm files in this directory are the versions used for the runs (identical to those of
`matlab/`, except `MarchingSimplex.m`, taken from the GCMH repository).

```sh
cd book/common/octave_runs
octave --no-gui --eval "jobsA"        # then copy *.pol and timings.txt to book/runs/hypercubes/
python ../make_table.py               # writes book/tables/octave_runs.tex
```
