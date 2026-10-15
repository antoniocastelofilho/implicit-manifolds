# Python code

* `fch_lib/` — the Python interface of Fast Continuation Hypercubes (`fch_wrapper.py`) with the prebuilt shared library
  `lib/liblib_fch_cpu.so` (Linux x86-64) and examples.
* `implicit_manifolds/` — the tools used throughout the book: reading the output format, Euler characteristic and
  homology, rendering and slicing, the adaptive triangulation $J_1^a$, implicit modelling and the partition of unity
  implicits of Chapter 14.

See [../docs/python_interface.md](../docs/python_interface.md). Install the dependencies with
`pip install -r ../requirements.txt`.
