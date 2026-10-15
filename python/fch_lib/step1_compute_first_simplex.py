import numpy as np
from fch_wrapper import FCH_CPU_Runner

# Importing the parameters of the domain and the implicit function that defines the manifold.
from my_math import *

runner = FCH_CPU_Runner("./lib/")

# Output file path.
out_path        = "out/params.bin"

# Call function to find the first simplex and write the result to a binary file.
runner.find_first_simplex(
    ndim=ndim,
    kdim=kdim,
    math_func=my_math,
    domain_min=dmin,
    domain_max=dmax,
    domain_div=ddiv,
    out_path=out_path,
    PRE_ROOT=pre_root,
    max_iter_per_exponent=3000
)
