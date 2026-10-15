import numpy as np
from fch_wrapper import FCH_CPU_Runner

# Importing *only* the implicit function that defines the manifold.
from my_math import my_math

runner = FCH_CPU_Runner("./lib/")

# Load parameters calculated by "find first simplex".
params = runner.read_find_first_simplex_binary("out/params.bin")


# Output file path.
out_path        = "out/results.bin"

# Run FCH with threadpool.
runner.run_threadpool(
    math_func=my_math,
    domain_min=params["dmin"],
    domain_max=params["dmax"],
    domain_div=params["ddiv"],
    first_coord=params["first_coord"],
    first_eidx=params["first_eidx"],
    out_path=out_path,
    print_interval=100
)

# Run batched FCH.
# runner.run_batched(
#     math_func=my_math,
#     domain_min=params["dmin"],
#     domain_max=params["dmax"],
#     domain_div=params["ddiv"],
#     first_coord=params["first_coord"],
#     first_eidx=params["first_eidx"],
#     out_path=out_path,
# )
