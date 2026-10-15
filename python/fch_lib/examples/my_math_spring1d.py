import math
import numpy as np

# Dimensions.
ndim = 3    # domain dimension
kdim = 2    # codomain dimension (f: R^n -> R^k)

# Bounding box and grid.
dmin = np.full((ndim,), -1.4, dtype=np.float64)
dmin[2] = -1.
dmax = np.full((ndim,),  1.4, dtype=np.float64)
dmax[2] = 1
ddiv = np.full((ndim,), 17, dtype=np.uint32)
ddiv[2] = 41

# Pre-existing root to use as a starting point.
pre_root = None

# User math (Numba-friendly: use math.*, no Python objects)
def my_math(x, y):
    f = 1
    w = 2 * math.pi * f
    tr = 0.45
    y[0] = x[0] - math.cos(w * x[2]) * tr
    y[1] = x[1] - math.sin(w * x[2]) * tr