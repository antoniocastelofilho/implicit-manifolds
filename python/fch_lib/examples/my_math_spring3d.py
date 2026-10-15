import math
import numpy as np

# Dimensions.
ndim = 4    # domain dimension
kdim = 1    # codomain dimension (f: R^n -> R^k)

# Bounding box and grid.
dmin = np.full((ndim,), -1.4, dtype=np.float64)
dmin[2] = -1.
dmin[3] = 0.
dmax = np.full((ndim,),  1.4, dtype=np.float64)
dmax[2] = 1
dmax[3] = 1
ddiv = np.full((ndim,), 17, dtype=np.uint32)
ddiv[2] = 41
ddiv[3] = 8

# Pre-existing root to use as a starting point.
pre_root = None

# User math (Numba-friendly: use math.*, no Python objects)
def my_math(x, y):
    f = 1
    w = 2 * math.pi * f
    tr = 0.45
    ar = 0.75
    y[0] = math.sqrt(math.pow(x[0] - math.cos(w * x[2]) * tr, 2) + math.pow(x[1] - math.sin(w * x[2]) * tr, 2)) - (x[3] * ar)