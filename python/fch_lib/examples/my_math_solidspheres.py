import math
import numpy as np

# Dimensions.
ndim = 4    # domain dimension
kdim = 1    # codomain dimension (f: R^n -> R^k)

# Bounding box and grid.
dmin = np.full((ndim,), -1.3, dtype=np.float64)
dmin[3] = 0.8
dmax = np.full((ndim,),  1.3, dtype=np.float64)
dmax[3] = 1
ddiv = np.full((ndim,), 10, dtype=np.uint32)
ddiv[3] = 3

# Pre-existing root to use as a starting point.
pre_root = None

# User math (Numba-friendly: use math.*, no Python objects)
def my_math(x, y):
    c1 = np.array([1,0,0], dtype=np.float64)
    c2 = np.array([-1,0,0], dtype=np.float64)
    y[0] = min((np.linalg.norm(x[:3]-c1) - x[3]), (np.linalg.norm(x[:3]-c2) - x[3]))