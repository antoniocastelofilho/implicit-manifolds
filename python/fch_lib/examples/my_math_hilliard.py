import math
import numpy as np

# Dimensions.
ndim = 100                  # domain dimension
kdim = ndim - 2             # codomain dimension (f: R^n -> R^k)

# Bounding box and grid.
dmin = np.full((ndim,), -0.3, dtype=np.float64)
dmin[0] = 39.0
dmax = np.full((ndim,),  0.3, dtype=np.float64)
dmax[0] = 43.0
ddiv = np.full((ndim,), 50, dtype=np.uint32)

# Pre-existing root to use as a starting point.
# pre_root = None
pre_root = np.array([3.9409e+01, -1.7640e-02, -1.6312e-02, -6.6755e-04, -1.4528e-04,
                     -9.9982e-04, -1.6667e-03, -2.0001e-03, -2.3333e-03, -2.6667e-03],
                    dtype=np.float64)

# User math (Numba-friendly: use math.*, no Python objects)
def my_math(x, y):
    """
    x layout:
        x[0] = lambda (lam)
        x[1] = c
        x[2:] = a[0 : la]  where la = NDIM - 2
    y layout:
        y[0 : la] = f[0 : la]
    """
    la = ndim - 2                  # la = NDIM - 2
    k = la                           # from C++: k = (2*la - 1 + 1)/2 -> la
    lfb = 2 * la - 1

    lam = x[0]
    c = x[1]

    # temporaries (match C++ sizes)
    fa = np.empty(2 * la - 1, dtype=np.float64)
    tb = np.empty((2 * k - 2) + lfb, dtype=np.float64)  # length = 4*la - 3
    aa = np.empty(k, dtype=np.float64)
    a3 = np.empty(k, dtype=np.float64)

    # ---- block 1 ----
    # build 'fa' with mirrored a[]
    # fa = [a[la-1], ..., a[1], a[0], a[0], a[1], ..., a[la-1]]
    for i in range(la):
        fa[i] = x[2 + (la - 1 - i)]      # a[la - 1 - i]
        fa[i + la - 1] = x[2 + i]        # a[i]

    # zero-pad ends of tb
    tb_len = (2 * k - 2) + lfb
    for i in range(k - 1):
        tb[i] = 0.0
        tb[tb_len - 1 - i] = 0.0

    # middle of tb gets mirrored a[]
    # index = abs((lfb - 1 - i) - (k - 1)) = abs((2*la - 2 - i) - (la - 1)) = abs(la - 1 - i)
    for i in range(lfb):
        idx = (lfb - 1 - i) - (k - 1)
        if idx < 0:
            idx = -idx
        tb[i + (k - 1)] = x[2 + idx]     # a[idx]

    # convolution to compute aa
    T = 2 * k - 1                        # = 2*la - 1
    for i in range(k - 1, 2 * k - 1):
        si = 0.0
        base = (T - 1) + i               # = (2*k - 2) + i
        for j in range(T):
            tb_idx = base - j            # (i - j) + (T - 1)
            si += tb[tb_idx] * fa[j]
        aa[i - (k - 1)] = si

    # ---- block 2 ----
    # refill tb middle with mirrored aa
    for i in range(lfb):
        idx = (lfb - 1 - i) - (k - 1)
        if idx < 0:
            idx = -idx
        tb[i + (k - 1)] = aa[idx]

    # convolution to compute a3
    for i in range(k - 1, 2 * k - 1):
        si = 0.0
        base = (T - 1) + i
        for j in range(T):
            tb_idx = base - j
            si += tb[tb_idx] * fa[j]
        a3[i - (k - 1)] = si

    # ---- outputs ----
    # f[0] = (a[0] - a3[0]) - c
    y[0] = (x[2 + 0] - a3[0]) - c

    # for i >= 1:
    # f[i] = (1 - ((pi^2 / lambda) * i^2)) * a[i] - a3[i]
    pi2_over_lambda = (math.pi * math.pi) / lam
    for i in range(1, la):
        y[i] = (1.0 - (pi2_over_lambda * (i * i))) * x[2 + i] - a3[i]