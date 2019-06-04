import numpy as np
from scipy.optimize import minimize
from xyz_zmatrix import *


def rosen(x):
    """The Rosenbrock function"""
    return np.sum(100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0, axis=0)


x0 = np.array([34.974239616979354, 119.89345671591812, 10.927791425828554,
                  23.794593254343177, -2.420176407853105, 136.34921807599954, -153.38899809503656])
res = minimize(my_func, x0, method='nelder-mead',
               options={'xtol': 1e-2, 'disp': True})
print(res.x)
writeres(my_func(res.x, True))

"""
[  74.66041347  209.52697611  -55.77540883   64.80466944  -41.54420781
  157.43080618 -166.38363659]
"""
