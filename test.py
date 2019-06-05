import numpy as np
from scipy.optimize import minimize

def rosen(x):
    """The Rosenbrock function"""
    return np.sum(100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0, axis=0)


x0 = np.array([2, 3, 4])
res = minimize(rosen, x0, method='nelder-mead',
               options={'xtol': 1e-9, 'disp': True})
print(res.x)
