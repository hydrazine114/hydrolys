import sys

import numpy as np
from math import sqrt
import math


def radius(p):
    return np.sum((p[0] - p[1]) ** 2) ** 0.5


def angle(p, isvector=False):
    vectors = np.array(p)
    if not isvector:
        vectors = np.array([p[0] - p[1], p[2] - p[1]])
    ch = np.sum(vectors[0] * vectors[1])
    zn = np.sum(vectors[0] ** 2) ** 0.5 * np.sum(vectors[1] ** 2) ** 0.5
    a = np.arccos(ch / zn)
    return a * 180 / np.pi


def dihedral(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return -np.degrees(np.arctan2(y, x))


def calc_coord(coordinations, radius=0.1, angle1=140, angle2=90):
    """
    Calculating xyz-coordinates from:
    :param coordinates: 3 atom in molec
    :param radius: distance from first atom and target atom
    :param angle1: angle target at---first at---secondat
    :param angle2: dihedral angle
    """
    coordinations = np.array(coordinations)
    na = 0
    nb = 1
    nc = 2
    z = np.zeros(3)
    x = coordinations[nc] - coordinations[nb]
    y = coordinations[na] - coordinations[nb]
    z[0] = x[1] * y[2] - x[2] * y[1]
    z[1] = x[2] * y[0] - x[0] * y[2]
    z[2] = x[0] * y[1] - x[1] * y[0]
    x[0] = y[1] * z[2] - y[2] * z[1]
    x[1] = y[2] * z[0] - y[0] * z[2]
    x[2] = y[0] * z[1] - y[1] * z[0]

    temp = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) ** 0.5
    help = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) ** 0.5
    work = (z[0] * z[0] + z[1] * z[1] + z[2] * z[2]) ** 0.5

    x = x / temp
    y = y / help
    z = z / work
    radh = np.pi / 180
    R = radius * np.sin(angle1 * radh)
    D = - radius * np.cos(angle1 * radh)
    E = R * np.cos(angle2 * radh)
    H = R * np.sin(angle2 * radh)

    vector = E * x + D * y + H * z
    result = coordinations[na] + vector
    return np.array(result)


def read_gro(input_file):
    system = []
    with open(input_file) as file:
        for line in file:
            if len(line) >= 44:
                system.append([int(line[0:5]), line[5:10].strip(),
                               line[10:15].strip(), int(line[15:20]), float(line[20:28]),
                               float(line[28:36]), float(line[36:44])])
    return system


def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = sqrt(mag2)
        v = tuple(n / mag for n in v)
    return np.array(v)


"""
Quat
"""


class Vector3f:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "Vector = " + str(self.x) + " " + str(self.y) + " " + str(self.z);

    @property
    def new_vector(self):
        return np.array([self.x, self.y, self.z])


class Quaternion:
    def __init__(self, w=1.0, x=0.0, y=0.0, z=0.0):
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "Quaternion = " + str(self.w) + " " + str(self.x) + " " + str(self.y) + " " + str(self.z)


def normal(v):
    normal = Vector3f()
    length = (v.x ** 2 + v.y ** 2 + v.z ** 2) ** 0.5
    normal.x = v.x / length
    normal.y = v.y / length
    normal.z = v.z / length
    return normal


def create_quat(rotate_vector, rotate_angle):
    quat = Quaternion()
    rotate_vector = normal(rotate_vector)
    quat.w = math.cos(rotate_angle / 2)
    quat.x = rotate_vector.x * math.sin(rotate_angle / 2)
    quat.y = rotate_vector.y * math.sin(rotate_angle / 2)
    quat.z = rotate_vector.z * math.sin(rotate_angle / 2)
    return quat


def quat_scale(q, val):
    q.w = q.w * val
    q.x = q.x * val
    q.y = q.y * val
    q.z = q.z * val
    return q


def quat_length(q):
    quat_length = (q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z) ** 0.5
    return quat_length


def quat_normalize(q):
    n = quat_length(q)
    return quat_scale(q, 1 / n)


def quat_invert(q):
    res = Quaternion()
    res.w = q.w
    res.x = -q.x
    res.y = -q.y
    res.z = -q.z
    return quat_normalize(res)


def quat_mul_quat(a, b):
    res = Quaternion()
    res.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    res.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y
    res.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x
    res.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w
    return res


def quat_mul_vector(a, b):
    res = Quaternion()
    res.w = -a.x * b.x - a.y * b.y - a.z * b.z
    res.x = a.w * b.x + a.y * b.z - a.z * b.y
    res.y = a.w * b.y - a.x * b.z + a.z * b.x
    res.z = a.w * b.z + a.x * b.y - a.y * b.x
    return res


def quat_transform_vector(q, v):
    t = Quaternion()
    t = quat_mul_vector(q, v)
    t = quat_mul_quat(t, quat_invert(q))
    ret = Vector3f()
    ret.x = t.x
    ret.y = t.y
    ret.z = t.z
    return ret.new_vector


"""
End quat
"""


def multy_vec(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def rotatepd(axial_vector, guide_vector, figure):
    """
    :param axial_vector: vector a is rotated to vector b
    :param guide_vector:
    :return: Rotated figure c like vector a
    """
    # figure = np.array(figure)
    an = angle([axial_vector, guide_vector], isvector=True) / 180 * math.pi
    v = multy_vec(axial_vector, guide_vector)
    qua = create_quat(Vector3f(*v), an)
    try:
        for i in figure.index:
            vec = Vector3f(*figure.loc[i])
            figure.loc[i] = quat_transform_vector(qua, vec)
    except Exception:
        for i in range(len(figure)):
            vec = Vector3f(figure[i])
            figure[i] = quat_transform_vector(qua, vec)

    return figure


def to_origin(var_atoms, first=None, second=None, third=None, params=None):
    if params is not None:
        var_atoms -= params[0]
        var_atoms = rotatepd(axial_vector=params[1], guide_vector=params[2], figure=var_atoms)
        var_atoms = rotatepd(axial_vector=params[3], guide_vector=params[4], figure=var_atoms)
        return var_atoms
    params = []
    rad = var_atoms.loc[first]
    params.append(rad.values.copy())
    var_atoms.loc[:] -= rad
    params.append(var_atoms.loc[second].values - var_atoms.loc[first].values)
    params.append([1, 0, 0])
    var_atoms = rotatepd(axial_vector=var_atoms.loc[second] - var_atoms.loc[first], guide_vector=[1, 0, 0],
                         figure=var_atoms)
    ax = list(var_atoms.loc[third])
    ax[0] = 0
    var_atoms = rotatepd(axial_vector=list(ax),
                         guide_vector=[0, 1, 0], figure=var_atoms)
    params.append(list(ax))
    params.append([0, 1, 0])
    return var_atoms, params


def go_back(var_atoms, params):
    var_atoms = rotatepd(axial_vector=params[-1], guide_vector=params[-2], figure=var_atoms)
    var_atoms = rotatepd(axial_vector=params[-3], guide_vector=params[-4], figure=var_atoms)
    var_atoms.loc[:] += params[-5]
    return var_atoms


if __name__ == '__main__':
    atoms = np.array([[0, 0, 0], [0.224464, 1.961145e-17, -7.596433e-19], [0.108703, 8.097870e-02, -1.734723e-17]])
    a = angle(atoms)
    print(a)
    r = radius(atoms[1:])
    print(0.224464 - math.cos(a / 180 * np.pi) * r, math.sin(a / 180 * np.pi) * r)
    pass
