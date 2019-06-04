import pandas as pd
from coolfuncs import *

sle = {'C3', 'O10', 'H11'}
sls = {'C2', 'O1', 'H10'}


def get_atoms(file='system.gro', res_num=15):
    system = read_gro(file)
    new_sys = []
    points_names = {'H5', 'H7', 'H8', 'H9', 'C6', 'O4'}
    points_values = []
    var_atoms = pd.DataFrame(np.zeros(shape=(6, 3)), columns=['x', 'y', 'z'], index=list(sls) + list(sle))
    for line in system:
        if line[0] == res_num + 1 and line[2] in sls:
            var_atoms.loc[line[2]] = line[4:]
        elif line[0] == res_num and line[2] in sle:
            var_atoms.loc[line[2]] = line[4:]
        else:
            new_sys.append(line)
            points_values.append(line[4:])
    points_values = pd.DataFrame(points_values, columns=['x', 'y', 'z'])
    return var_atoms, points_values, new_sys


def get_fake():
    return pd.DataFrame(np.random.rand(18).reshape(6, 3), columns=['x', 'y', 'z'], index=list(sls) + list(sle))


def make_zmatr(matr, atoms):
    for line in matr:
        if len(line) >= 3:
            line[2] = radius([atoms.loc[line[0]].values, atoms.loc[line[1]].values])
        if len(line) >= 5:
            line[4] = angle([atoms.loc[line[3]].values, atoms.loc[line[1]].values, atoms.loc[line[0]].values])
        if len(line) >= 7:
            line[6] = dihedral(np.array([atoms.loc[line[5]].values, atoms.loc[line[3]].values,
                                         atoms.loc[line[1]].values, atoms.loc[line[0]].values]))
    return matr


def z_matr2xyz(z_matrix):
    system = pd.DataFrame(np.zeros(shape=(len(sls) + len(sle), 3)), columns=['x', 'y', 'z'],
                          index=list(sls) + list(sle))
    system.loc[z_matrix[0][0]] = [0, 0, 0]
    system.loc[z_matrix[1][0]] = [z_matrix[1][2], 0, 0]
    system.loc[z_matrix[2][0]] = [z_matrix[1][2] - math.cos(z_matrix[2][4] / 180 * np.pi) * z_matrix[2][2],
                                  math.sin(z_matrix[2][4] / 180 * np.pi) * z_matrix[2][2], 0]
    for i in range(3, 6):
        system.loc[z_matrix[i][0]] = calc_coord([system.loc[z_matrix[i][1]], system.loc[z_matrix[i][3]],
                                                 system.loc[z_matrix[i][5]]], radius=z_matrix[i][2],
                                                angle1=z_matrix[i][4], angle2=z_matrix[i][6])
    return system


def lennard(mass):
    e = 1e-8
    q = 0.2
    energy = 4 * e * ((q / mass) ** 12 - (q / mass) ** 6)
    energy = np.sum(energy)
    return energy


def my_func(x, get_coord=False):
    z_matrix[2][-1] = x[0]
    z_matrix[3][-3] = x[1]
    z_matrix[3][-1] = x[2]
    z_matrix[4][-3] = x[3]
    z_matrix[4][-1] = x[4]
    z_matrix[5][-3] = x[5]
    z_matrix[5][-1] = x[6]
    xyz = np.array(z_matr2xyz(z_matrix))
    if get_coord:
        return z_matr2xyz(z_matrix)
    xyz = np.vstack((xyz, points))
    mass = np.sum((xyz[:, np.newaxis, :] - xyz[np.newaxis, :, :]) ** 2, axis=-1)
    return lennard(mass[np.triu_indices(len(mass), k=1)])


def writeres(xyz):
    xyz = go_back(xyz, params=params)

    with open('opt.gro', 'w') as file:
        count = 1
        file.write('!comment\n{}\n'.format(len(xyz.index) + len(system)))
        for i in xyz.index:
            s = '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % (1, 'SLS', i, count, *xyz.loc[i].values)
            count += 1
            file.write(s)
        for i in system:
            i[3] = count
            s = '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % (i[0], i[1], i[2], i[3], i[4], i[5], i[6])
            file.write(s)
            count += 1
        file.write('10 10 10\n')


var_atoms, points, system = get_atoms()
var_atoms, params = to_origin(var_atoms, 'C3', 'C2', 'O1')
points = to_origin(points, params=params)
for i in range(len(system)):
    system[i][4:] = points.loc[i]
z_matrix = [['C3'],
            ['C2', 'C3', 'Rcc'],
            ['O1', 'C2', 'Rco', 'C3', 'Acco'],
            ['H10', 'O1', 'Roh', 'C2', 'Acoh', 'C3', 'Dccoh'],
            ['O10', 'C3', 'Rco', 'C2', 'Acco', 'O1', 'Docco'],
            ['H11', 'O10', 'Roh', 'C3', 'Acoh', 'C2', 'Dccoh']]
z_matrix = make_zmatr(z_matrix, var_atoms)

if __name__ == '__main__':
    x0 = np.array([34.974239616979354, 119.89345671591812, 10.927791425828554,
                   23.794593254343177, -2.420176407853105, 136.34921807599954, -153.38899809503656])
    print(my_func(x0))
