import pandas as pd
import numpy as np
from coolfuncs import *
from scipy.optimize import minimize


def lennard(mass): # Func for calc energy
    e = 1e-8
    q = 0.2
    energy = 4 * e * ((q / mass) ** 12 - (q / mass) ** 6)
    energy = np.sum(energy)
    return energy


class Optimizator:
    def __init__(self, variable_atoms, points):
        self.sle = {'C3', 'O10', 'H11'}  # sets with names of atoms (variable) for optimization
        self.sls = {'C2', 'O1', 'H10'}
        self.xyzcoords = pd.DataFrame(np.zeros(shape=(len(self.sls) + len(self.sle), 3)), columns=['x', 'y', 'z'],
                                      index=list(self.sls) + list(self.sle))  # xyz coord for var atoms
        self.variable_atoms = pd.DataFrame(np.zeros(shape=(len(self.sls) + len(self.sle), 3)), columns=['x', 'y', 'z'],
                                           index=list(self.sls) + list(self.sle))
        self.points_value = []  # values with atoms, which are close to var atoms
        for line in variable_atoms:
            if line[2] in self.sls or line[2] in self.sle:
                self.variable_atoms.loc[line[2]] = line[4:]
        for i in points:
            self.points_value.append(i[4:])
        self.z_matrix = [['C3'],
                         ['C2', 'C3', 'Rcc'],
                         ['O1', 'C2', 'Rco', 'C3', 'Acco'],
                         ['H10', 'O1', 'Roh', 'C2', 'Acoh', 'C3', 'Dccoh'],
                         ['O10', 'C3', 'Rco', 'C2', 'Acco', 'O1', 'Docco'],
                         ['H11', 'O10', 'Roh', 'C3', 'Acoh', 'C2', 'Dccoh']]
        self.variable_atoms, self.params = to_origin(self.variable_atoms, 'C3', 'C2', 'O1')
        # move atoms to position for z matrix second param to 0,0,0 third to x, 0, 0 fourth to x, y, 0
        self.points_value = pd.DataFrame(self.points_value, columns=['x', 'y', 'z'])
        self.points_value = to_origin(self.points_value, params=self.params)
        self.opt_x = None  #

        for line in self.z_matrix:
            if len(line) >= 3:
                line[2] = radius([self.variable_atoms.loc[line[0]].values, self.variable_atoms.loc[line[1]].values])
                # set radius
            if len(line) >= 5:
                line[4] = angle([self.variable_atoms.loc[line[3]].values, self.variable_atoms.loc[line[1]].values,
                                 self.variable_atoms.loc[line[0]].values])
                # set angle between 3 atoms
            if len(line) >= 7:
                line[6] = dihedral(np.array([self.variable_atoms.loc[line[5]].values,
                                             self.variable_atoms.loc[line[3]].values,
                                             self.variable_atoms.loc[line[1]].values,
                                             self.variable_atoms.loc[line[0]].values]))
                # set dihedral atoms between 4 atoms

    def z_matrix2xyz(self):  # translate z matrix to xyz coordination
        self.xyzcoords.loc[self.z_matrix[0][0]] = [0, 0, 0]
        self.xyzcoords.loc[self.z_matrix[1][0]] = [self.z_matrix[1][2], 0, 0]
        self.xyzcoords.loc[self.z_matrix[2][0]] = [self.z_matrix[1][2] -
                                                   math.cos(self.z_matrix[2][4] / 180 * np.pi) * self.z_matrix[2][2],
                                                   math.sin(self.z_matrix[2][4] / 180 * np.pi) * self.z_matrix[2][2], 0]
        for i in range(3, 6):
            self.xyzcoords.loc[self.z_matrix[i][0]] = calc_coord([self.xyzcoords.loc[self.z_matrix[i][1]],
                                                                  self.xyzcoords.loc[self.z_matrix[i][3]],
                                                                  self.xyzcoords.loc[self.z_matrix[i][5]]],
                                                                 radius=self.z_matrix[i][2],
                                                                 angle1=self.z_matrix[i][4], angle2=self.z_matrix[i][6])

    def calc_energy(self, x):
        self.z_matrix[2][-1] = x[0]
        self.z_matrix[3][-3] = x[1]
        self.z_matrix[3][-1] = x[2]
        self.z_matrix[4][-3] = x[3]
        self.z_matrix[4][-1] = x[4]
        self.z_matrix[5][-3] = x[5]
        self.z_matrix[5][-1] = x[6]
        self.z_matrix2xyz()
        # set variable in z matrix and translate it to xyz
        xyz = np.array(self.xyzcoords)
        xyz = np.vstack((xyz, self.points_value))  # add points
        mass = np.sum((xyz[:, np.newaxis, :] - xyz[np.newaxis, :, :]) ** 2, axis=-1)
        return lennard(mass[np.triu_indices(len(mass), k=1)])  # return energy

    def optimaze(self):
        x0 = [self.z_matrix[2][-1],
              self.z_matrix[3][-3],
              self.z_matrix[3][-1],
              self.z_matrix[4][-3],
              self.z_matrix[4][-1],
              self.z_matrix[5][-3],
              self.z_matrix[5][-1]]
        # set start point
        res = minimize(self.calc_energy, x0, method='nelder-mead',
                       options={'xtol': 1e-2, 'disp': True})
        self.opt_x = res.x

    @property
    def opt_struc(self):
        self.z_matrix[2][-1] = self.opt_x[0]
        self.z_matrix[3][-3] = self.opt_x[1]
        self.z_matrix[3][-1] = self.opt_x[2]
        self.z_matrix[4][-3] = self.opt_x[3]
        self.z_matrix[4][-1] = self.opt_x[4]
        self.z_matrix[5][-3] = self.opt_x[5]
        self.z_matrix[5][-1] = self.opt_x[6]
        self.z_matrix2xyz()  # returns optimization structure is returned to begin coord
        return go_back(self.xyzcoords, params=self.params)

    def get_part(self, res_num):
        system = []
        for i in self.opt_struc.index:
            if i in self.sle:
                system.append([res_num, 'SLE', i, 0, *self.opt_struc.loc[i]])
        for i in self.opt_struc.index:
            if i in self.sls:
                system.append([res_num + 1, 'SLS', i, 0, *self.opt_struc.loc[i]])
        return system  # return optimization part of molecule


if __name__ == '__main__':
    system = read_gro('structures\\system.gro')
    sle = {'C3', 'O10', 'H11'}
    sls = {'C2', 'O1', 'H10'}
    var_atoms = []
    points = []
    for line in system:
        if line[0] == 16 and line[2] in sls:
            var_atoms.append(line)
        elif line[0] == 15 and line[2] in sle:
            var_atoms.append(line)
        else:
            points.append(line)
    opt = Optimizator(var_atoms, points)
    opt.z_matrix2xyz()
    opt.optimaze()
    system = opt.get_part()
    for i in system:
        print(i)
