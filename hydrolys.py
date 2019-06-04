import math
import numpy as np
from coolfuncs import *


class Hydrolis:
    def __init__(self, input_file, distant, num_of_wat, out_file):
        self.out_file = out_file
        self.atoms = []
        self.changed_system = []
        self.system = []
        self.not_water = []
        self.water = []
        self.input_file = input_file
        self.distant = distant
        self.distant = self.distant ** 2
        self.num_of_wat = num_of_wat

    def read_file(self):
        with open(self.input_file) as file:
            for line in file:
                if len(line) >= 44:
                    self.system.append([int(line[0:5]), line[5:10].strip(),
                                        line[10:15].strip(), int(line[15:20]), float(line[20:28]),
                                        float(line[28:36]), float(line[36:44])])
        return self.system

    def split_system(self):
        for line in self.system:
            if line[1] == 'SOL':
                self.water.append(line[4:])
            else:
                self.not_water.append(line)
        self.water = np.array(self.water)

    def choose_res(self):
        self.chosen = set()
        terminator = set()
        for atom in self.not_water:
            if atom[1] == 'SLM' and atom[2] == 'O1':
                if self.close2wat(atom[4:]):
                    self.chosen.add(atom[0])
            if atom[1] == 'SLS' or atom[1] == 'SLE':
                terminator.add(atom[0])

        # previous = -2
        # for i in self.chosen.copy():
        #     if i - 1 in terminator or i + 1 in terminator or previous + 1 == i:
        #         self.chosen.remove(i)
        #     else:
        #         previous = i

    def close2wat(self, atom):
        a = np.sum(np.sum((self.water - np.array(atom)) ** 2, axis=-1) < self.distant)
        return a >= self.num_of_wat

    def change_res(self):
        for line in self.system:
            if line[0] in self.chosen:  # to sle
                line[1] = 'SLE'
                self.changed_system.append(line)
                if line[2][0] == 'C':
                    self.atoms.append(line)
                if line[2] == 'H9':
                    reses = self.create_res()
                    for i in reses:
                        self.changed_system.append(i)
            elif line[0] - 1 in self.chosen:  # to sls
                line[1] = 'SLS'
                self.changed_system.append(line)
                if line[2] == 'C2' or line[2][0] == 'O':
                    self.atoms.append(line)
                if line[2] == 'H9':
                    reses = self.create_res()
                    for i in reses:
                        self.changed_system.append(i)
            else:
                self.changed_system.append(line)

    def create_res(self):
        if self.atoms[0][2] == 'C2':
            res = []
            distant = 0.09
            angle1 = 240
            angle2 = 60
            # angle1 = 0
            self.atoms[0], self.atoms[1] = self.atoms[1], self.atoms[0]
            coords = []
            for i in range(3):
                coords.append(self.atoms[i][4:])
            coords_o10 = calc_coord(coords, radius=distant, angle1=angle1, angle2=angle2)
            res.append([self.atoms[0][0], 'SLE', 'O10', 1, *coords_o10])
            coords[0] = coords_o10
            angle1 = -100
            coords_H11 = calc_coord(coords, radius=distant, angle1=angle1, angle2=angle2)
            res.append([self.atoms[0][0], 'SLE', 'H11', 1, *coords_H11])
        else:
            res = []
            distant = 0.09
            angle1 = 120
            angle2 = 10
            # angle1 = -60
            coords = []
            for i in range(3):
                coords.append(self.atoms[i][4:])
            coords_H10 = calc_coord(coords, radius=distant, angle1=angle1, angle2=angle2)
            res.append([self.atoms[0][0], 'SLS', 'H10', 1, *coords_H10])
        self.atoms = []
        return res

    def write(self):
        count = 1
        with open(self.out_file, 'w') as file:
            file.write('!comment\n{}\n'.format(len(self.changed_system)))
            for line in self.changed_system:
                # if line[0] <= 15:
                #     line[-3] += 0.4
                if 14 < line[0] <= 16:
                    line[3] = count
                    file.write('{:5d}{:5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(*line))
                    count += 1
            file.write('10 10 10\n')

    @property
    def ch_system(self):
        return self.changed_system.copy()


def optimization(system):
    pass


if __name__ == '__main__':
    input_file = 'polylac50.gro'
    hydro = Hydrolis(input_file, distant=0.3, num_of_wat=3, out_file='out.gro')
    hydro.read_file()
    hydro.split_system()
    hydro.choose_res()
    hydro.change_res()
    hydro.write()
