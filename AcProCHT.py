import math
import random
import sys


def calc_degrees(system, ad=0., pd=0.):
    '''
    Calculating protonation and acetylation degrees with considering terminal residues
    '''
    cht0 = 0.  # 17
    chtn = 0.  # 16
    chtr = 0.  # 14
    chtp = 0.  # 16
    cht = 0.  # 15
    ace = 0  # 3
    for atom in system:
        if atom[1] == 'CHT0':
            cht0 = cht0 + 1.
        elif atom[1] == 'CHTR':
            chtr = chtr + 1.
        elif atom[1] == 'CHT':
            cht = cht + 1.
        elif atom[1] == 'CHTN':
            chtn = chtn + 1.
        elif atom[1] == 'CHTP':
            chtp = chtp + 1.
        elif atom[1] == 'ACE2':
            ace += 1

    if cht % 15 == 0 and cht0 % 17 == 0 and chtn % 16 == 0 and chtp % 16 == 0 and chtr//14 == ace // 3:
        pass
    else:
        print('Fatal Error!!! Some residues lost atom')
        print(chtp % 16, chtr / 14, ace / 3)
    max_degree = (cht / 15 + chtp / 16 + chtr / 14) / (cht / 15 + chtp / 16 + chtr / 14 + cht0 / 17 + chtn / 16)
    if ad + pd > max_degree:
        print(
            'System has only {:.1f}% changeable residues. '
            'Script will change less residues than you requested.'.format(max_degree * 100))
        return ad / max_degree, 1
    if ad + pd == max_degree:
        return ad / max_degree, 1.
    return ad / max_degree, pd / (max_degree - ad)


def print_info(system):
    cht0 = 0.  # 17
    chtn = 0.  # 16
    chtr = 0.  # 14
    chtp = 0.  # 16
    cht = 0.  # 15
    for atom in system:
        if atom[1] == 'CHT0':
            cht0 = cht0 + 1.
        elif atom[1] == 'CHTR':
            chtr = chtr + 1.
        elif atom[1] == 'CHT':
            cht = cht + 1.
        elif atom[1] == 'CHTN':
            chtn = chtn + 1.
        elif atom[1] == 'CHTP':
            chtp = chtp + 1.
    cht0 = cht0 / 17.
    chtn = chtn / 16.
    chtr = chtr / 14
    chtp = chtp / 16.
    cht = cht / 15
    allres = cht0 + chtn + chtr + chtp + cht
    prot_degree = chtp / allres
    acil_d = chtr / allres
    print(' Chitosan residues: {}\n'
          ' Available for Acetylation: {}\n'
          ' Acetylated: {}\n'
          ' Degree of Acetylation (DA): {:.2f} %\n'
          ' Available for Protonation: {}\n'
          ' Protonated: {}\n'
          ' Degree of Protonation (DP): {:.2f} %\n'.format(
        allres, cht + chtp + chtr, chtr, acil_d * 100.,
                cht0 + cht + chtp, chtp, prot_degree * 100))
    return chtr / allres, chtp / allres


def calc_coord(coordinations, radius=0.102103, angle=139., angle2=79.):
    """
    Calculating xyz-coordinates from:
    :param coordinates: 3 atom in chitosan
    :param radius: distance from first atom and target atom
    :param angle: angle target at---first at---secondat
    :param angle2: dihedral angle
    """
    na = 0
    nb = 1
    nc = 2
    X = [0 * x for x in range(3)]
    Y = [0 * x for x in range(3)]
    Z = [0 * x for x in range(3)]
    for i in range(3):
        X[i] = coordinations[nc][i] - coordinations[nb][i]
        Y[i] = coordinations[na][i] - coordinations[nb][i]
    Z[0] = X[1] * Y[2] - X[2] * Y[1]
    Z[1] = X[2] * Y[0] - X[0] * Y[2]
    Z[2] = X[0] * Y[1] - X[1] * Y[0]
    X[0] = Y[1] * Z[2] - Y[2] * Z[1]
    X[1] = Y[2] * Z[0] - Y[0] * Z[2]
    X[2] = Y[0] * Z[1] - Y[1] * Z[0]

    temp = math.pow(X[0] * X[0] + X[1] * X[1] + X[2] * X[2], 0.5)
    help = math.pow(Y[0] * Y[0] + Y[1] * Y[1] + Y[2] * Y[2], 0.5)
    work = math.pow(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2], 0.5)

    for i in range(3):
        X[i] = X[i] / temp
        Y[i] = Y[i] / help
        Z[i] = Z[i] / work

    R = radius * math.sin(angle * 0.0174533)
    D = - radius * math.cos(angle * 0.0174533)
    E = R * math.cos(angle2 * 0.0174533)
    H = R * math.sin(angle2 * 0.0174533)

    vector = [0 for x in range(3)]
    result = [0 for x in range(3)]
    for i in range(3):
        vector[i] = E * X[i] + D * Y[i] + H * Z[i]
        result[i] = coordinations[na][i] + vector[i]
    return result


class ReadGro:
    def __init__(self, input_file):
        self._allAtoms = []
        self.input_file = input_file

    def read(self):
        """
        Reading gro file and cutting water and ions.
        if your atom has many empty strings or strings with many spaces
        can be a problem.
        """
        count = 0
        for line in open(self.input_file):
            count += 1
            try:
                if count > 2 and len(line) > 43 and line[5:10].strip() != 'SOL' \
                        and line[5:10].strip() != 'HOH' and line[5:10].strip() != 'CL':
                    self.all_atoms.append([2 * int(line[0:5]), line[5:10].strip(), line[10:15].strip(),
                                           int(line[15:20]), float(line[20:28]), float(line[28:36]),
                                           float(line[36:44])])
            except ValueError:
                print('Check your output file, because one error occurred while program was reading input file.\n'
                      'May be there are many spaces in the end of file. Or string length about box size is longer '
                      'than 44')
            if len(line) > 3:
                self._boxsize = line
        self._deacil_all()
        self._deprotonated_all()

    def _deacil_all(self):
        """
        Go for system, cut ACE2 and rename residue
        """
        for atom in range(len(self._allAtoms) - 1, 0, -1):
            if self._allAtoms[atom][1] == 'CHTR':
                self._allAtoms[atom][1] = 'CHT'
            if self._allAtoms[atom][1] == 'ACE2' and self._allAtoms[atom][2] == 'C':
                self._allAtoms[atom][1] = 'CHT'
                self._allAtoms[atom][2] = 'H22'
                coord = [self._allAtoms[atom - 3][4:7], self._allAtoms[atom - 2][4:7],
                         self._allAtoms[atom - 10][4:7]]
                self._allAtoms[atom][4:7] = calc_coord(coord, 0.102, 109, 159)
                self._allAtoms[atom][0] = self._allAtoms[atom - 2][0]
            if self._allAtoms[atom][1] == 'ACE2':
                del self._allAtoms[atom]

    def _deprotonated_all(self):
        for atom in range(len(self._allAtoms) - 1, 0, -1):
            if self._allAtoms[atom][1] == 'CHTP':
                self._allAtoms[atom][1] = 'CHT'
                if self._allAtoms[atom][2] == 'H23':
                    del self._allAtoms[atom]

    @property
    def all_atoms(self):
        return self._allAtoms

    @property
    def box_size(self):
        return self._boxsize


class WriteGro:
    def __init__(self, output_file, all_atom, box_size, degrees):
        self.outputFile = output_file
        self.allAtom = all_atom
        self.boxSize = box_size
        self.degrees = degrees

    def write(self):
        file = open(self.outputFile, 'w')
        file.write('!Acilation degree = {:.3f}; Protonation degree = {:.3f}\n'.format(*self.degrees))
        file.write(str(len(self.allAtom)) + '\n')
        count = 1
        tps = set()
        for i in self.allAtom:
            tps.add(i[0])
        tps = list(tps)
        tps.sort()
        for atom in range(len(self.allAtom) - 1, -1, -1):
            if self.allAtom[atom][0] in tps:
                self.allAtom[atom][0] = tps.index(self.allAtom[atom][0]) + 1

        for atom in self.allAtom:
            line = '%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % \
                   (atom[0], atom[1], atom[2], count, atom[4], atom[5], atom[6])
            file.write(line)
            if count == 99999:
                count = 0
            count += 1
        file.write(self.boxSize)  # box size from input system
        file.close()


class AddProton:
    def __init__(self, all_atoms, pd):
        self._allAtoms = all_atoms
        self._protonationDegree = pd

    def calc(self):  # main method
        self.choose_cht()
        x = len(self._allCht)
        while len(self._allCht) > round(self._protonationDegree * x):
            del self._allCht[random.randint(0, len(self._allCht) - 1)]
        for residue in self._allCht:
            self._protonated(residue)

    def choose_cht(self):
        self._allCht = set()
        for atom in self._allAtoms:
            if atom[1] == 'CHT':
                self._allCht.add(atom[0])

        self._allCht = list(self._allCht)

    def _protonated(self, res):
        residue = []  # list with one residue
        coord = [0, 1, 2]
        for i in range(len(self._allAtoms)):  # rename residue and append it to list 'residue'
            if self._allAtoms[i][0] == res and self._allAtoms[i][1] == 'CHT':
                self._allAtoms[i][1] = 'CHTP'
                residue.append(self._allAtoms[i])

        count = self._allAtoms.index(residue[0])

        for atom in residue:  # add 3 atom for calc xyz coords for H23
            if atom[2] == 'N2':
                coord[0] = (atom[4], atom[5], atom[6])
            elif atom[2] == 'H21':
                coord[1] = (atom[4], atom[5], atom[6])
            elif atom[2] == 'H22':
                coord[2] = (atom[4], atom[5], atom[6])

        h23 = [res, 'CHTP', 'H23', 0]
        h23 += calc_coord(coord)
        residue.insert(8, h23)
        self._allAtoms[count:count + 15] = residue  # rewrite cht to chtp

    @property
    def protonated_system(self):
        return self._allAtoms


class AddAcil:

    def __init__(self, all_atoms, ad):
        self._allAtoms = all_atoms
        self._acilDegree = ad

    def calc(self):  # main method
        self.choose_cht()
        x = len(self._allCht)
        while len(self._allCht) > round(self._acilDegree * x):
            del self._allCht[random.randint(0, len(self._allCht) - 1)]
        for residue in self._allCht:
            self._acil_res(residue)

    def choose_cht(self):
        self._allCht = set()
        for atom in self._allAtoms:  # choose free cht for acilation
            if atom[1] == 'CHT':
                self._allCht.add(atom[0])
        self._allCht = list(self._allCht)

    def _acil_res(self, res):
        meet = False
        N2 = 0
        H21 = 0
        O3 = 0
        H22 = 0
        note = 0
        residue = []
        for atom in range(len(self._allAtoms)):
            if self._allAtoms[atom][0] == res and self._allAtoms[atom][1] == 'CHT':
                note = self._allAtoms[atom]
                self._allAtoms[atom][1] = 'CHTR'  # rename residue
                if self._allAtoms[atom][2] == 'H22':
                    H22 = atom
                elif self._allAtoms[atom][2] == 'H21':
                    H21 = atom
                elif self._allAtoms[atom][2] == 'N2':
                    N2 = atom
                elif self._allAtoms[atom][2] == 'C1':
                    O3 = atom
                residue.append(self._allAtoms[atom])
        count = self._allAtoms.index(note) - 14
        del self._allAtoms[H22]
        coord = [residue[N2 - count][4:7], residue[H21 - count][4:7],
                 residue[O3 - count][4:7]]
        ace = [[res + 1, 'ACE2', 'C', 0]]
        ace[0] += calc_coord(coord, 0.142, 112.191, 161.107)

        ace.append([res + 1, 'ACE2', 'O', 0])
        ace[1] += calc_coord(coord, 0.231, 137.728, 172.218)

        ace.append([res + 1, 'ACE2', 'CA', 0])
        ace[2] += calc_coord(coord, 0.249, 80.828, 150.604)
        self._allAtoms.insert(count + 14, ace[2])
        self._allAtoms.insert(count + 14, ace[1])
        self._allAtoms.insert(count + 14, ace[0])

    @property
    def acil_system(self):
        return self._allAtoms


if __name__ == '__main__':
    # input_file = sys.argv[sys.argv.index('-f')+1]
    # outputFile = sys.argv[sys.argv.index('-o')+1]
    # protonation_degree = float(sys.argv[sys.argv.index('-dp')+1])/100.
    # acilation_degree = float(sys.argv[sys.argv.index('-da')+1])/100.

    input_file = 'files\\cht.gro'  # 'chitin4-ace.gro'

    acilation_degree = 0.4
    protonation_degree = 0.4
    outputFile = 'files\\cht_{}_ad_{}_pd.gro'.format(acilation_degree, protonation_degree)
    readFile = ReadGro(input_file)
    readFile.read()
    degrees = calc_degrees(readFile.all_atoms, acilation_degree, protonation_degree)

    acilation = AddAcil(readFile.all_atoms, degrees[0])
    acilation.calc()

    protonation = AddProton(readFile.all_atoms, degrees[1])
    protonation.calc()

    calc_degrees(acilation.acil_system)

    w = WriteGro(outputFile, readFile.all_atoms, readFile.box_size, print_info(readFile.all_atoms))
    w.write()
