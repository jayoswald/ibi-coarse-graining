#!/usr/bin/env python
from operator import mul
from scipy import constants

class LmpData:
    def __init__(self, path):
        self.path = path
        fid = open(self.path, 'r')
        self.box = [[]]*3
        while True:
            line = fid.readline()
            if line == '': break
            line = line.split()
            if len(line) == 0: continue

            # Read number of atoms and interactions.
            elif line[-1] == 'atoms':     self.natoms = int(line[0])
            elif line[-1] == 'bonds':     self.nbonds = int(line[0])
            elif line[-1] == 'angles':    self.nangle = int(line[0])
            elif line[-1] == 'dihedrals': self.ndihed = int(line[0])
            elif line[-1] == 'impropers': self.nimpro = int(line[0])

            # Read number of types.
            if line[-1] == 'types':
                if line[-2] == 'atom':     self.natomtypes  = int(line[0])
                if line[-2] == 'bond':     self.nbondtypes  = int(line[0])
                if line[-2] == 'angle':    self.nangletypes = int(line[0])
                if line[-2] == 'dihedral': self.ndihedtypes = int(line[0])
                if line[-2] == 'improper': self.nimprotypes = int(line[0])

            # Get box size.
            if   line[-1] == 'xhi' and line[-2] == 'xlo':
                self.box[0] = [float(line[0]), float(line[1])]
            elif line[-1] == 'yhi' and line[-2] == 'ylo':
                self.box[1] = [float(line[0]), float(line[1])]
            elif line[-1] == 'zhi' and line[-2] == 'zlo':
                self.box[2] = [float(line[0]), float(line[1])]
    
            if line[0] == 'Masses':
                self.mass = {}
                while len(self.mass) < self.natomtypes:
                    line = fid.readline()
                    if line == '': 
                        print 'Unexpected end-of-file in Masses section'
                        break
                    line = line.split()
                    if len(line) == 2:
                        self.mass[int(line[0])] = float(line[1])
            # Read atoms.
            if line[0] == 'Atoms':
                ct = 0
                self.system_mass = 0.0
                while ct < self.natoms:
                    line = fid.readline()
                    if line == '':
                        print 'Unexpected end-of-file in Atoms section'
                        break
                    line = line.split()
                    if len(line) > 4:
                        ct += 1
                        atom_type = int(line[2])
                        self.system_mass += self.mass[atom_type]
    # Returns the density of the system.
    def density(self):
        v = reduce(mul, [b[1] - b[0] for b in self.box])
        rho = self.system_mass / v # in g/mol / A^3
        rho /= 1e-24*constants.N_A # converts to g/cc
        return rho

if __name__ == '__main__':                                    
    import sys
    lmp = LmpData(sys.argv[1])
    print 'density is', lmp.density(), 'g/cc'
    print 'mass is ', lmp.system_mass, 'g/mol'
