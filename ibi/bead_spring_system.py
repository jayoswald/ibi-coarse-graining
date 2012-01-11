#!/usr/bin/env python
from numpy import * 
from scipy import constants
import glob

"""
  Class to hold the bead positions and connectivies.
"""
class BeadSpringSystem:

    # Initialize an empty system.
    def __init__(self):
        self.beads      = []
        self.bonds      = []
        self.box_length = 1.0
        self.bead_masses = []

    # Writes an input file to LAMMPS following format at
    # http://lammps.sandia.gov/doc/read_data.html
    def write_to_lammps(self, lmpdatafile):
        f = open(lmpdatafile, 'w')

        # Get number of unique bead and bond types.
        num_bead_types = len(set([b[1] for b in self.beads]))
        num_bond_types = len(set([b[0] for b in self.bonds]))

        # Write header.
        f.write('LAMMPS 2005 data file for soft beads\n\n')
        f.write('% 7d atoms\n'        %len(self.beads))
        f.write('% 7d bonds\n\n'      %len(self.bonds))
        f.write('% 4d atom types\n'   %num_bead_types)
        f.write('% 4d bond types\n\n' %num_bond_types)
        # Write box size.
        f.write(' %15.9f %15.9f xlo xhi\n'%   (0.0, self.box_length))
        f.write(' %15.9f %15.9f ylo yhi\n'%   (0.0, self.box_length))
        f.write(' %15.9f %15.9f zlo zhi\n\n'% (0.0, self.box_length))

        f.write('Masses\n\n')
        for m in self.bead_masses:
            f.write('1 %f\n'%m)
        f.write('\n')

        f.write('Atoms\n\n')
        for index,bead in enumerate(self.beads):
            bead = tuple([index+1] + bead[0:2] + list(bead[2]))
            f.write(' % 6d % 6d % 3d %.3f %.3f %.3f\n' %bead)
        f.write('\n')

        f.write('Bonds\n\n')
        for index,bond in enumerate(self.bonds):
            bond = tuple([index+1] + list(bond))
            f.write('% 6d % 3d % 6d % 6d\n' %bond)
        f.close()

"""
  Creates a bead spring system. 
    Input parameters:
        block:       definition of a chain subblock (e.g. 'AAB').
        nchain:      the number of chains in the system.
        nblk:        the number of blocks in a chain.
        bond_length: the distance between bonded atoms.
        box_length:  the size of the sides of the cubic system box.
"""
def create(block, nchain, nblk, bond_length, box_length):

    # Define atom_types such that for each bead there is a unique number.
    atom_types = {}
    for i in block:
        if i not in atom_types:
            atom_types[i] = len(atom_types) + 1

    # Gets bond types for each combination of pairs of beads.
    bond_types = {}
    for i in range(len(block)-1):
        pair = min(block[i:i+2]) +  max(block[i:i+2])
        if pair not in bond_types:
            bond_types[pair] = len(bond_types)+1
    # If there are more than one blocks then we need to check last, first.
    if nblk > 1:
        pair = ''.join(sorted([block[0], block[-1]]))
        if pair not in bond_types:
            bond_types[pair] = len(bond_types)+1

    ## FOR DEBUGGING, REMOVE LATER
    print 'Atom types are:', atom_types
    print 'Bond types are:', bond_types 
    ##

    # Given the index of an atom in a chain, return the type.
    def get_atom_type(i):
        atom = block[i%len(block)]
        return atom_types[atom]

    # Given the index, i of the first atom in a chain, 
    # return the bond type between i and i+1.
    def get_bond_type(i):
        itype = block[(i)  %len(block)]
        jtype = block[(i+1)%len(block)]
        pair = min(itype,jtype) + max(itype,jtype)
        return bond_types[pair]

    system = BeadSpringSystem()
    system.box_length = box_length
    system.atom_types = atom_types
    chain_size = len(block)*nblk

    def random_unit():
        r = random.rand(3) - 0.5
        r *= 1.0/linalg.norm(r)
        return r
        

    for molecule in range(nchain):

        first_bead = system.box_length*random.rand(3)
        system.beads.append([molecule, get_atom_type(0), first_bead])

        for z in range(chain_size - 1):
            ran_vect  = bond_length * random_unit()
            next_bead = system.beads[-1][2] + ran_vect

            # Check if angle ijk is not near zero degrees.
            if z > 1:
                while True:
                    rij = system.beads[-2][2] - system.beads[-1][2]
                    rjk = next_bead - system.beads[-1][2]
                    if dot(rij,rjk) < 0.9*bond_length**2: break
                    ran_vect  = bond_length * random_unit()
                    next_bead = system.beads[-1][2] + ran_vect

            btype = get_bond_type(z)
            system.bonds.append([btype, len(system.beads), len(system.beads)+1])

            # Add the bead to the system.
            atype = get_atom_type(z+1)
            system.beads.append([molecule, get_atom_type(z+1), next_bead])
    return system


"""
  Makes chains of atoms in a box.
  Each chain is lined up along the x-direction (i.e. not random).

  Parameters:
      path:    where to write the lammps input file.
      nchain:  number of chains in the model.
      block:   one character per bead, example for soft and hard beads: SSSSHH
      nblk:    number of blocks in a chain.
      r0:      equilibrium bond length between atoms (doesn't have to be exact)
      density: target density of the system (in g/cc, used to compute box size).
"""
# TODO: need to support more than one type of bead.
#   Each bead needs to have its own mass.
#   Also need to define r0 for each bead combination (A-A, A-B, B-B).
def make_system(path, nchain, block, nblk, bond_length, density):
    
    # Beads A,B
    # A = ....
    # B = ....

    mass_per_bead = 72.10776  # in kg/mol (soft bead)
    #mass_per_bead = 506.518  # in kg/mol (hard bead)

    #TODO: won't be true when blocks aren't the same bead.
    mass_per_block = len(block)*mass_per_bead
    mass   = mass_per_block*nblk*nchain

    print 'density is', density, 'g/cc'
    density *= 1e-24*constants.N_A 
    print 'density is', density, 'g/mol / A^3'
    volume = mass / density 
    print 'volume is', volume, 'A^3'

    box_length = volume**(1.0/3.0)

    system = create(block, nchain, nblk, bond_length, box_length)

    # TODO - needs to be fixed.
    system.bead_masses = [mass_per_bead]
    system.write_to_lammps(path)

# Standalone mode (called from command window).
# chain_maker.py <nchain=20> <block=SSSSSSS> <nblock=1> <r0=5.0> <density=1.0>
def main(args):
    if args==[]:
        print 'Using default values:'
        print '    nchain=20, block=14*S, nblock=1, r0=5.0, rho=1.0.'
        make_system('beadsystem.lammps', 20, 14*'S', 1, 5.0, 1.0)
    elif args==['-h']:
        help_msg  = 'Usage:\n\t'
        help_msg += 'bead_spring_system <nchain=20> <block=SSSSSSS> '
        help_msg += '<nblock=1> <r0=5.0> <rho=1.0>'
        print help_msg
        return None
    else:
        n,nblk = 20,1
        blk    = 8*'S'
        r0     = 0.5
        rho    = 1.0
        if len(args) > 0: n    = int(args[0])
        if len(args) > 1: blk  = args[1]
        if len(args) > 2: nblk = int(args[2])
        if len(args) > 3: r0   = float(args[3])
        if len(args) > 4: rho  = float(args[4])
        if len(args) > 5: print 'Warning too many arguments.'
        make_system(n, blk, nblk, r0, rho)


# If called at top level.
if __name__ == '__main__':
    import sys
    main(sys.argv[1::])

