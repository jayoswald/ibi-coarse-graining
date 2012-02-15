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
 
    # return the bond type pair 'HH,HS.....'
    def get_bond_pair(i):
        return [k for k, v in bond_types.iteritems() if v == i][0]

    system = BeadSpringSystem()
    system.box_length = box_length
    system.atom_types = atom_types
    chain_size = len(block)*nblk
    system.bond_types = bond_types

    def random_unit():
        r = random.rand(3) - 0.5
        r *= 1.0/linalg.norm(r)
        return r
        

    for molecule in range(nchain):

        first_bead = system.box_length*random.rand(3)
        system.beads.append([molecule, get_atom_type(0), first_bead])

        for z in range(chain_size - 1):

            # What bond type am I and what is the length
            # bond_type =
            btype     = get_bond_type(z)
            bpair     = get_bond_pair(btype)     
            ran_vect  = bond_length[bpair] * random_unit()
            next_bead = system.beads[-1][2] + ran_vect

            # Check if angle ijk is not near zero degrees.
            if z > 1:
                while True:
                    btypeij = get_bond_type(z-1)
                    btypejk = get_bond_type(z)
                    bpairij = get_bond_pair(btypeij)
                    bpairjk = get_bond_pair(btypejk)
                    rij     = system.beads[-2][2] - system.beads[-1][2]
                    rjk     = next_bead - system.beads[-1][2]

                    cos_angle_ijk = dot(rij,rjk) / (bond_length[bpairjk] * bond_length[bpairij])

                    if cos_angle_ijk < 0.8: break
                    ran_vect  = bond_length[bpairjk] * random_unit()
                    next_bead = system.beads[-1][2] + ran_vect

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
    
    #TODO: only works for H,S beads - probably should use a file to store this. 
    bead_mass = {'H': 253.261245, 'S':72.10776}

    mass_per_block = 0.0
    for bead in block:
        mass_per_block += bead_mass[bead]

    mass   = mass_per_block*nblk*nchain
    print '# of beads:   ', nblk*nchain*len(block)
    print 'Density is:    %0.4f g/cc' % density
    density *= 1e-24*constants.N_A 
    print 'Density is:    %0.4f g/mol/A^3' % density
    volume = mass / density 

    box_length = volume**(1.0/3.0)
    print 'Volume is:     %0.1f A^3'% volume
    print 'box length is: %0.3f A' % box_length

    # TODO: we need to pass this in from input somehow.
    random.seed(1237)
    system = create(block, nchain, nblk, bond_length, box_length)

    # Get the masses of each type of bead in the order that they appear.
    system.bead_masses = []
    bead_types = []
    for bead in block:
        if not bead in bead_types:
            bead_types.append(bead)
            system.bead_masses.append(bead_mass[bead])

    print zip(bead_types, system.bead_masses) 
    system.write_to_lammps(path)
    return system.bond_types
# Standalone mode (called from command window).
def main():


    parser = OptionParser()
    parser.add_option('', '--num_chains', type='int', dest='nchains', default=40,
                      help='sets number of chains')
    parser.add_option('', '--num-blocks', type='int', dest='nblocks', default=14,
                      help='sets number of blocks')
    parser.add_option('', '--block', dest='blockstr', default='S',
                      help='sets beads in block')
    parser.add_option('', '--filename',  dest='path', default='bead_system.lammps',
                      help='specifies output file')
    parser.add_option('','--bond_length', type='float', dest='r0', default= 5,
                      help='sets the bond length')
    parser.add_option('','--density', type='float', dest='rho', default=1.0,
                      help='sets the system density')
    opt,args = parser.parse_args()
    make_system(opt.path, opt.nchains, opt.blockstr, opt.nblocks, opt.r0, opt.rho)

# If called at top level.
if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    main()
