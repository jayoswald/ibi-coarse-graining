#!/usr/bin/env python
from numpy import * 
from scipy import constants
import glob

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
def make_system(path, nchain, block, nblk, r0, density):
    print 'Making', nchain, 'chains of', nblk*block
    chain_size = len(block)*nblk
    nbeads     = chain_size*nchain
    nbonds     = (chain_size-1)*nchain

    # Counts the number of unique species in block.
    atomtypes = list(set(block))
    print 'Bead types are: {', ' '.join(atomtypes),'}'
    natomtypes = len(atomtypes)
    # Number of bonds is 1 + 2 + ... + natomtypes.
    nbondtypes = sum([i+1 for i in range(natomtypes)]) 

    f = open(path, 'w')
    f.write('LAMMPS 2005 data file for soft beads\n\n')
    f.write('% 7d atoms\n'        %nbeads)
    f.write('% 7d bonds\n\n'      %nbonds)
    f.write('% 4d atom types\n'   %natomtypes)
    f.write('% 4d bond types\n\n' %nbondtypes)

    nx = ceil(sqrt(nchain))
    ny = ceil(float(nchain)/nx)
    nz = chain_size 

    # Compute size of system based on total mass.
    # TODO - this is valid for soft bead only!
    mass_per_bead = 72.10776  # in kg/mol
    mass   = mass_per_bead * nbeads
    # convert density from g/cc to g/mol / A^3.
    density *= 1e-24*constants.N_A 
    print 'density is', density, 'g/mol / A^3'
    volume = mass / density 

    lz = nz * r0
    dr = sqrt(volume/lz/(nx*ny))
    lx = nx*sqrt(volume/lz/(nx*ny))
    ly = ny * lx / nx

    f.write(' %15.9f %15.9f xlo xhi\n'%   (0.0, lx))
    f.write(' %15.9f %15.9f ylo yhi\n'%   (0.0, ly))
    f.write(' %15.9f %15.9f zlo zhi\n\n'% (0.0, lz))

    f.write('Masses\n\n')
    f.write('1 %f\n\n'%mass_per_bead) 


    f.write('Atoms\n\n')
    print 'packing in a', nx, 'by', ny, 'grid'

    # Returns the atom type for position z in a chain.
    def atom_type(z):
        blockidx = z % len(block)
        return atomtypes.index(block[blockidx])+1 
    
    # Make all atoms.
    for m in range(nchain):
        x = m % nx
        y = floor(m / nx)
        x *= dr
        y *= dr
        for z in range(nz):
            tag   = z + m*nz + 1
            atype = atom_type(z)
            z *= r0
            atom = (tag, m, atype, x, y, z)
            f.write(' % 6d % 6d % 3d %.3f %.3f %.3f\n' %atom)

    f.write('\nBonds\n\n')
    ct = 0

    bond_types = []
    for m in range(nchain):
        for z in range(nz-1):
            ct  += 1

            atypes = [atom_type(z), atom_type(z+1)]
            btype = [min(atypes), max(atypes)]
            if not btype in bond_types:
                bond_types.append(btype)
            btype = bond_types.index(btype)+1

            tag   = z + m*nz + 1
            bond  = (ct, btype, tag, tag+1)
            f.write('% 6d % 3d % 6d % 6d\n' %bond)
    print 'Bond types are:', bond_types
    f.close()


# Standalone mode (called from command window).
# chain_maker.py <nchain=20> <block=SSSSSSS> <nblock=1> <r0=5.0> <density=1.0>
def main(args):
    if args==[]:
        print 'Using default values:'
        print '    nchain=20, block=14*S, nblock=1, r0=5.0, rvdw=5.0.'
        make_chains(20, 14*'S')
    elif args==['-h']:
        print 'Usage:'
        print '\tchain_maker.py <nchain=20> <block=SSSSSSS> <nblock=1> <r0=5.0>'
        return None
    else:
        n,nblk = 20,1
        blk    = 8*'S'
        r0     = 0.5
        if len(args) > 0: n    = int(args[0])
        if len(args) > 1: blk  = args[1]
        if len(args) > 2: nblk = int(args[2])
        if len(args) > 3: r0   = float(args[3])
        if len(args) > 4: rho  = float(args[4])
        if len(args) > 5: print 'Warning too many arguments.'
        make_chains(n, blk, nblk, r0, rho)


# If called at top level.
if __name__ == '__main__':
    import sys
    main(sys.argv[1::])

