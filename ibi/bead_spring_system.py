#!/usr/bin/env python
from numpy import * 
import glob

"""
  Makes chains of atoms.
  Parameters:
      Block:   should be like SSSSSHH
      nblk:    number of blocks in a chain.
      nchain:  number of chains in the model.
"""
def make_system(path, nchain, block=14*'S', nblk=1, r0=5.0, rvdw=5.0):
    print 'Making', nchain, 'chains of', nblk*block
    chain_size = len(block)*nblk
    natoms     = chain_size*nchain
    nbonds     = (chain_size-1)*nchain

    # Counts the number of unique species in block.
    atomtypes = list(set(block))
    print 'Bead types are: {', ' '.join(atomtypes),'}'
    natomtypes = len(atomtypes)
    # Number of bonds is 1 + 2 + ... + natomtypes.
    nbondtypes = sum([i+1 for i in range(natomtypes)]) 

    f = open(path, 'w')
    f.write('LAMMPS 2005 data file for soft beads\n\n')
    f.write('% 7d atoms\n'        %natoms)
    f.write('% 7d bonds\n\n'      %nbonds)
    f.write('% 4d atom types\n'   %natomtypes)
    f.write('% 4d bond types\n\n' %nbondtypes)

    nx = ceil(sqrt(nchain))
    ny = ceil(float(nchain)/nx)
    nz = chain_size 

    f.write(' %15.9f %15.9f xlo xhi\n'%   (0.0, nx*rvdw))
    f.write(' %15.9f %15.9f ylo yhi\n'%   (0.0, ny*rvdw))
    f.write(' %15.9f %15.9f zlo zhi\n\n'% (0.0, nz*r0))

    f.write('Masses\n\n')
    f.write('1 72.0\n\n') # TODO: not exact


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
        x *= rvdw
        y *= rvdw
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
# chain_maker.py <nchain=20> <block=SSSSSSS> <nblock=1> <r0=5.0>
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
        if len(args) > 4: print 'Warning too many arguments.'
        make_chains(n, blk, nblk, r0)


# If called at top level.
if __name__ == '__main__':
    import sys
    main(sys.argv[1::])

