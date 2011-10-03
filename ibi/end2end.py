#!/usr/bin/env python

from numpy import array, zeros, dot
from numpy.linalg import norm

# Reads the dump file and 
def compute_end_to_end(dump, ifile):
    f = open(ifile, 'r')

    # Skip to atoms.
    line = f.readline()
    while len(line)>0 and not line.startswith('Atoms'): 
        line = f.readline()

    # Read all atoms in input. 
    molecule = {}
    while True: 
        line = f.readline()
        if len(line)==0 or line.startswith('Bonds'): break
        cols = line.split()
        if len(cols) == 0: continue
        tag = int(cols[0])
        m   = cols[1]
        xyz = array([float(x) for x in cols[-3::]])
        if not m in molecule: molecule[m] = []
        molecule[m].append([tag, xyz])
    f.close()
    # Read output file.    
    
    f = open(dump, 'r')
    
    step = []
    rr   = []
    while True:
        line = f.readline()
        if line=='': break
        if line.startswith('ITEM: TIMESTEP'):
            step.append(int(f.readline()))
        elif line.startswith('ITEM: NUMBER OF ATOMS'):
            num_atoms = int(f.readline())
        elif line.startswith('ITEM: BOX BOUNDS'):
            bounds = [[float(x) for x in f.readline().split()] for i in [0,1,2]]
        elif line.startswith('ITEM: ATOMS'):
            var = line.split()[2::]
            mi  = var.index('mol')
            xi  = var.index('xu')
            yi  = var.index('yu')
            zi  = var.index('zu')
            xyz = zeros((num_atoms, 3))
            for i in range(num_atoms):
                cols = f.readline().split()
                tag  = int(cols[0])
                xyz[tag-1,:] = [float(cols[j]) for j in [xi,yi,zi]]
            r0t = 0.0
            for m in molecule:                
                beg = molecule[m][0][0]
                end = molecule[m][-1][0]
                r0  = molecule[m][-1][1]-molecule[m][0][1]
                rt  = xyz[end-1,:] - xyz[beg-1,:]

                r0t += dot(r0,rt)/norm(r0)/norm(rt)
            rr.append(r0t/len(molecule))
    from matplotlib import pyplot
    pyplot.plot(step, rr)
    pyplot.show()






import sys

dump  = 'dump-equil.lammpstrj'
ifile = 'soft.lammps'

    

compute_end_to_end(dump, ifile)
