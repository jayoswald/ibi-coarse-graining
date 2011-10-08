#!/usr/bin/env python
import os, paths

# Global params
use_pbs = False

# Starts the initial lammps.
def run_coarse(param):
    ifile  = ['in.cg-ibi-1', 'in.cg-ibi-2']
    lmpin  = [in_cg1 %param, in_cg2 %param]

    for i,s in zip(ifile, lmpin):
        f = open(i, 'w')
        f.write(s)
        f.close()
        run_lammps(i, param['nproc'])
    
# Runs LAMMPS via the command line.
def run_lammps(inp, nproc):
    lmp = paths.lammps
    if use_pbs: cmd = 'mpiexec %s -in %s' %(lmp, inp)
    else:       cmd = 'mpiexec -n %d %s -in %s' %(nproc, lmp, inp)
    os.system(cmd)

# Writes a table to a file designated by path.
def write_table(path, keyword, r, e, f):
    fid = open(path, 'w')
    fid.write(keyword+'\n')
    fid.write('N %d R %f %f\n\n' %(len(r), min(r), max(r)))
    for i in range(len(r)):
        fid.write('%d %f %f %f\n' %(i, r[i], e[i], f[i]))


# Initial input script used to equilibriate the system.
in_cg1 = """
units         real
atom_style    bond
boundary      p p p 
read_data     %(data)s 
bond_style    harmonic
bond_coeff    1 %(k)f %(r0)f
pair_style    table linear 1000 
pair_coeff    1 1 pair.table.%(iteration)d SS
special_bonds lj/coul 0.0 1.0 1.0 
#dump          1 all custom 100 dump-equil.lammpstrj id mol xu yu zu
thermo        200
thermo_style  custom step temp press ke pe etotal

fix           1 all nvt temp %(T)f %(T)f 100
velocity      all create %(T)f 123456 
timestep      50
run           10000
write_restart restart.equil
"""

# Input script used to generate samples.
in_cg2 = """
read_restart  restart.equil
fix           1 all nvt temp %(T)f %(T)f 200
pair_coeff    1 1 pair.table.%(iteration)d SS
thermo        100
dump          1 all custom 250 %(dump)s id mol xs ys zs
run           20000
write_restart restart.samples
"""

