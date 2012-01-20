#!/usr/bin/env python
import os, paths

# Global params
use_pbs = False

# Starts the initial lammps.
def run_coarse(param):
    ifile  = ['in.cg-ibi-1', 'in.cg-ibi-2']
    lmpin  = [in_cg1 %param, in_cg2 %param]
    modes  = ['equil',  'sample']

    for i,s,mode in zip(ifile, lmpin, modes):
	log = 'log-%s-%d.lammps' % (mode, param['iteration'])
        f   = open(i, 'w')
        f.write(s)
        f.close()
        run_lammps(i, param['nproc'], log)
    
# Runs LAMMPS via the command line.
def run_lammps(inp, nproc, log):
    lmp = paths.lammps
    if use_pbs: cmd = 'mpiexec %s -in %s -l %s' %(lmp, inp, log)
    else:       cmd = 'mpiexec -n %d %s -in %s -l %s' %(nproc, lmp, inp, log)
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
pair_coeff    * * pair.table.%(iteration)d SS
special_bonds lj/coul 0.0 1.0 1.0 
#dump          1 all custom 1 dump-equil.lammpstrj id mol xs ys zs
thermo        200
thermo_style  custom step temp press ke pe etotal
neighbor      5.0 bin
neigh_modify  every 1 delay 0

# Will move atoms to reasonable distances apart
fix           1 all nve/limit 2.0
velocity      all create 100.0 1234
fix           2 all temp/berendsen 500 500 15
timestep      10
run           2000

minimize     1e-6 1e-6 1000 1000

unfix         1
unfix         2
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

