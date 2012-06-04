#!/usr/bin/env python
import os, paths
import numpy

# Global params
use_pbs = False


# Finds the bond/pair coeffs needed for this run.
def get_forcefield_coeff(param):
     # to get the bond type for lammps input file 'HH' 'SS' 'SH'
    def get_bond_pair(i):
        return [k for k, v in btype.iteritems() if v == i][0]
    # adding bond coeff
    bond_coeff = ''
    for b in param['btype']:
        bond_coeff += 'bond_coeff    ' + str(param['btype'][b]) + ' '
        bond_coeff += str(param['k'][b]) + ' ' + str(param['r0'][b]) + '\n'
    # adding available pair coeffs
    pair_coeff = ''
    for b in param['pair_known']:
        a1 = param['atype'][b[0]]
        a2 = param['atype'][b[1]]
        pair_coeff += 'pair_coeff    ' + str(a1) + ' ' 
        pair_coeff +=  str(a2) + ' ' + param['pair_known'][b] + ' ' + b + '\n'
 
    # adding missing pair
    m1 = param['atype'][param['misspair'][0]]
    m2 = param['atype'][param['misspair'][1]]

    pair_coeff += 'pair_coeff %s %s pair.table.%d' %(m1,m2,param['iteration']) 
    pair_coeff += ' ' + param['misspair'] + '\n'
    return {'bond_coeff': bond_coeff, 'pair_coeff':pair_coeff}
   
# Initializes the coarse system.
def init_coarse(param):
    ifile = 'in.cg-ibi-init'
    fid   = open(ifile, 'w')
    # TODO param1 seems to be unused?
    param1 = param.update(get_forcefield_coeff(param))
    fid.write(in_cg_init %param)
    fid.close()

    log = 'log-init-%d.lammps' %param['iteration']
    run_lammps(ifile, param['nproc'], log)

    return log

# Runs the coarse system to determine the pressure.
def pressure_run(param):
    ifile = 'in.cg-ibi-press'
    fid   = open(ifile, 'w')
    # TODO param1 seems to be unused?
    param1 = param.update(get_forcefield_coeff(param))    
    fid.write(in_cg_press %param)
    fid.close()

    log = 'log-press-%d.lammps' %param['iteration']
    run_lammps(ifile, param['nproc'], log, quiet=1)
    return log

# Runs the coarse system to generate samples.
def sample_run(param):
    ifile = 'in.cg-ibi-sample'
    fid   = open(ifile, 'w')
    # TODO param1 seems to be unused?
    param1 = param.update(get_forcefield_coeff(param))    
    fid.write(in_cg_sample %param)
    fid.close()

    log = 'log-sample-%d.lammps' %param['iteration']
    run_lammps(ifile, param['nproc'], log, quiet=0)

    return log

# Runs LAMMPS via the command line.
def run_lammps(inp, nproc, log, quiet=0):
    lmp = paths.lammps
    if use_pbs: cmd = 'mpiexec %s -in %s -l %s' %(lmp, inp, log)
    else:       cmd = 'mpiexec -n %d %s -in %s -l %s' %(nproc, lmp, inp, log)

    if quiet: cmd += ' > /dev/null'

    os.system(cmd)

# Writes a table to a file designated by path.
def write_table(path, keyword, r, e, f):
    fid = open(path, 'w')
    fid.write(keyword+'\n')
    fid.write('N %d R %f %f\n\n' %(len(r), min(r), max(r)))
    for i in range(len(r)):
        fid.write('%d %f %f %f\n' %(i, r[i], e[i], f[i]))


# Extracts the pressure from the lammps log file.
def read_pressure_from_log(log):
    try:
        fid = open(log, 'r')
    except IOError:
        print 'Failed to open LAMMPS log file:', log
    while fid:
        line = fid.readline()
        if line == '': 
            print 'Previous LAMMPS run,', log, 'failed.'
            sys.exit(1)
        if line.startswith('Step'): 
            break
    pcol = line.split().index('Press')
    p = []
    while fid:
        line = fid.readline()
        if line == '' or line.startswith('Loop'):
            break
        p.append(float(line.split()[pcol]))
    return numpy.mean(p)

# Initial input script used to equilibriate the system.
in_cg_init = """
units         real
atom_style    bond
boundary      p p p 
read_data     %(data)s 
bond_style    harmonic
%(bond_coeff)s
pair_style    table linear 1000 
%(pair_coeff)s
special_bonds lj/coul 0.0 1.0 1.0 
#dump          1 all custom 1 dump-equil.lammpstrj id mol xs ys zs
thermo_style  custom step temp press ke pe etotal
neighbor      5.0 bin
neigh_modify  every 1 delay 0

# Will move atoms to reasonable distances apart
fix           1 all nve/limit 2.0
velocity      all create 100.0 1234
fix           2 all temp/berendsen 500 500 15
timestep      10
thermo        100
run           500

# Minimizes energy.
velocity      all scale 100
minimize      1e-6 1e-6 5000 5000

# Equilibriate at temperature.
unfix         1
unfix         2
fix           1 all nvt temp %(T)f %(T)f 200
timestep      40
thermo        500
run           20000
write_restart restart.equil
"""

# Input script used to generate samples.
in_cg_sample = """
read_restart  restart.equil
fix           1 all nvt temp %(T)f %(T)f 200
%(pair_coeff)s
compute       msd     all msd
thermo_style  custom step temp press vol c_msd[4] pe ke
dump          1 all custom 250 %(dump)s id type mol xs ys zs vx vy vz
thermo        1000
run           40000
"""

# Input script used for pressure control
in_cg_press = """
read_restart restart.equil
fix           1 all nvt temp %(T)f %(T)f 200
%(pair_coeff)s
thermo_style  custom step temp press
thermo        10 
run           30000 
write_restart restart.equil
"""
