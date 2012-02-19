"""
        self.o

"""

from bead_spring_system import make_system
from boltzmann_invert   import PairTable
from distribution       import compute_rdf, compare
import lammps
import os

"""
  
"""
class InverseBoltzmannIterator:
    #
    def __init__(self, opt):
        self.iteration_ct = 0
        self.options      = opt
        self.pair_table   = PairTable(opt.mdtemp, 'rdf', False)
        self.lmp_data     = 'coarse_system.lammps'
        self.bond_types   = {}
        self.atom_types   = {}
        self.bond_r0      = {}
        self.bond_k       = {}

        # Reads the bond properities from the
        for p in opt.bonds:
            self.bond_r0[p] = opt.bonds[p]['r0']
            self.bond_k[p]  = opt.bonds[p]['k']
        print 'Bond lengths are:    ', self.bond_r0
        print 'Bond stiffnesss are: ', self.bond_k
        print 'Available pairs are: ', opt.pairs

        # make system and return dict of bondtype for lammps input
        self.bond_types,self.atom_types = make_system(self.lmp_data, opt.nchain, opt.block,
                    opt.nblock, self.bond_r0, opt.density,opt.masses)

    # Makes an estimate of the coarse-grain potentials by Boltzmann inversion.
    def iterate(self, max_iterations):
        for i in range(max_iterations):
            tag = 'cg-%02d' % self.iteration_ct
            dump_file = tag + '.lammpstrj'
            param = {'k': self.bond_k, 'r0':self.bond_r0,'btype': self.bond_types, 
                     'T':self.options.mdtemp,'nproc':self.options.nproc, 'data':self.lmp_data, 
                     'dump':dump_file,'iteration':self.iteration_ct,'atype':self.atom_types,
                     'pair_known':self.options.pairs,'misspair':self.options.missing_pair}
       
            self.pair_table.write_lammps('pair.table.%d' % self.iteration_ct,
                                         self.options.missing_pair, self.iteration_ct)

            # Don't rerun lammps if the output file already exisits.
            # This means you have to delete output files if you want lammps to rerun.
            if self.options.clean_run or not os.path.exists(dump_file): 
                lammps.run_coarse(param)

            # Now we need to analyze the coarse-grain rdf function.
            r = self.pair_table.distance[-1]
            r_range = (min(r), max(r), 0.1)
            b_range = (0.0,    15.0,   0.1)
            a_range = (0.0,    180.0,  1.0)

            compute_rdf(self.lmp_data, self.atom_types, dump_file, tag, r_range, b_range, a_range)
            compare(self.iteration_ct, self.options.missing_pair, 'rdf-%d.png'%i)
            self.pair_table.correction(self.iteration_ct, self.options.missing_pair)
            self.iteration_ct += 1
         
