"""
        self.o

"""

from bead_spring_system import make_system
from boltzmann_invert   import PairTable
import lammps

"""
  
"""
class InverseBoltzmannIterator:
    #
    def __init__(self, opt):
        self.iteration_ct = 0
        self.options      = opt
        self.pair_table   = PairTable(opt.mdtemp, 'rdf', False)
        self.lmp_data     = 'coarse_system.lammps'
        self.bond_r0 = 5.0 # TODO fix me
        self.bond_k  = 5.0        

        make_system(self.lmp_data, opt.nchains, opt.blockstr,
                    opt.nblocks, self.bond_r0)

    # Makes an estimate of the coarse-grain potentials by Boltzmann inversion.
    def iterate(self):
        dump_file = 'coarse_samples.lammpstrj'

        param = {'k': self.bond_k, 'r0':self.bond_r0, 'T':self.options.mdtemp,
                 'data':self.lmp_data, 'dump':dump_file, 'iteration':self.iteration_ct}
        
        self.pair_table.write_lammps('pair.table.%d' % self.iteration_ct,
                                     'SS', self.iteration_ct)
        lammps.run_coarse(param)
        self.iteration_ct += 1
         

"""
NOT USED ANYMORE!!!
def main():

    # Make the initial coarse system (this should only need to be one once.
    r_range = (min(r), max(r), 0.1)
    b_range = (0.0,    15.0,   0.1)
    a_range = (0.0,    180.0,  1.0)
    for i in range(5):
        distribution.compute_rdf(read_data, dump, 'coarse', r_range, b_range, a_range)
        distribution.compare('rdf', '.', 'rdf-comparison-it%d.png'%i)
        
        r,e,f = force_table.corrected_pair_table(temp, f)
        lammps.write_table('pair.table.%d'%(i+1), 'SS', r,e,f)
        param['iteration'] = i+1
        lammps.run_coarse(param)
"""     
