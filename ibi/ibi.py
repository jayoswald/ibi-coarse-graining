"""

"""

from bead_spring_system import make_system
from boltzmann_invert   import PairTable

"""
  
"""
class InverseBoltzmannIterator:
    #
    def __init__(self, opt):
        self.iteration_ct = 0
        self.options      = opt
        self.pair_table   = PairTable(opt.mdtemp, 'rdf', True)
        self.lmp_data     = 'coarse_system.lammps'

        r0 = 5.0 # TODO fix me        
        make_system(self.lmp_data, opt.nchains, opt.blockstr, opt.nblocks, r0)
        

    # Makes an estimate of the coarse-grain potentials by Boltzmann inversion.
    def iterate(self):
        self.iteration_ct += 1

                    

# Main program.
def main():
    # TODO: get these settings from elsewhere.
    dump         = 'coarse_samples.lammpstrj'
    param        = {'k':bond_k, 'r0':bond_r0, 'T':temp, 
                    'data':read_data, 'dump':dump, 'iteration':0}

    # Make the initial coarse system (this should only need to be one once.
    r,e,f = force_table.pair_table(temp)
    lammps.write_table('pair.table.0','SS',r,e,f)

    lammps.run_coarse(param)
    
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
        
