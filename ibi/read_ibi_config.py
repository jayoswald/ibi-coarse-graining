"""
  Reads the ...
"""

class IbiConfig:
    def __init__(self, config_file):
        f = open(config_file, 'r')

        self.nchain  = 0
        self.block   = ''
        self.nblock  = 0
        self.density = 0.0
        self.nproc   = 4
        self.mdtemp  = 300.0 
        self.iterations = 5
        self.clean_run = 1

        self.bonds = {}
        self.pairs = {}
        self.masses  = {}
        try:
            for line in f:
                if line == '': break
                line = line.split('#')[0].strip().split()
                if len(line) == 0: continue

                if line[0] == 'nchain':
                    self.nchain = int(line[1])
                elif line[0] == 'block':
                    self.block = line[1]
                elif line[0] == 'nblock':
                    self.nblock = int(line[1])
                elif line[0] == 'density':
                    self.density = float(line[1])
                # bond types r0 value k value
                elif line[0] == 'bond':
                    bond = {}
                    for i in range(2,len(line),2):
                        if line[i] == 'r0':
                            bond['r0'] = float(line[i+1])
                        elif line[i] == 'k':
                            bond['k'] = float(line[i+1])
                        else: 
                            print 'Unknown keyword in', line[i]
                    self.bonds[''.join(sorted(line[1]))] = bond
                # bead masses
                elif line[0] == 'mass':
                    self.masses[''.join(sorted(line[1]))] = float(line[2])
                # pair types ptable
                elif line[0] == 'pair':
                    self.pairs[''.join(sorted(line[1]))] = line[2] 
                elif line[0] == 'iterations':
                    self.iterations = int(line[1])
                elif line[0] == 'md_temperature':
                    self.mdtemp = float(line[1])
                elif line[0] == 'num_cpu':
                    self.nproc = int(line[1])
                elif line[0] == 'clean_run':
                    self.clean_run = 1
                else:
                    print 'Unknown input in ini file', line
        except:
                print 'Syntax was invalid:', ' '.join(line)
                import sys
                sys.exit(1)

        self.check()

    # Determines if the required pair and bond data is known and 
    # if there is one unknown pair to fit.
    def check(self):
        block = self.block
        types = set(block)
        all_pairs = set()
        missing_masses = [p for p in types if p not in self.masses]
        if len(missing_masses) == 0:
            print 'All masses defined.'
        else:
            print 'mass undefined:', missing_masses

        for a in types:
            for b in types:
                all_pairs.add(min(a,b) + max(a,b))

        missing_pairs = [p for p in all_pairs if p not in self.pairs]
        if len(missing_pairs) == 1:    
            print '1 pair is not defined:', missing_pairs[0]
            self.missing_pair = missing_pairs[0]
        elif len(missing_pairs) == 0:  
            print 'No pairs are missing, nothing to do'
            self.missing_pair = None 
        else:
            print 'Too many pairs are undefined:', missing_pairs
        
        all_bonds = set()
        for b1,b2 in zip(block, block[1::]+block[0]):
            all_bonds.add(min(b1,b2) + max(b1,b2))
        missing_bonds = [b for b in all_bonds if b not in self.bonds]
        if len(missing_bonds) == 0:
            print 'All bonds defined.'
        else: 
            print 'Bonds undefined:', missing_bonds
         

if __name__ == '__main__':
    config = IbiConfig('ibi.ini')
