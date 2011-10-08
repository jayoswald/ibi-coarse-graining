#!/usr/bin/python
"""
   cg-ibi.py
"""
import ibi, ibi.ibi
from optparse import OptionParser


usage = """ %prog [options] """

help = {
    'sample': 'test entry'
}


#
def main():
    parser = OptionParser(usage=usage, version='%prog ' + ibi.__version__)
    parser.add_option('-n', '--num-cpu', 
                      type='int',   dest='nproc' ,  default=4)
    parser.add_option('',   '--md-temperature', 
                      type='float', dest='mdtemp',  default=300.0)
    parser.add_option('',   '--num_chains',        
                      type='int',   dest='nchains', default=40)
    parser.add_option('',   '--num-blocks',
                      type='int',   dest='nblocks', default=14)
    parser.add_option('',   '--block',
                      dest='blockstr', default='S')

    options, args = parser.parse_args()

    iterator = ibi.ibi.InverseBoltzmannIterator(options)

    iterator.iterate(5)

if __name__ == '__main__': main() 

