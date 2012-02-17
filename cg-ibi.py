#!/usr/bin/python
"""
   cg-ibi.py
"""
import ibi, ibi.ibi
import sys
from ibi.read_ibi_config import IbiConfig
#
def main():
    
    inifile = 'ibi.ini'
    if len(sys.argv) > 1:  inifile = sys.argv[1]
    options = IbiConfig(inifile)

    iterator = ibi.ibi.InverseBoltzmannIterator(options)
    iterator.iterate(options.iterations)

if __name__ == '__main__': main() 

