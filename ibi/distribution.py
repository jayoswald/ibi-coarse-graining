#!/usr/bin/env python
from numpy import * 
from scipy import optimize
from glob import glob
import os
import paths
from matplotlib import pyplot as py

"""
    This module provides functions that are used to calculate the 
    radial density, bond length, and angle distribution functions.
"""

# Computes the rdf/bdf/adf functions by calling cg-post.
#  lmp_data is the path to the LAMMPS read_data input file.
#  out is what we want cg-post to write to e.g. rdf_out_01.txt
#  type are the bead types to compute (expecting a list)
def compute_rdf(lmp_data, types, dump, out, r_range, b_range, a_range):
    cgini  = 'input    %s\n' % lmp_data
    cgini += 'dump     %s\n' % dump
    cgini += 'output   %s\n' % out
    cgini += 'rdf      %f %f %f\n' % tuple(r_range)
    cgini += 'bdf      %f %f %f\n' % tuple(b_range)
    cgini += 'adf      %f %f %f\n' % tuple(a_range)

    for t in types: 
        cgini += 'type %d %s\n' %(types[t], t)

    for t in types:
        cgini += '\nbead %s %s' %(t,t)
        cgini += '\ncenter %s 0' %t
    
    f = open('cg.ini', 'w')
    f.write(cgini)
    f.close()
    os.system(paths.cgpost)

# Gets specific iteration rdf files.
def iteration_rdf_files(it, pair):
    pair_swap = pair[1] + pair[0]
    rdf_files  = glob('rdf-cg-%02d_%s.txt' % (it, pair))
    rdf_files += glob('rdf-cg-%02d_%s.txt' % (it, pair_swap))
    return rdf_files

# Returns all rdf files in the rdf folder..
def md_rdf_files():
    return glob('rdf/rdf*.txt')

# Reads the rdf/bdf/adf files and computes the average.
def average_rdf(rdf_files):
    if len(rdf_files) == 0:
        print 'No rdf files found'
        import sys
        sys.exit(1)
    rdf = [genfromtxt(f) for f in rdf_files]
    for i in range(1,len(rdf)):
        rdf[0][:,1] += rdf[i][:,1]
    rdf[0][:,1] *= 1.0/float(len(rdf))
    return rdf[0]

# Compares an iteration to the md rdf.
def compare(it, pair, figpath=''):
    md_rdf = average_rdf(md_rdf_files())
    cg_rdf = average_rdf(iteration_rdf_files(it, pair))

    py.clf()
    py.plot(md_rdf[:,0], md_rdf[:,1])
    py.plot(cg_rdf[:,0], cg_rdf[:,1])
    py.xlabel('Radial distance ($\AA$)')
    py.ylabel('Radial density function g(r)')
    py.legend(['All atom model', 'coarse grain model'], loc='lower right')

    if figpath == '': py.show()
    else:             py.savefig(figpath)


