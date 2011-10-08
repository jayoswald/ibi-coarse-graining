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
def compute_rdf(lmp_data, dump, out, r_range, b_range, a_range):
    types = ['S']
    cgini  = 'input    %s\n' % lmp_data
    cgini += 'dump     %s\n' % dump
    cgini += 'output   %s\n' % out
    cgini += 'rdf      %f %f %f\n' % tuple(r_range)
    cgini += 'bdf      %f %f %f\n' % tuple(b_range)
    cgini += 'adf      %f %f %f\n' % tuple(a_range)

    for i in range(len(types)):
        cgini += 'type %d %s\n' %(i+1, types[i])
    cgini += 'bead S\ncenter 0'
    
    f = open('cg.ini', 'w')
    f.write(cgini)
    f.close()
    os.system(paths.cgpost)

# Gets specific iteration rdf files.
def iteration_rdf_files(it):  return glob('rdf-cg-%02d*.txt' % it)
def md_rdf_files():           return glob('rdf/rdf*.txt')


# Reads the rdf/bdf/adf files and computes the average.
def average_rdf(rdf_files):
    rdf = [genfromtxt(f) for f in rdf_files]
    for i in range(1,len(rdf)):
        rdf[0][:,1] += rdf[i][:,1]
    rdf[0][:,1] *= 1.0/float(len(rdf))
    return rdf[0]

# Compares an iteration to the md rdf.
def compare(it, figpath=''):
    md_rdf = average_rdf(md_rdf_files())
    cg_rdf = average_rdf(iteration_rdf_files(it))

    py.clf()
    py.plot(md_rdf[:,0], md_rdf[:,1])
    py.plot(cg_rdf[:,0], cg_rdf[:,1])
    py.xlabel('Radial distance ($\AA$)')
    py.ylabel('Radial density function g(r)')
    py.legend(['All atom model', 'coarse grain model'], loc='lower right')

    if figpath == '': py.show()
    else:             py.savefig(figpath)

