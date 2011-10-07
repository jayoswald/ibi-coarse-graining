#!/usr/bin/env python
from numpy import * 
from scipy import optimize
import glob
import os
import paths

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


# Reads the rdf/bdf/adf files and computes the average.
def average_rdf(t='r', folder='rdf'):
    s = os.path.join(folder, '%sdf' %t)
    df_files = glob.glob(s+'*.txt')
    df = None
    for f in df_files:
        data = open(f, 'r').readlines()
        data = array([[float(x) for x in s.split()] for s in data])
        if df==None: df = data
        else:        df[:,1] += data[:,1]
    df[:,1] /= len(df_files)
    return df


from matplotlib  import pyplot as py

""" 
    Compares two distributions.
"""
def compare(path1, path2, path=''):
    md_rdf = average_rdf('r')
    cg_rdf = average_rdf('r', '')

    md_rdf[:,1] /= md_rdf[-1,1]
    py.clf()
    py.plot(md_rdf[:,0], md_rdf[:,1])
    py.plot(cg_rdf[:,0], cg_rdf[:,1])
    py.xlabel('Radial distance ($\AA$)')
    py.ylabel('Radial density function g(r)')
    py.legend(['All atom model', 'coarse grain model'], loc='lower right')

    if path == '': py.show()
    else:          py.savefig(path)


# For now, debugging will call the compare function.
if __name__ == '__main__':
    compare()


