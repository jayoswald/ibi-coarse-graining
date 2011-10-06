#!/usr/bin/env python
import sys
import matplotlib.pyplot as py
from scipy import *
from scipy.optimize import leastsq

"""
  Main entry point fo the code - passes x,y values to the fitting function.
"""
def main():
    # If called with no arguments, print error message.
    if len(sys.argv)==1:
       print 'fitting.py file1 [file2] [...]'
    # If called with 1 argument (then compute fit for single file).
    elif len(sys.argv)==2:
        x,y = read_xy(sys.argv[1])
    # If called with multiple arguments, then average the distribution functions.
    # Note they must have the same x values.
    else:
        ysum = None
        for f in sys.argv[1::]:
            x,y = read_xy(f)
            if ysum==None: ysum = y
            else:          ysum += y
        y = ysum/(len(sys.argv)-1)
    fit_distribution(x,y,sys.argv[1][0])

# Returns two column data.
def read_xy(filename):
    lines = open(filename, 'r').readlines()
    xy = array([[float(s) for s in line.split()] for line in lines])
    return (xy[:,0], xy[:,1])

# Fits the distribution for bond lengths or angles.
def fit_distribution(x, y, term):
    kBT = 300.0 * 8.6173324e-5 # in eV
    # Fit bond distribution.
    if term=='b':        
        xlab = 'Bond length'
        p    = [0.6, 5.0, 1]
        g_x  = lambda p: p[2]*x*x*exp(-p[0]*(x-p[1])**2 / (2.0*kBT)) 
    # Fit angle distribution.
    elif term=='a':
        xlab = 'Bond angle'
        p    = [0.05, 120.0, 1]
        q    = x*pi/180.0
        g_x  = lambda p: p[2]*sin(q)*exp(-p[0]*(q-p[1]*pi/180.0)**2 / (2.0*kBT))

    # Compute parameters, p, that minimize residual, r.
    r = lambda p: y-g_x(p)
    p = leastsq(r, p, maxfev=2000)[0]

    # Output parameters, and plot fit.
    print 'g(x) = f(x) * exp(-a(x-x0)^2)'
    print 'a = %f, x0 = %f' % (p[0],p[1])
    py.plot(x, y, '.', x, g_x(p),'-')
    py.xlabel(xlab)
    py.ylabel('Count')
    py.show()

if __name__ == '__main__': main()

