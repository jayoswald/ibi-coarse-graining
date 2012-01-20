#!/usr/bin/env python
from numpy import *
from scipy.optimize import leastsq
from matplotlib.pyplot import *
import sys

"""
    fitpair.py:
     This script reads a LAMMPS pair table and attempts to fit the energy
     as a combination of a LJ96 potential and Gaussians: a*e^-(x-x0)^2/k.

    To run:
      ./fitpair.py [pair_table_file]

"""

# Number of Gaussian functions to use.
num_gaussian = 9


# Computes the energy of a Lennard Jones 96 potential.
def lj96(v,x):
    sig = v[0]
    A,B = v[1], v[2]
    xs = x-sig
    return A*xs**-9 + B*xs**-6

# Computes the force from a Lennard Jones 96 potential.
def lj96force(v,x):
    sig = v[0]
    A,B = v[1], v[2]
    xs = x-sig
    return 9.0*A*xs**-10 + 6.0*B*xs**-7  

# Computes the energy of a Gaussian potential.
def gauss(v, x):
    e = 0.0
    for i in range(len(v)/3):
        x0 = v[3*i]
        a  = v[3*i+1]
        k  = v[3*i+2]
        xs = (x-abs(x0))**2
        e += a*exp(-xs/(k**2))
    return e

# Computes the force from a Gaussian potential.
def gaussforce(v, x):
    f = 0.0
    for i in range(len(v)/3):
        x0 = v[3*i]
        a  = v[3*i+1]
        k  = v[3*i+2]
        xs = (x-abs(x0))**2
        f += a*exp(-xs/(k**2))*(k**-2)*2.0*(x-abs(x0))
    return f

# Read the data file.
try:    
    data = genfromtxt(sys.argv[1], skip_header=2)
except:
    print 'Reading pair.table.3'
    data = genfromtxt('pair.table.3', skip_header=2)
pr   = data[:,3] < 5.0

r = data[pr,1]
e = data[pr,2]
f = data[pr,3]

error1 = lambda v: lj96(v,r) - e

vlj = [1.0, 1.0, 1.0]  
vlj,status = leastsq(error1, vlj)

error2 = lambda v: lj96(vlj, r) + gauss(v, r) - e

v = []
for x in linspace(2,10,num_gaussian):
    v += [x, 1.0, 1.0]

v,status = leastsq(error2, v, maxfev=8000)
efit = gauss(v, r) + lj96(vlj, r)
print 'LJ parameters'
print reshape(vlj,(1,3))
print 'Gaussian parameters'
print reshape(v,(-1,3))

figure(1)
plot(r, e, 'g-', r, efit, 'b-')


ffit = gaussforce(v,r) + lj96force(vlj,r)
figure(2)
plot(r,f, 'g-', r, ffit, 'b-')
show()
