#!/usr/bin/env python
from numpy import * 
from scipy import optimize
import distribution 
from interpolator import Interpolator
from matplotlib import pyplot as py

# Boltzmann constant in kcal/mol/K
kB = 0.0019872041 
# Derivative of a 9/6 Lennard Jones potential.
# as in http://lammps.sandia.gov/doc/pair_lj96.html
# v[0]: Sigma, v[1]: Epsilon
lj96_force = lambda v,r: v[1]*(36.0*v[0]**9/r**10 - 24.0*v[0]**6/r**7)

"""
  Computes potential energy and force vs. distance from the radial density function.
"""
class PairTable:
    # Initializes the pair table from all-atom data.
    def __init__(self, md_temp, md_rdf, plot=False):

        # Sets the temperature of the all-atom simulation.
        self.temperature  = md_temp 

        # Sets the smallest distance in the pair table.
        # If LAMMPS simulations crash with pair cutoff error, this needs to be smaller.
        self.min_distance = 0.00001
        self.npts = 1000

        # Computes the average distribution functions for the all-atom case.
        self.all_atom_rdf = distribution.average_rdf(distribution.md_rdf_files())

        # Initializes the pair tables as empty lists.        
        self.distance  = []
        self.force     = []
        self.energy    = []

        # Computes the initial potential with the all-atom data.
        self.compute(self.all_atom_rdf)
        
        # Plot the potential.
        if plot: 
            self.plot_force()
            self.plot_energy()
            py.show()

    # Computes the pair table and appends it to the current step.
    def compute(self, rdf):
        cut_beg    = 13.0
        cut_end    = 15.0

        # Removes values with a density of zero, they will get added back later.
        rdf = rdf[nonzero(rdf[:,1])[0], :]
        # Finds the first index where the density is greater than 0.25*rho.
        i0 = nonzero(rdf[:,1] > 0.25)[0][0]

        # Get distance (r) and energy (e).
        r =  rdf[:,0]
        e = -kB*self.temperature*log(rdf[:,1])


        dr_final = (cut_end-self.min_distance)/self.npts
        # Figure out what range we need

        # Compute derivative from splines.
        rr     = arange(cut_end, r[i0], -dr_final)[::-1]

        interp = Interpolator(r,e)
        ff = [-interp.derivative(ri) for ri in rr]

        # Subtract 96 part away and smooth.
        v0 = [5.0, 0.01] 
        lj96_err = lambda v: ff - lj96_force(v,rr)
        v = optimize.leastsq(lj96_err, v0)[0]

        ff -= lj96_force(v, rr)
        w,K = 1,2
        for k in range(K):
            for i in range(w,len(ff)-w):
                ff[i] = mean(ff[i-w:i+w])

        ff += lj96_force(v, rr)


        # Add more values to short distance to make sure that 
        # LAMMPS run won't fail when pair distance < table min.
        rpad = arange(self.min_distance, rr[0]-0.5*dr_final, dr_final)
        fpad = lj96_force(v, rpad) + ff[0] - lj96_force(v, rr[0])
        rr   = concatenate((rpad, rr))
        ff   = concatenate((fpad, ff))
    
        # Now we cut off forces smoothly at any point past max_attraction.
        ff[rr>cut_beg] *= 0.5*(1.0+cos(pi*(rr[rr>cut_beg]-cut_beg)/(cut_end-cut_beg)))
        ff[rr>cut_end] = 0.0

        # Compute energy by integrating forces.
        # Integrating backwards reduces noise.
        ee = -simpson_integrate(rr[::-1], ff[::-1])[::-1]
        ee -= ee[-1]

        self.distance.append(rr)
        self.force.append(ff)
        self.energy.append(ee)

    # Writes the pair table data for iteration, it.
    def write_lammps(self, path, key, it):
        r = self.distance[it]
        f = self.force[it]
        e = self.energy[it]
        fid = open(path, 'w')
        fid.write(key+'\n')
        fid.write('N %d R %f %f\n\n' %(len(r), min(r), max(r)))
        for i in range(len(r)):
            fid.write('%d %f %f %f\n' %(i, r[i], e[i], f[i]))
        
    # Plots the forces at an iteration.
    def plot_force(self, it=-1):
        r = self.distance[it]
        f = self.force[it]

        py.figure()
        py.hold(1)
        py.plot(r, f, 'b', linewidth=2)
        py.axis((min(r), max(r), min(f)-0.2, min(f) + 1.0))
        py.hold(0)
        py.xlabel('Pair distance (A)')
        py.ylabel('Force (kcal/mol/Angstrom)')

    # Plots the forces at an iteration.
    def plot_energy(self, it=-1):
        r = self.distance[it]
        e = self.energy[it]
    
        py.figure()
        py.plot(r, e, linewidth=2, color='b')
        py.axis((min(r), max(r), min(e)-0.2, min(e) + 1.0))
        py.xlabel('Pair distance (A)')
        py.ylabel('Energy (kcal/mol)')

    # Computes the corrections to the pair table.
    def correction(self, it, pair):
        # Compute force table based on current iteration.
        rdf = distribution.iteration_rdf_files(it, pair)
        # Appends new force, energy, distance, table.
        self.compute(distribution.average_rdf(rdf))
        # Computes the correction to the force.
        df = self.force[-1] - self.force[0]
        self.force[-1] = self.force[-2] - df

        # Need to reintegrate energy
        rr = self.distance[-1]
        ff = self.force[-1]
        self.energy[-1] = -simpson_integrate(rr[::-1], ff[::-1])[::-1]
        self.energy[-1] -= self.energy[-1][-1]
    
# Cumulative integration of f using Simpson's rule.
def simpson_integrate(x,f):
    F = zeros((len(f)))
    F[0] = 0.0
    F[1] = 0.5*(f[0]+f[1]) * (x[1]-x[0])
    for i in range(2,len(f)):
        # Integral is from i-2 to here + F[i-2]
        F[i] = F[i-2] + (f[i-2]+4.0*f[i-1]+f[i])*(x[i]-x[i-2])/6.0
    return F

