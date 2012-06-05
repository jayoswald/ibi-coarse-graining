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
def lj96_force(v,r):
    f = v[1]*(36.0*v[0]**9/r**10 - 24.0*v[0]**6/r**7)
    for i in range(2,len(v),3):
        f += v[i+2]*exp(-v[i+1]*v[i+1]*(r-v[i])**2)
    return f

"""
  Computes potential energy and force vs. distance from the radial density function.
"""
class PairTable:
    # Initializes the pair table from all-atom data.
    def __init__(self, md_temp, md_rdf):

        # Sets the temperature of the all-atom simulation.
        self.temperature  = md_temp 

        # Sets the smallest distance in the pair table.
        # If LAMMPS simulations crash with pair cutoff error, this needs to be smaller.
        self.min_distance = 0.00001
        self.npts = 1000

        # Pair force correct per unit pressure.
        self.pfactor = 1.0/5000.0
        self.last_p_error = 0.0

        # Computes the average distribution functions for the all-atom case.
        self.all_atom_rdf = distribution.average_rdf(distribution.md_rdf_files())

        # Initializes the pair tables as empty lists.        
        self.distance  = []
        self.force     = []
        self.energy    = []

        # Computes the initial potential with the all-atom data.
        self.allatomcompute(self.all_atom_rdf)
       
    # Computes the pair table and appends it to the current step.


    def allatomcompute(self,rdf):

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

        ff, v = smooth_force(rr, array(ff), len(self.force))

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
    
        v0 = [5.0, 0.01]
        lj96_err = lambda v: ff - lj96_force(v,rr)
        v = optimize.leastsq(lj96_err, v0)[0]


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
    def plot_force(self, name, it=-1):
        r = self.distance[it]
        f = self.force[it]

        # Only use first iteration to compute force
        # so all plots share the same range.
        fmin = 3.0*min(self.force[0])
        fmax = -4.0*fmin

        py.clf()
        py.hold(1)
        py.plot(r, f, 'b')
        py.axis((min(r), max(r), fmin, fmax))
        py.xlabel('Pair distance (A)')
        py.ylabel('Force (kcal/mol/Angstrom)')
        py.savefig(name)

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
        ff = smooth_force(rr, ff, len(self.force))[0]
        self.energy[-1] = -simpson_integrate(rr[::-1], ff[::-1])[::-1]
        self.energy[-1] -= self.energy[-1][-1]

    # Makes a pressure correction to the pair potential.
    def pressure_correction(self, p_error):
        r  = self.distance[-1]
        cut_end = r[-1]

        incr = 1.0
        # Dynamically adjust ratio so that each iteration reduces error by 95%.
        if self.last_p_error != 0.0:
            ratio = (self.last_p_error - p_error)/self.last_p_error
            if ratio > 0.0:
                incr /= ratio + 0.05
                self.pfactor *= incr

        self.last_p_error = p_error
        A = -p_error * self.pfactor

        print 'Applying force correction of:', A
        print 'Correction factor is        :', self.pfactor, '(%.2fx)'%incr

        dV = A*kB*self.temperature*(1.0-r/cut_end)
        dF = A*kB*self.temperature*(1.0/cut_end)
        self.energy[-1] += dV
        self.force[-1]  += dF 
    
# Cumulative integration of f using Simpson's rule.
def simpson_integrate(x,f):
    F = zeros((len(f)))
    F[0] = 0.0
    F[1] = 0.5*(f[0]+f[1]) * (x[1]-x[0])
    for i in range(2,len(f)):
        # Integral is from i-2 to here + F[i-2]
        F[i] = F[i-2] + (f[i-2]+4.0*f[i-1]+f[i])*(x[i]-x[i-2])/6.0
    return F

# Smooths out the computed force.
def smooth_force(r, f, it):
    mask = f < -min(f)*8.0
    rm,fm = r[mask],f[mask]

    error = lambda p: fm - lj96_force(p, rm)
    p0 = [5.0,0.01]
#    for d in [5.0,6.0,8.0,10]:
#        p0 += [d, 2.0, 0.0]

    fit = optimize.leastsq(error, p0, maxfev=4000, full_output=1)
    while len(p0) > 2 and fit[2]['nfev'] == 4000:
        p0 = p0[0:-3]
        fit = optimize.leastsq(error, p0, maxfev=4000, full_output=1)
    p = fit[0]
    
    print 'Gaussian peaks    at  :', p[2::3]
    print 'Gaussian widths  are:', array(p[3::3])**-2
    print 'Gaussian heights are:', array(p[4::3])

    py.clf()
    py.plot(rm, fm, '.')
    py.plot(rm, lj96_force(p,rm))

    fm -= lj96_force(p,rm)
    w,K = 1,2
    for k in range(K):
        for i in range(w,len(fm)-w):
            fm[i] = mean(fm[i-w:i+w+1])

    f[mask]    = fm + lj96_force(p,rm)
    f[mask==0] = lj96_force(p, r[mask==0]) + fm[0] 

    py.plot(rm,f[mask])
    py.legend(['original', 'fitted', 'smoothed'], loc='upper right')
    py.savefig('smooth-%d.png' %it)

    return f,p


