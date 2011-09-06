#!/usr/bin/env python
from numpy import * 
from scipy import optimize
import interpolator
import distribution 
from matplotlib import pyplot as py

# Boltzmann constant in kcal/mol/K
kB = 0.0019872041 
lj96_force = lambda v,x: v[1]*(18.0*v[0]**9/x**10 - 18.0*v[0]**6/x**7)

"""
"""
class PairTable:
    # 
    def __init__(self, md_temp, md_rdf, plot=False, npts=1000):
        self.temperature     = md_temp 
        self.min_distance    = 2.0
        self.all_atom_rdf = distribution.average_rdf('r', md_rdf)
        
        self.distance  = []
        self.force     = []
        self.rawforce  = []
        self.energy    = []

        print 'Computing initial pair table'
        self.compute(self.all_atom_rdf)
        
        if plot: 
            self.plot_force()
            self.plot_energy()
            py.show()

    # Computes the pair table, returning [distance, energy, force]
    def compute(self, rdf, npts=1000):
        # Removes any values with a density of zero.
        # This gets added back later.
        rdf = rdf[nonzero(rdf[:,1])[0], :]
        # Get distance (r) and energy (e).
        r =  rdf[:,0]
        e = -kB*self.temperature*log(rdf[:,1])

        # Compute derivative from splines.
        rr = linspace(r[0], r[-1], npts)
        interp = interpolator.Interpolator(r,e)
        ff = [-interp.derivative(ri) for ri in rr]

        # Set minimum radial value where force is most repulsive.
        i0 = ff.index(max(ff))
        ff,rr = ff[i0::],rr[i0::]

        lj96_err = lambda v,x,y: y - lj96_force(v,x)

        # Subtract 96 part away and smooth.
        v0 = [5.0, 0.01] 
        v = optimize.leastsq(lj96_err, v0, args=(rr, ff))[0]
        self.rawforce.append(ff[:])
        ff -= lj96_force(v, rr)
        w,K = 2,4
        for k in range(K):
            for i in range(w,len(ff)-w):
                ff[i] = mean(ff[i-w:i+w])

        ff += lj96_force(v, rr)

        # Add more values to short distance to make sure that 
        # LAMMPS run won't fail when pair distance < table min.
        dr   = rr[1]-rr[0]
        rpad = arange(rr[0]-dr, self.min_distance, -dr)[::-1]
        fpad = lj96_force(v, rpad) + ff[0] - lj96_force(v, rr[0])
        rr   = concatenate((rpad, rr))
        ff   = concatenate((fpad, ff))
        self.rawforce[-1] = concatenate((fpad, self.rawforce[-1]))

        # Make ff die off smoothly at rcut.
        ff -= (ff[-1]/rr[-1]) / rr
        ff *= exp(-1.0/(rr[-1] - rr  + 1e-20))

        # Compute energy by integrating forces.
        # Integrating backwards reduces noise.
        ee = -simpson_integrate(rr[::-1], ff[::-1])[::-1]
        ee -= ee[-1]

        self.distance.append(rr)
        self.force.append(ff)
        self.energy.append(ee)

    # Plots the forces at an iteration.
    def plot_force(self, it=-1):
        r = self.distance[it]
        f = self.force[it]

        py.figure()
        py.hold(1)
        py.plot(r, self.rawforce[it], linewidth=2, color='#CCCCEE')
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


# Computes the corrected pair table.
def corrected_pair_table(T, ff0, npts=1000):
    rri,eei,ffi = pair_table(T, '', npts)
    rr, ee, ff  = pair_table(T, 'rdf', npts)

    ffi_int = interpolator.Interpolator(rri, ffi)
    ff_int  = interpolator.Interpolator(rr, ff)

    df = array([ffi_int(r) - f for r,f in zip(rr,ff)])
    # Remove any sharp kinks.
    for i in range(1,len(df)-1):
        if df[i] == max(df[i-1:i+1]) or df[i] == min(df[i-1:i+1]):
            df[i] = mean(df[i-1:i+1])

    force = ff0 + df

    energy = -simpson_integrate(rr[::-1], force[::-1])[::-1]
    energy -= energy[-1]

    return rr, force, energy

# Cumulative integration of f using Simpson's rule.
def simpson_integrate(x,f):
    F = zeros((len(f)))
    F[0] = 0.0
    F[1] = 0.5*(f[0]+f[1]) * (x[1]-x[0])
    for i in range(2,len(f)):
        # Integral is from i-2 to here + F[i-2]
        F[i] = F[i-2] + (f[i-2]+4.0*f[i-1]+f[i])*(x[i]-x[i-2])/6.0
    return F

# If this script was called as top level, run main.
if __name__=='__main__': 
    pair_table()


