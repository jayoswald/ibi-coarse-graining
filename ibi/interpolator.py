#!/usr/bin/env python
"""
    This module takes a pair of x and y vector data and creates an interpolation
    object that can return the interpolated value or its derivative at any point
    within [min(x), max(x)].
"""
from numpy import linspace, sin

class Interpolator:
    def __init__(self, x, y):
        self.intervals = intervals = []
        # Generate a cubic spline for each interpolation interval.
        for i in range(len(x)-1): 
            u,v   = x[i],x[i+1]
            FU,FV = y[i],y[i+1] 
            # adjust h as needed
            h = 0.001
            if i>0:
                DU = (y[i+1]-y[i-1]) / (x[i+1]-x[i-1])
            else:
                DU = (y[i+1]-y[i]) / (x[i+1]-x[i])

            if i<len(x)-2:
                DV = (y[i+2]-y[i]) / (x[i+2]-x[i])
            else:
                DV = (y[i+1]-y[i]) / (x[i+1]-x[i])

            denom = (u - v)**3

            A = (u-v)*(DV+DU) + 2*(FV-FU)
            B = -(v*v*(-DV - 2*DU) +
                  u*(v*(DU-DV) + 3*(FV-FU)) + 3*v*(FV-FU) + u*u*(2*DV + DU)) 
            C = (- DU*v**3 + u*(v*v*(-2*DV-DU)  + 6*v*(FV-FU)) +
                 v*u*u*(DV + 2*DU) + DV*u**3) 
            D = -(u*(-DU*v**3  - 3*FU*v*v) +
                  FU*v**3 + u*u*(v*v*(DU-DV) + 3*FV*v) + u**3 * (DV*v-FV)) 
            intervals.append((u, A/denom, B/denom, C/denom, D/denom))

    def __call__(self, x):
        u, A, B, C, D = self.get_interval(x)
        # Plug coefficients into polynomial.
        return ((A*x + B)*x + C)*x + D

    def derivative(self,x):        
        u, A, B, C, D = self.get_interval(x, 0, len(self.intervals))
        return 3*A*x*x + 2*B*x + C

    # Tree-search the intervals to get coefficients.
    def get_interval(self, x, beg, end):
        n = end-beg
        if n < 2: 
            return self.intervals[beg]
        mid = beg + n/2
        if x < self.intervals[mid][0]: 
            return self.get_interval(x, beg, mid)
        else: 
            return self.get_interval(x, mid, end)

