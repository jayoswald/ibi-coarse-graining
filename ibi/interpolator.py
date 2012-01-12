#!/usr/bin/env python
"""
    This module takes a pair of x and y vector data and creates an interpolation
    object that can return the interpolated value or its derivative at any point
    within [min(x), max(x)].
"""
from numpy import *

class Interpolator:
    # Initializes the spline coefficients.
    # see http://mathworld.wolfram.com/CubicSpline.html for details.
    def __init__(self, x, y):
        rhs       = zeros((len(x)))
        rhs[0]    = y[1]-y[0]
        rhs[1]    = y[2]-y[0]
        rhs[2:-2] = y[3:-1] - y[1:-3]
        rhs[-2]   = y[-1] - y[-3]
        rhs[-1]   = y[-1] - y[-2]

        n  = len(rhs)
        A  = diag([4.0]*n, 0)
        A += diag([1.0]*(n-1), -1) + diag([1.0]*(n-1),  1)
        A[0,0] = A[-1,-1] = 2.0
        self.d = linalg.solve(A, 3.0*rhs)
        self.x = copy(x)
        self.y = copy(y)
        self.npts = len(self.d)-1

    # Computes the interpolated function.
    def __call__(self, x):
        xi,a,b,c,d = self.get_interval(x, 0, self.npts)
        return a + b*xi + c*xi*xi + d*xi*xi*xi

    # Computes the derivative of the interpolated function.
    def derivative(self,x):        
        xi,a,b,c,d = self.get_interval(x, 0, self.npts)
        return (3.0*d*xi*xi + 2.0*c*xi + b)/(self.x[1]-self.x[0])

    # Binary search the intervals to get coefficients.
    def get_interval(self, x, beg, end):
        n = end-beg
        if n < 2: 
            i = beg
            xi = (x - self.x[i]) / (self.x[i+1] - self.x[i])
            a = self.y[i]
            b = self.d[i]
            c = 3.0*(self.y[i+1] - self.y[i])   - 2.0*b - self.d[i+1]
            d = 2.0*(self.y[i]   - self.y[i+1]) +     b + self.d[i+1]
            return (xi,a,b,c,d)

        mid = beg + n/2
        if x < self.x[mid]: 
            return self.get_interval(x, beg, mid)
        else: 
            return self.get_interval(x, mid, end)

