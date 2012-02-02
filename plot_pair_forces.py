#!/usr/bin/env python
from matplotlib.pyplot import *
from numpy import *
from glob import glob

pair_tables = glob('pair.table.*')



data = [ genfromtxt(f, skip_header=3) for f in pair_tables]
for i,d in enumerate(data):
    prange = d[:,3] < 5.0
    plot(d[prange,1], d[prange,2], label='step %d'%i)

legend()    
show()





