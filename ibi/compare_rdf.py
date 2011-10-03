#!/usr/bin/env python

from matplotlib  import pyplot as py
from distribution import average_rdf


md_rdf = average_rdf('r')
cg_rdf = average_rdf('r', '')

md_rdf[:,1] /= md_rdf[-1,1]

py.plot(md_rdf[:,0], md_rdf[:,1])
py.plot(cg_rdf[:,0], cg_rdf[:,1])
py.xlabel('Radial distance ($\AA$)')
py.ylabel('Radial density function g(r)')
py.legend(['All atom model', 'coarse grain model'], loc='lower right')
py.show()


