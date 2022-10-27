#!/usr/bin/env python

"""

an example of using the multigrid class to solve Laplace's equation.  Here, we
solve

u_xx = sin(x)
u = 0 on the boundary [0,1]

The analytic solution is u(x) = -sin(x) + x sin(1)

"""
#from io import *
import numpy
import multigrid
import pylab

# the analytic solution
def true(x):
    return -numpy.sin(x) + x*numpy.sin(1.0)


# the L2 error norm
def error(myg, r):

    # L2 norm of elements in r, multiplied by dx to
    # normalize
    return numpy.sqrt(myg.dx*numpy.sum((r[myg.ilo:myg.ihi+1]**2)))


# the righthand side
def f(x):
    return numpy.sin(x)

                
# test the multigrid solver
nx = 64


# create the multigrid object
a = multigrid.ccMG1d(nx,
                     xlBCtype="dirichlet", xrBCtype="dirichlet",
                     verbose=1)

# initialize the solution to 0
init = a.solnGrid.scratchArray()

a.initSolution(init)

# initialize the RHS using the function f
rhs = f(a.x)
a.initRHS(rhs)

# solve to a relative tolerance of 1.e-11
a.solve(rtol=1.e-11)

# alternately, we can just use smoothing by uncommenting the following
#a.smooth(a.nlevels-1,50000)

# get the solution 
v = a.getSolution()

# compute the error from the analytic solution
b = true(a.x)
e = v - b

print " L2 error from true solution = %g\n rel. err from previous cycle = %g\n num. cycles = %d" % \
      (error(a.solnGrid, e), a.relativeError, a.numCycles)

