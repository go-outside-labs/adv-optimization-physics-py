# example from Pang, Ch. 4
#
# Solve a Poisson equation via shooting
#
# u'' = -0.25*pi**2 (u + 1)
#
# with u(0) = 0, u(1) = 1
#
# this has the analytic solution: u(x) = cos(pi x/2) + 2 sin(pi x/2) - 1
#
# M. Zingale (2013-02-18)

import numpy
import math
import pylab


# for plotting different colors
colors = ["k", "r", "g", "b", "c"]

def rk4(y1_0, y2_0, rhs, xl = 0.0, xr = 1.0, n=100):
    """ R-K 4 integration:
           y1_0 and y2_0 are y1(0) and y2(0)
           rhs is the righthand side function
           xl and xr are the domain limits
           n is the number of integration points """

    x = numpy.linspace(xl, xr, n)
    h = x[1] - x[0]  # stepsize

    y1 = numpy.zeros(n)
    y2 = numpy.zeros(n)

    # left boundary initialization
    y1[0] = y1_0
    y2[0] = y2_0

    m = 0
    while (m < n-1):

        dy1dx_1, dy2dx_1 = rhs(y1[m], y2[m])
        dy1dx_2, dy2dx_2 = rhs(y1[m] + 0.5*h*dy1dx_1, y2[m] + 0.5*h*dy2dx_1)
        dy1dx_3, dy2dx_3 = rhs(y1[m] + 0.5*h*dy1dx_2, y2[m] + 0.5*h*dy2dx_2)
        dy1dx_4, dy2dx_4 = rhs(y1[m] + h*dy1dx_3, y2[m] + h*dy2dx_3)

        y1[m+1] = y1[m] + (h/6.0)*(dy1dx_1 + 2.0*dy1dx_2 + 2.0*dy1dx_3 + dy1dx_4)
        y2[m+1] = y2[m] + (h/6.0)*(dy2dx_1 + 2.0*dy2dx_2 + 2.0*dy2dx_3 + dy2dx_4)
    
        m += 1

    return y1, y2


def rhs(y1, y2):
    """ RHS function.  Here y1 = u, y2 = u' 
        This means that our original system is:
           y2' = u'' = -0.25*pi**2 (u+1) """

    dy1dx = y2
    dy2dx = -0.25*math.pi**2 * (y1 + 1.0)

    return dy1dx, dy2dx



def analytic(x):
    """ analytic solution """
    return numpy.cos(math.pi*x/2) + 2.0*numpy.sin(math.pi*x/2) - 1.0

    
# shoot from x = 0 to x = 1.  We will do this by selecting a boundary
# value for y2 and use a secant method to adjust it until we reach the
# desired boundary condition at y1(1)

# number of integration points
npts = 32

# desired tolerance
eps = 1.e-8

# initial guess
y1_0 = 0.0   # this is the correct boundary condition a x = 0
y2_0 = 0.0   # this is what we will adjust to get the desired y1(1)

# desired right BC, y1(1)
y1_1_true = 1.0

# integrate
y1_old, y2_old = rk4(y1_0, y2_0, rhs, xl = 0.0, xr = 1.0, n=npts)

x = numpy.linspace(0.0, 1.0, npts)
pylab.scatter(x, y1_old, label="initial guess", marker="x",c=colors[0])

# new guess -- we don't have any info on how to compute this yet, so
# just choose something
y2_m1 = y2_0   # store the old guess
y2_0 = -1.0

# Secant loop
dy = 1000.0    # fail first time through

# keep track of iteration for plotting
iter = 1

while (dy > eps):
    
    # integrate
    y1, y2 = rk4(y1_0, y2_0, rhs, xl = 0.0, xr = 1.0, n=npts)

    pylab.scatter(x, y1, label="iteration %d" % (iter), marker="x", 
                  c=colors[iter%len(colors)])

    # do a Secant method to get 
    # let eta = our current y2(0) -- this is what we control
    # we want to zero f(eta) = y1_1_true(1) - y1_1(eta)
    
    # derivative (for Secant)
    dfdeta = ( (1.0 - y1_old[npts-1]) - (1.0 - y1[npts-1]) ) / (y2_m1 - y2_0)

    # correction by f(eta) = 0 = f(eta_0) + dfdeta deta 
    deta = -(1.0 - y1[npts-1])/dfdeta

    y2_m1 = y2_0
    y2_0 += deta

    dy = abs(deta)

    y1_old = y1
    y2_old = y2

    iter += 1



pylab.plot(x, analytic(x), color="0.5", label="analytic")

leg = pylab.legend(loc=2)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.xlim(0.0, 1.0)

pylab.savefig("shoot.png")

