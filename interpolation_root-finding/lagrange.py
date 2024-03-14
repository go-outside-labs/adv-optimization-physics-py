# lagrange interpolation example
# author: MIA STEIN

import math
import numpy
import pylab

# globals to control some behavior
func_type = "tanh"   # can be sine or tanh
points = "fixed"  # can be variable or fixed

npts = 15

def fun_exact(x):
    """ the exact function that we sample to get the points to
    interpolate through """

    if func_type == "sine":
        return numpy.sin(x)
    elif func_type == "tanh":
        return 0.5*(1.0+numpy.tanh((x-1.0)/0.1))


def get_interp_points(N, xmin, xmax):
    """ get the x points that we interpolate at """
    if points == "fixed":
        x = numpy.linspace(xmin, xmax, N)

    elif points == "variable":
        # the Chebyshev nodes
        x = 0.5*(xmin + xmax) + \
            0.5*(xmax - xmin)*numpy.cos(2.0*numpy.arange(N)*math.pi/(2*N))

    return x
    
    
def lagrange_poly(x, xp, fp):
    """ given points (xp, fp), fit a lagrange polynomial and return
        the value at point x """

    f = 0.0
    
    # sum over points
    m = 0
    while (m < len(xp)):

        # create the Lagrange basis polynomial for point m        
        l = None

        n = 0
        while (n < len(xp)):
            if n == m:
                n += 1
                continue

            if l == None:
                l = (x - xp[n])/(xp[m] - xp[n])
            else:
                l *= (x - xp[n])/(xp[m] - xp[n])

            n += 1

        
        f += fp[m]*l

        m += 1

    return f


if func_type == "sine":
    xmin = 0.0
    xmax = 2.0*math.pi
elif func_type == "tanh":
    xmin = 0.0
    xmax = 2.0




# xp, fp are the points that we build the interpolant from
xp = get_interp_points(npts, xmin, xmax)
fp = fun_exact(xp)


# xx are the finely grided data that we will interpolate at to get
# the interpolated function values ff
xx = numpy.linspace(xmin, xmax, 200)
ff = numpy.zeros(len(xx))

n = 0
while (n < len(xx)):
    ff[n] = lagrange_poly(xx[n], xp, fp)
    n += 1


# exact function values at the interpolated points
fexact = fun_exact(xx)


# error
e = fexact-ff


pylab.subplot(211)

pylab.scatter(xp, fp, marker="x", color="r", s=30)
pylab.plot(xx, ff, color="k")

pylab.plot(xx, fexact, color="0.5")

pylab.xlim(xmin, xmax)


pylab.subplot(212)

pylab.plot(xx, e)

pylab.xlim(xmin, xmax)


pylab.savefig("lagrange.png")

