# Simpson's rule


import math
import numpy
import sys

# function we wish to integrate
def fun(x):
    return numpy.exp(-x)


# analytic value of the integral
def I_exact(a,b):
    return -math.exp(-b) + math.exp(-a)


# do a Simpson's integration by breaking up the domain [a,b] into N
# slabs.  Note: N must be even, because we do a pair at a time
def simp(a,b,f,N):

    xedge = numpy.linspace(a,b,N+1)

    integral = 0.0

    if not N%2 == 0:
        sys.exit("ERROR: N must be even")

    delta = (xedge[1] - xedge[0])

    n = 0
    while n < N:
        integral += (1.0/3.0)*delta*(f(xedge[n]) + 
                                     4.0*f(xedge[n+1]) + 
                                     f(xedge[n+2]))
        n += 2

    return integral


a = 0.0
b = 1.0

N = 2
while (N <= 128):
    t = simp(a,b,fun,N)
    e = t - I_exact(a,b)
    print N, t, e

    N *= 2


 
