# Simpson's rule -- odd intervals

import math
import numpy
import pylab

# function we wish to integrate
def fun(x):
    return numpy.sin(math.pi*x)


# analytic value of the integral of fun() above
def I_exact(a,b):
    return (1.0/math.pi)*(-math.cos(math.pi*b) + math.cos(math.pi*a))


# do a Simpson's integration by breaking up the domain [a,b] into N
# slabs.  Note: we consider the N odd case here:
def simp(a,b,f,N):

    xedge = numpy.linspace(a,b,N+1)

    integral = 0.0

    if N%2 == 0:
        M = N
        odd = 0
    else:
        M = N - 1
        odd = 1

    delta = (xedge[1] - xedge[0])

    n = 0
    while n < M:
        integral += (1.0/3.0)*delta*(f(xedge[n]) + 
                                     4.0*f(xedge[n+1]) + 
                                     f(xedge[n+2]))
        n += 2

    # if we had an odd # of bins, do the last one:
    if (odd):
        integral += (delta/12.0)*(-f(xedge[N-2]) + 
                                   8.0*f(xedge[N-1]) + 
                                   5.0*f(xedge[N]))

    return integral


a = 0.0
b = 1.0

N = [3, 7, 15, 31]

delta = []
error = []

for n in N:
    t = simp(a,b,fun,n)
    e = t - I_exact(a,b)
    dx = (b - a)/n

    print dx, t, e
    
    # for plotting
    delta.append(dx)
    error.append(e)


d = numpy.array(delta)
e = numpy.array(error)
pylab.scatter(d, e)

# plot a line representing delta**4 scaling
pylab.plot(d, e[0]*(d/d[0])**4)

ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.xlabel(r"$\delta$")
pylab.ylabel("absolute error")

pylab.savefig("simp-odd.eps")



 
