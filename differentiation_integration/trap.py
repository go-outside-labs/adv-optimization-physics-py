# trapezoidal rule


import math
import numpy

# function we wish to integrate
def fun(x):
    return numpy.exp(-x)


# analytic value of the integral
def I_exact(a,b):
    return -math.exp(-b) + math.exp(-a)


# do a trapezoid integration by breaking up the domain [a,b] into N
# slabs
def trap(a,b,f,N):

    xedge = numpy.linspace(a,b,N+1)

    integral = 0.0

    n = 0
    while n < N:
        integral += 0.5*(xedge[n+1] - xedge[n])*(f(xedge[n]) + f(xedge[n+1]))
        n += 1

    return integral


N = 3
a = 0.0
b = 1.0

N = 2
while (N <= 128):
    t = trap(a,b,fun,N)
    e = t - I_exact(a,b)
    print N, t, e

    N *= 2


 
