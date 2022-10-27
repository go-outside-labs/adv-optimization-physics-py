import math
import numpy

# perform 5-point Gauss-Legendre quadrature on a degree 9 (2 N - 1)
# polynomial, and compare to the compound Simpson's rule.

degree = 9


# integrand
def f(x):
    p = 0.0

    n = 0
    while (n <= degree):
        p += x**n
        n += 1

    return p


# analytic (true) integral
def true(a, b):
    t = 0.0

    n = 0
    while (n <= degree):
        t += (b**(n+1) - a**(n+1))/(n+1)
        n += 1

    return t


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


def x(z, a, b):
    """ convert from [-1, 1] (the integration range of Gauss-Legendre)
        to [a, b] (our general range) through a change of variables z
        -> x """

    return 0.5*(b + a) + 0.5*(b - a)*z


# integration limits
a = 0.0
b = 1.0

# we are doing 5-point quadrature for all methods.  delta is the width
# of the slab (5 points = 4 slabs)
delta = (b - a)/4


# Simpson's
I_S = simp(a, b, f, 4)


# Gauss-Legendre

# we need to convert from [-1, 1] (the range in which the roots are
# found) to [a, b] (the range in which our integrand is defined), so
# convert the roots z1, z2, ...

z1 = -math.sqrt(5.0 + 2.0*math.sqrt(10.0/7.0))/3.0
x1 = x(z1, a, b)
w1 = (322.0-13.0*math.sqrt(70.0))/900.0

z2 = -math.sqrt(5.0 - 2.0*math.sqrt(10.0/7.0))/3.0
x2 = x(z2, a, b)
w2 = (322.0+13.0*math.sqrt(70.0))/900.0

z3 = 0.0
x3 = x(z3, a, b)
w3 = 128.0/225.0

z4 = math.sqrt(5.0 - 2.0*math.sqrt(10.0/7.0))/3.0
x4 = x(z4, a, b)
w4 = (322.0+13.0*math.sqrt(70.0))/900.0

z5 = math.sqrt(5.0 + 2.0*math.sqrt(10.0/7.0))/3.0
x5 = x(z5, a, b)
w5 = (322.0-13.0*math.sqrt(70.0))/900.0


# 5-point Gauss-Legendre quadrature -- note the factor in the front
# is a result of the change of variables from x -> z
integral = 0.5*(b-a)*( w1*f(x1) + w2*f(x2) + w3*f(x3) + w4*f(x4) + w5*f(x5) )

print "exact:                  ", true(a,b)
print "5-point Simpson's:      ", I_S, I_S-true(a,b)
print "5-point Gauss-Legendre: ", integral, integral-true(a,b)

