#########################################################################
#            GAUSSIAN QUADRATURE                   
#                                  
# This module (based on Mike Zingale's code)
#    1) consider a 5-point quadrature and a 9th polynomial,
#    2) compute the integral using the gauss-legendre method and the
# compound version of the simpsons rules (two pairs of interval)
#    3) compute the errors against the exact integral.
#
#    (MIA STEIN, spring 2013)
#########################################################################
import math
import numpy
import pylab



""" Functions """
def func(x):
    """ function to integrate: a 9th degree polynomial"""
    fx = x**9 + x**8 + x**7 + x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
    return fx


def exact_int(a,b):
    """ analytic value of the integral """
    I = (1/10.0)*(b-a)**10 + (1/9.0)*(b-a)**9 + (1/8.0)*(b-a)**8 + (1/7.0)*(b-a)**7 + (1/6.0)*(b-a)**6 + (1/5.0)*(b-a)**5 + (1/4.0)*(b-a)**4 + (1/3.0)*(b-a)**3 + (1/2.0)*(b-a)**2 +  (b-a)
    return I


def interval_converter(z, a, b):
    """ convert from [-1, 1] (the integration range of Gauss-Legendre)
        to [a, b] (our general range) through a change of variables z -> x """
    z1 = 0.5*(b + a) + 0.5*(b - a)*z
    return z1


def gauss_legen_int_5_points(x, w, a, b, f):
    """ calculate the gauss legendre integral, where 0.5 is from the interval transformation"""
    Igs = 0.5*(b-a)*( w[0]*f(x[0]) + w[1]*f(x[1]) + w[2]*f(x[2]) +  w[3]*f(x[3]) +  w[4]*f(x[4]) )
    return  Igs


def simpson_int(a,b,f,N):
    """"do a Simpson's integration by breaking up the domain [a,b] into N, considering N is even."""
    xedge = numpy.linspace(a,b,N)
    delta = xedge[1] - xedge[0]
    Is = 0.0
    n = 0
    while n < N-1:
        Is += (delta/3.0)*(f(xedge[n]) + 4.0*f(xedge[n+1]) + f(xedge[n+2]))
        n += 2  
    
    print "N = ", N
    # you only need to add the odd if you are indeed odd!
    Is_odd = (delta/12.0)*(-f(xedge[b - 2.0*delta]) + 8.0*f(xedge[b - delta]) + 5.0*f(xedge[b]))    # this only helps with larger Ns
    Is_with_odd = Is + Is_odd 
    
    Is_chopped = 0.0
    n = 0
    while n < N-1:
        Is_chopped  += (delta/3.0)*(f(xedge[n]) + 4.0*f(xedge[n+1]) + f(xedge[n+2]))
        n += 2
    
    return Is_chopped, Is_with_odd


def trap_int(a,b,f,N):
    xedge = numpy.linspace(a,b,N)
    delta = xedge[1] - xedge[0]
    It = 0.5*delta*(f(a) + f(0.5*(a+b))) + 0.5*delta*(f(0.5*(a+b)) + f(b))
    return It


def print_results(eags, eas_odd, eas_chopped,Is_chopped,I, Is_with_odd, Igs, It):
    print "\nAnalytical value for the integral: ", I
    print "Numerical Integration by the Simpson's rule chopping the odd part: ", Is_chopped
    print "Numerical Integration by the Simpson's rule with odd approximation: ", Is_with_odd
    print "Numerical Integration by the Gauss-Legendre method: ", Igs
    print "Numerical Integration by Trapezoidal method: ", It
    print "\nAbsolute error to the Simpson's rule integral (without odd): ", eas_chopped
    print "Absolute error to the Simpson's rule integral (with odd): ", eas_odd
    print "Absolute error to the Gauss-Legendre integral: ", eags
    return 0





""" Variables """
CONST_A = 0.0
CONST_B = 1.0

CONST_NUMBER_POINTS_QUAD = 5.0            # for 5-point quadrature
CONST_DELTA = (CONST_B - CONST_A)/(CONST_NUMBER_POINTS_QUAD-1.0)      # witdth if slab

ARRAY_ROOTS = [0, (1.0/3.0)*math.sqrt( 5.0 - 2.0*math.sqrt( 10.0/7.0 ) ), - (1.0/3.0)*math.sqrt( 5.0 - 2.0*math.sqrt( 10.0/7.0) ), (1.0/3.0)*math.sqrt( 5.0 + 2.0*math.sqrt( 10.0/7.0 ) ), - (1.0/3.0)*math.sqrt( 5.0 + 2.0*math.sqrt( 10.0/7.0) ) ]
ARRAY_WEIGHTS = [ 128.0/225.0 , (322.0 + 13.0*math.sqrt(70.0))/(900.0) , (322.0 + 13.0*math.sqrt(70.0))/(900.0),  (322.0 - 13.0*math.sqrt(70.0))/(900.0), (322.0 - 13.0*math.sqrt(70.0))/(900.0)]

new_array_roots = []



""" Main Code """
# convert all the roots to the interval [-1,1]
for i in range(len(ARRAY_ROOTS)):
    x = interval_converter(  ARRAY_ROOTS [i], CONST_A , CONST_B )
    new_array_roots.append(x)


# calculate the gauss legendre, the simpson, and the exact integral
Igs = gauss_legen_int_5_points(new_array_roots, ARRAY_WEIGHTS, CONST_A , CONST_B, func)
I = exact_int( CONST_A , CONST_B)
Is_chopped, Is_with_odd = simpson_int(CONST_A , CONST_B, func, CONST_NUMBER_POINTS_QUAD)
It = trap_int(CONST_A , CONST_B, func, CONST_NUMBER_POINTS_QUAD)


# calculates the errors
eags = abs(I - Igs)
eas_chopped = abs(I - Is_chopped)
eas_odd = abs(I - Is_with_odd)

#far from the analytical?
print_results(eags, eas_odd, eas_chopped,  Is_chopped,I, Is_with_odd, Igs, It)



print "\nDone!"
