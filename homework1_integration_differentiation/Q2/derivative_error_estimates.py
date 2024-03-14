#################################################################################
#        DERIVATIVE ERROR ESTIMATION AND THE RICHARDSON EXTRAPOLATION METHOD
#    This program computes:
#    1) the Richardson derivative of f(x) = sin(x) at x = 1 using an adaptative scheme,
#    2) for a relative error < 1e-7,
#                                                                          
#    (MIA STEIN, spring 2013)
#
#################################################################################



import math
import numpy
import pylab



""" Functions """

def func(x):
    """ the function to be plotted """
    fx = numpy.sin(x)
    return fx


def fprime(x):
    """ the analytic derivative of func(x) """
    fp = numpy.cos(x)
    return fp


def machine_precision():
    """ calculates the machine precision, 2*eps"""
    x = 1.0
    eps = 1.0
    while (not x + eps == x):
        eps = eps/2.0
    return 2.0*eps


def delta1(h, x):            
    """ calculates the second-order centered difference """
    d1 = (func(x+h) - func(x-h))/(2.0*h)
    return d1


def relative_error(fa, fn):            
    """ calculates the relative error """
    e = abs( ( fa - fn ) / fa )
    return e


def absolute_error(fa, fn):            
    """ calculates the relative error """
    e = abs(  fa - fn  )
    return e


def richardson_derivative(h, x):
    """ calculate the derivative using delta1 for h and h/2 """
    d = delta1(h, x)
    h = h/2
    d1 = delta1(h, x)
    f1 =  (4*d1-d)/3

    # MZ -- why is there an h**2 here?  
    epsilon = h*h*abs(d-d1) 
    return f1, epsilon


def printing(n, h, f1, fp, e, er, ea):
    print "\n________________________________________________________________________________"
    print "Iteration n = ", int(n)
    print "The distance between two points in this grid (h) is:    ", h            
    print "Calculated (Richardsob) Derivative of f(x)=sin(x) at x = 1:    ", f1
    print "Analytic First Derivative of f(x)=sin(x) at x = 1:    ", fp
    print "Truncated Error (from the Taylor approx):    ", e
    print "Relative Error:    ", er
    print "Absolute Error:    ", ea
    return 0
    
    
""" Variables """

CONST_EPS = 10**(-6)            # relative error
CONST_XL = 0.0              # most left point
CONST_XR = math.pi          # most right point: 3.14 > x=1,
CONST_X = 1.0               # point where we are calculating the derivative

h_array = []
er_array = []
ea_array = []
epsilon_array = []


epsilon = 1.1*10**(-6) 
n = 0.0



""" Main Program """
while (epsilon >= CONST_EPS):
    n += 1
    h = abs((CONST_XR - CONST_XL)/(2*float(n)))        # distance in the grid

    # MZ -- the error you are asked to monitor is the error between 
    # the difference with h and the one with h/2
    f1, epsilon = richardson_derivative(h, CONST_X)
    error_r = relative_error(f1, fprime(CONST_X))
    error_a = absolute_error(f1, fprime(CONST_X))

    printing(int(n), h, f1, fprime(CONST_X), epsilon, error_r, error_a)
    
    h_array.append(h)
    er_array.append(error_r)
    ea_array.append(error_a)
    epsilon_array.append(epsilon)




""" Plotting the Errors vs h """
pylab.loglog(h_array, er_array, 'go',  label="Relative Error")
pylab.loglog(h_array, ea_array, 'r*',  label="Absolute Error")
pylab.loglog(h_array, epsilon_array,'b--',  label="Truncated Error $O(h^4)$")

leg = pylab.legend(loc=4,labelspacing=0.1)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(1)

pylab.xlabel('Log( Griding Distance h)')
pylab.ylabel('Log (Error)')
pylab.grid(True)

pylab.savefig("der_error.png")



""" Plotting the Errors vs n """
pylab.clf()
pylab.cla()
nlist = list(range(1, int(n)+1, 1))
pylab.loglog(nlist, er_array, 'go',  label="Relative Error")
pylab.loglog(nlist, ea_array, 'r*',  label="Absolute Error")
pylab.loglog(nlist, epsilon_array,'b--',  label="Truncated Error $O(h^4)$")

leg = pylab.legend(loc=1,labelspacing=0.1)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(1)

pylab.xlabel('Log (Number of Steps n)')
pylab.ylabel('Log (Error)')
pylab.grid(True)

pylab.savefig("der_error2.png")





print "\nMachine precision: ", machine_precision()
print "Done!"
