#########################################################################
#            SIMPSON'S RULES                    
#                                  
# This module (based on Mike Zingale's code)
#    1) calculates the integral of f(x) = sin (pi x) over the 
#     interval [0,1], using N = 3,7,15 and 31 slabs/intervals, for odd Ns.
#    2) Plots the absolute error vs delta = (b-a)/N.
#
#    (MIA STEIN, spring 2013)
#########################################################################

import math
import numpy
import pylab



""" Functions """

def func(x):
    """ function to integrate"""
    fx = numpy.sin(math.pi*x)
    return fx


def exact_int(a,b):
    """ analytic value of the integral """
    I = (numpy.cos(math.pi*a) - numpy.cos(math.pi*b ))/math.pi
    return I


def simpson_int_odd(a,b,f,N): 
    """ calculates remaining odd slabs """
    xedge = numpy.linspace(a,b,N+1)
    delta = xedge[1] - xedge[0]
    Is_odd = (delta/12.0)*(-f(xedge[N-2]) + 8.0*f(xedge[N-1]) + 5.0*f(xedge[N]))
    return Is_odd


def simpson_int(a,b,f,N):
    """"do a Simpson's integration by breaking up the domain [a,b] into N, considering N is even."""
    # MZ -- you were passing in the wrong N, so the xedge was wrong
    xedge = numpy.linspace(a,b,N+1)
    delta = xedge[1] - xedge[0]
    Is = 0.0
    n = 0
    # MZ: with the proper N, this loop executed too many times
    while n < N-1:
        Is += (delta/3.0)*(f(xedge[n]) + 4.0*f(xedge[n+1]) + f(xedge[n+2]))
        n += 2

    return Is


def printing(slab, I, Is, ea):
    """ print output """
    print "\nNumber of slabs to be Integrated: ", slab
    print "Simpson's Integral value: ", Is
    print "Analytic Integral value: ", I
    print "Absolute error: ", ea
    return 0




""" Variables """
CONST_A = 0.0
CONST_B = 1.0
CONST_N_SLABS = [3, 7, 15, 31]

ea_array = []
delta_array = []



""" Main Function""" 
""" Since the intervals are odd, we need to do extra calculations in the edges """
for i in range (len(CONST_N_SLABS)):
    slab = CONST_N_SLABS[i]
    
    if not slab%2 == 0:
        Is_odd = simpson_int_odd(CONST_A, CONST_B, func, slab)
        edge = slab - 1
        
    else:
        edge = slab 
        Is_odd = 0.0

    # MZ: you don't call this with edge -- your linspace needs to know the
    # correct N (not N-1) to get the right xedge values, so you call
    # with slab
    Is = simpson_int(CONST_A, CONST_B, func, slab) 
    Is = Is + Is_odd    
    I = exact_int(CONST_A, CONST_B)

    ea =   abs(Is - I)
    delta = abs(CONST_B -CONST_A)/slab
    
    printing(slab, I, Is, ea)
    
    ea_array.append(ea)
    delta_array.append(delta)

    
    


""" Plotting ea vs delta"""
pylab.loglog(delta_array, ea_array,  'bo')
pylab.xlabel('$\delta$ = (b-a)/N')
pylab.ylabel('Absolute Error')
pylab.grid(True)
pylab.savefig("simp.png")


""" Plotting ea vs N"""
pylab.clf()
pylab.cla()
pylab.loglog(CONST_N_SLABS, ea_array,  'go')
pylab.xlabel('N')
pylab.ylabel('Absolute Error')
pylab.grid(True)
pylab.savefig("simp2.png")






print "\n\nDone!"
