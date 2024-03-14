#########################################################################
#            NUMERICAL DERIVATIVES                                        
#                                                                          
# This module explore compare the difference between 2-points, 3-points                        
# and 4-points derivatives, with the difference centered or one-sided.
#
#    Copied from Mike Zingale's, with addition of d2r and d2l.
#    (MIA STEIN, spring 2013)
#
#########################################################################

import math
import numpy
import pylab




""" Functions """

def func(x):
    """ the function to be plotted """
    f = numpy.sin(x)
    return f


def fprime(x):
    """ the analytic derivative of func(x) """
    fp = numpy.cos(x)
    return fp


def d1l(dx, fc, i):
    """ first-order, left-sided derivative at index i """
    D = (fc[i] - fc[i-1])/dx
    return D


def d1r(dx, fc, i):
    """ first-order, right-sided derivative at index i """
    D = (fc[i+1] - fc[i])/dx
    return D


def d2(dx, fc, i):
    """ second-order centered derivative at index i """
    D = fc[i+1] - fc[i-1]
    D = D/(2.0*dx)
    return D


def d2r(dx, fc, i):
    """ second-order right-sided derivative at index i """
    D = -fc[i] + 4.0*fc[i+1] - 3.0*fc[i+2] 
    D = D/(2.0*dx)
    return D


def d2l(dx, fc, i):
    """ second-order left-sided derivative at index i """
    D = -fc[i+2] + 4.0*fc[i+1] -3.0*fc[i] 
    D = D/(2.0*dx)
    return D


def d4(dx, fc, i):
    """ fourth-order centered derivative at index i """
    D = -fc[i+2] + 8.0*fc[i+1] - 8.0*fc[i-1] + fc[i-2]
    D = D/(12.0*dx)
    return D


def line(x, slope, x0, y0):
    return y0 + slope*(x - x0)




""" Constants """

xl = 0.0
xr = math.pi     
coarse_point_chosen = 3     # which for this griding is pi/3       




""" Main Program """

# fine and  grid (to show exact function and for differencing)
# linspace(start, stop, number_of_points)
fine = numpy.linspace(xl, xr, 500)      # 500 points, from 0 to pi
coarse = numpy.linspace(xl, xr, 10)     # 10 points, from 0 to pi


# define the analytic function for fine and for coarse griding, and the derivative for coarse
f = func(fine)
c = func(coarse)
analytic = fprime(coarse[coarse_point_chosen])


# equally spaced dx for coarse points
dx = coarse[1] - coarse[0]


# get the derivative approximations
Dl = d1l(dx, c, coarse_point_chosen)
Dr = d1r(dx, c, coarse_point_chosen)
D2 = d2(dx, c, coarse_point_chosen)
D2r = d2r(dx, c, coarse_point_chosen)
D2l = d2l(dx, c, coarse_point_chosen)
D4 = d4(dx, c, coarse_point_chosen)



""" Plotting """
# plot the fine and discrete gridded analytic function
pylab.plot(fine, f, color="0.5", lw=2)
pylab.scatter(coarse, c)


# function and point values
x0 = coarse[coarse_point_chosen]
y0 = func(x0)
xplot = numpy.linspace(coarse[1], coarse[5], 50)

pylab.plot(xplot, line(xplot, Dl, x0, y0), "b--", label="left-sided first-order approx")
pylab.plot(xplot, line(xplot, Dr, x0, y0), "bp", label="right-sided first-order approx")
pylab.plot(xplot, line(xplot, D2, x0, y0), "y", lw=3,  label="centered second-order approx")
pylab.plot(xplot, line(xplot, D2l, x0, y0), "r--", label="left-sized second-order approx")
pylab.plot(xplot, line(xplot, D2r, x0, y0), "rp", label="right-sized second-order approx")
pylab.plot(xplot, line(xplot, D4, x0, y0), "m*", lw=3, label="centered fourth-order approx")
#pylab.plot(xplot, line(xplot, analytic, x0, y0), color="0.5", ls=":", label="analytic")

print "analytic:          ", analytic
print "left-sided O(dx):  ", Dl
print "right-sided O(dx): ", Dr
print "centered O(dx**2): ", D2
print "left-sided O(dx**2): ", D2l
print "right-sided O(dx**2): ", D2r
print "centered O(dx**4): ", D4

leg = pylab.legend(loc=2,labelspacing=0.1)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(1)

# axes run through 0
# http://matplotlib.org/examples/pylab_examples/spine_placement_demo.html
ax = pylab.gca()
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['top'].set_color('none')
ax.spines['left'].set_smart_bounds(True)
ax.spines['bottom'].set_smart_bounds(True)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

pylab.xlabel("x (radians)")
pylab.ylabel("Sin(x)")
pylab.grid(True)
pylab.xlim(0,2.0)

pylab.savefig("fprime.png")



print "\nDone!"
