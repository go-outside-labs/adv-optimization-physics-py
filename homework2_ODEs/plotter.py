########################################################################
# This program illustrate many approaches to solve ODEs. In this case
#    we use the simple pendulum as the system. 
#
# This module implements all the classes and functions for the problem.
#
#
#    1) We compare the solutions for theta = 10o and 100o.
#    2) We make a plot of both theta(t) vs t and E vs t.
#
#    MIA STEIN, spring/2013 
#    (based on Mike Zingale's codes)
#
########################################################################

import numpy as np
import pylab
import math
from classes import *


# global parameters and constants
CONST_l = 10.0
CONST_g = 10.0
CONST_m = 1.0
CONST_dt = [0.005, 0.01, 0.05, 0.1, 0.5]
CONST_Tmax = 10.0
CONST_theta0 = [10.0, 100.0]


# calculating for the two values of theta0
theta10 = pendulumTrajectory(CONST_theta0[0], CONST_l, CONST_g, CONST_m)  
theta100 = pendulumTrajectory(CONST_theta0[1], CONST_l, CONST_g, CONST_m)  

T10 = theta10.pendulumPeriod()
T100 = theta10.pendulumPeriod()
E10 = theta10.pendulumEnergy()
E100 = theta100.pendulumEnergy()

# printing results
print "The approx Period of the simple pendulum for theta0 = 10 degrees is %f s-1." % (T10)
print "The approx Period of the simple pendulum for theta0 = 100 degrees is %f s-1." % (T100)
print "The Absolute Total Energy of the simple pendulum for theta0 = 10 degrees is %f joules." % (E10)
print "TheAbsolute Total Energy of the simple pendulum for theta0 = 100 degrees is %f joules." % (E100)


# iterating for many time steps
for i in range(len(CONST_dt)):
    
    # calculates for euler and for euler-cromer
    histEuler10 = theta10.intEuler(CONST_dt[i], CONST_Tmax)
    histEC10 = theta10.intEulerCromer(CONST_dt[i], CONST_Tmax)
    histEuler100 = theta100.intEuler(CONST_dt[i], CONST_Tmax)
    histEC100 = theta100.intEulerCromer(CONST_dt[i], CONST_Tmax)


    # plotting theta vs t theta0 = 10
    time = histEC10.t
    ana = CONST_theta0[0]*np.cos((2*math.pi/T10)*time)
                    
    pylab.clf()
    pylab.cla()
    
    pylab.plot(histEuler10.t, histEuler10.theta, label="Euler's method", color="m", lw = 2)
    pylab.plot(histEC10.t, histEC10.theta, label="Euler-Cromer's method", color="b", lw = 2)
    pylab.plot(time, ana, 'g*', label="Analytic Approximation")
   
    pylab.xlabel('Time (s)')
    pylab.ylabel('$\Theta$ (degrees)')
    pylab.title("Simple Pendulum: $\Theta_0$ = $10^o$, $\Delta$t = %s" % (str(CONST_dt[i])))   
    pylab.grid(True)
    leg = pylab.legend(loc=3,labelspacing=0.005)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(1)
    ax = pylab.gca()
    pylab.savefig("plots/theta_vs_t_for_theta0-10_dt-%s.png" % (str(CONST_dt[i])))



    # plotting theta vs t, theta0 = 100
    time = histEC100.t
    ana = CONST_theta0[0]*np.cos((2*math.pi/T100)*time)
                    
    pylab.clf()
    pylab.cla()
    
    pylab.plot(histEuler100.t, histEuler100.theta, label="Euler's method", color="m", lw = 2)
    pylab.plot(histEC100.t, histEC100.theta, label="Euler-Cromer's method", color="b", lw = 2)
    pylab.plot(time, ana, 'g*', label="Analytic Approximation")
   
    pylab.xlabel('Time (s)')
    pylab.ylabel('$\Theta$ (degrees)')
    pylab.title("Simple Pendulum: $\Theta_0$ = $100^o$, $\Delta$t = %s" % (str(CONST_dt[i])))   
    pylab.grid(True)
    leg = pylab.legend(loc=3,labelspacing=0.005)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(1)
    ax = pylab.gca()
    pylab.savefig("plots/theta_vs_t_for_theta0-100_dt-%s.png" % (str(CONST_dt[i])))


   
    # plotting E(t) vs t, for theta0 = 10
    time = histEC10.t
    fun_en = []
    for j in range (len(time)):
                    fun_en.append(E10)
    pylab.clf()
    pylab.cla()
    pylab.plot(histEuler10.t, histEuler10.energy, label="Euler's method", color="r", lw=2)
    pylab.plot(histEC10.t, histEC10.energy, label="Euler-Cromer's method", color="b", lw =2)
    pylab.plot(time, fun_en, label="Approximation from initial values", color="g", lw = 2)
    leg = pylab.legend()
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(0)
    pylab.grid(True)
    pylab.xlabel('Time (s)')
    pylab.ylabel('Absolute Energy (J)')
    pylab.title("Simple Pendulum: $\Theta_0$ = $10^o$, $\Delta$t = %s" % (str(CONST_dt[i])))   
    pylab.grid(True)
    leg = pylab.legend(loc=0,labelspacing=0.005)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(1)
    ax = pylab.gca()
    pylab.savefig("plots/energy_vs_t_for_theta0-10_dt-%s.png" % (str(CONST_dt[i])))
    

    # plotting E(t) vs t, for theta0 = 100
    time = histEC100.t
    fun_en = []
    for j in range (len(time)):
                    fun_en.append(E100)
    pylab.clf()
    pylab.cla()
    pylab.plot(histEuler100.t, histEuler100.energy, label="Euler's method", color="r", lw=2)
    pylab.plot(histEC100.t, histEC100.energy, label="Euler-Cromer's method", color="b", lw =2)
    pylab.plot(time, fun_en, label="Approximation from initial values", color="g", lw = 2)
    leg = pylab.legend()
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(0)
    pylab.grid(True)
    pylab.xlabel('Time (s)')
    pylab.ylabel('Absolute Energy (J)')
    pylab.title("Simple Pendulum: $\Theta_0$ = $100^o$, $\Delta$t = %s" % (str(CONST_dt[i])))   
    pylab.grid(True)
    leg = pylab.legend(loc=0,labelspacing=0.005)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(1)
    ax = pylab.gca()
    pylab.savefig("plots/energy_vs_t_for_theta0-100_dt-%s.png" % (str(CONST_dt[i])))
    
    
    
# plotting E(t) vs t, for theta0 = 10 for the verlet method
CONST_dt = [0.0001, 0.001, 0.01, 0.1]
pylab.clf()
pylab.cla()
for j in range(len(CONST_dt)):
    histVer10 = theta10.intVerlet(CONST_dt[j], CONST_Tmax)
    pylab.loglog(histVer10.t, histVer10.dE, 'ro')
pylab.setp(ltext, fontsize='small')
pylab.grid(True)
pylab.ylim(10,0.001)
pylab.xlabel('Time Step ($\Delta t$)')
pylab.ylabel('Absolute Error in Energy')
pylab.title("Verlet Method for the Simple Pendulum, $\Theta_0$ = $10^o$")
ax = pylab.gca()
pylab.savefig("plots/verlet_energy_vs_t_for_theta0-10.png" )


# plotting E(t) vs t, for theta0 = 100 for the verlet method
CONST_dt = [0.0001, 0.001, 0.01, 0.1]
pylab.clf()
pylab.cla()
for k in range(len(CONST_dt)):
    histVer100 = theta100.intVerlet(CONST_dt[k], CONST_Tmax)
    pylab.loglog(histVer100.t, histVer100.dE, 'bo')
    print CONST_dt[k], histVer100.dE
pylab.setp(ltext, fontsize='small')
pylab.grid(True)
pylab.xlabel('Time Step ($\Delta t$)')
pylab.ylabel('Absolute Error in Energy')
pylab.title("Verlet Method for the Simple Pendulum, $\Theta_0$ = $100^o$")
ax = pylab.gca()
pylab.savefig("plots/verlet_energy_vs_t_for_theta0-100.png" )
    



print "\nDone!"

