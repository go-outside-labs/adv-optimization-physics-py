# integrate the equations of motion of a pendulum, w/o the small angle
# approximation

import numpy
import pylab
import math

# global parameters 
g = 9.81     # gravitational acceleration [m/s]
L = 9.81     # length of pendulum [m]


class pendulumHistory:
    """ simple container to store the pendulum history """

    def __init__(self):
        self.t = None
        self.theta = None
        self.omega = None

    def energy(self):
        """ return the energy (per unit mass) """
        return 0.5*L**2*self.omega**2 - g*L*numpy.cos(self.theta)
        


def rhs(theta, omega):
    """ equations of motion for a pendulum
        dtheta/dt = omega
        domega/dt = - (g/L) sin theta """

    return omega, -(g/L)*numpy.sin(theta)


def intEuler(theta0, dt, tmax, rhs):
    """ integrate the equations of motion using Euler's method """
        
    
    # initial conditions
    t = 0.0
    theta = theta0
    omega = 0.0    # at the maximum angle, the angular velocity is 0


    # store the history for plotting
    tPoints = [t]
    thetaPoints = [theta]
    omegaPoints = [omega]

    while (t < tmax):

        # get the RHS
        thetadot, omegadot = rhs(theta, omega)

        # advance
        thetanew = theta + dt*thetadot
        omeganew = omega + dt*omegadot

        t += dt

        # store
        tPoints.append(t)
        thetaPoints.append(thetanew)
        omegaPoints.append(omeganew)

        # set for the next step
        theta = thetanew; omega = omeganew

    # return a pendulumHistory object with the trajectory
    H = pendulumHistory()
    H.t = numpy.array(tPoints)
    H.theta = numpy.array(thetaPoints)
    H.omega = numpy.array(omegaPoints)

    return H



def intEC(theta0, dt, tmax, rhs):
    """ integrate the equations of motion using Euler-Cromer """
        
    
    # initial conditions
    t = 0.0
    theta = theta0
    omega = 0.0    # at the maximum angle, the angular velocity is 0


    # store the history for plotting
    tPoints = [t]
    thetaPoints = [theta]
    omegaPoints = [omega]

    while (t < tmax):

        # get the RHS
        thetadot, omegadot = rhs(theta, omega)

        # advance
        omeganew = omega + dt*omegadot
        thetanew = theta + dt*omeganew


        t += dt

        # store
        tPoints.append(t)
        thetaPoints.append(thetanew)
        omegaPoints.append(omeganew)

        # set for the next step
        theta = thetanew; omega = omeganew

    # return a pendulumHistory object with the trajectory
    H = pendulumHistory()
    H.t = numpy.array(tPoints)
    H.theta = numpy.array(thetaPoints)
    H.omega = numpy.array(omegaPoints)

    return H


def intVVerlet(theta0, dt, tmax, rhs):
    """ integrate the equations of motion using Euler-Cromer """
        
    
    # initial conditions
    t = 0.0
    theta = theta0
    omega = 0.0    # at the maximum angle, the angular velocity is 0


    # store the history for plotting
    tPoints = [t]
    thetaPoints = [theta]
    omegaPoints = [omega]

    while (t < tmax):

        # get the RHS at time-level n
        thetadot, omegadot = rhs(theta, omega)

        thetanew = theta + dt*thetadot + 0.5*dt**2*omegadot

        # get the RHS with the updated theta -- omega doesn't matter
        # here, since we only need thetadot and omega doesn't affect
        # that.
        thetadot_np1, omegadot_np1 = rhs(thetanew, omega)

        omeganew = omega + 0.5*dt*(omegadot + omegadot_np1)

        t += dt

        # store
        tPoints.append(t)
        thetaPoints.append(thetanew)
        omegaPoints.append(omeganew)

        # set for the next step
        theta = thetanew; omega = omeganew

    # return a pendulumHistory object with the trajectory
    H = pendulumHistory()
    H.t = numpy.array(tPoints)
    H.theta = numpy.array(thetaPoints)
    H.omega = numpy.array(omegaPoints)

    return H




# 10 degree pendulum
theta0 = 10.0*math.pi/180.0
dt = 0.1
tmax = 30.0

HEuler = intEuler(theta0, dt, tmax, rhs)
HEC = intEC(theta0, dt, tmax, rhs)
HVVerlet = intVVerlet(theta0, dt, tmax, rhs)

pylab.plot(HEuler.t, HEuler.theta, label="Euler")
pylab.plot(HEC.t, HEC.theta, label="Euler-Cromer")
pylab.plot(HVVerlet.t, HVVerlet.theta, label="velocity Verlet")

pylab.xlabel("t")
pylab.ylabel(r"$\theta$(t)")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.savefig("pendulum-theta-10.png")


pylab.clf()

pylab.subplot(211)

pylab.plot(HEuler.t, HEuler.energy(), label="Euler")

pylab.xlabel("t")
pylab.ylabel("E(t)")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.subplot(212)

pylab.plot(HEC.t, HEC.energy(), label="Euler-Cromer")
pylab.plot(HVVerlet.t, HVVerlet.energy(), label="velocity Verlet")

pylab.xlabel("t")
pylab.ylabel("E(t)")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.tight_layout()

pylab.savefig("pendulum-energy-10.png")




# 100 degree pendulum
pylab.clf()

theta0 = 100.0*math.pi/180.0
dt = 0.1
tmax = 30.0

HEuler = intEuler(theta0, dt, tmax, rhs)
HEC = intEC(theta0, dt, tmax, rhs)
HVVerlet = intVVerlet(theta0, dt, tmax, rhs)

pylab.plot(HEuler.t, HEuler.theta, label="Euler")
pylab.plot(HEC.t, HEC.theta, label="Euler-Cromer")
pylab.plot(HVVerlet.t, HVVerlet.theta, label="velocity Verlet")

pylab.xlabel("t")
pylab.ylabel(r"$\theta$(t)")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.savefig("pendulum-theta-100.png")


pylab.clf()

pylab.subplot(211)
pylab.plot(HEuler.t, HEuler.energy(), label="Euler")

pylab.xlabel("t")
pylab.ylabel("E(t)")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.subplot(212)
pylab.plot(HEC.t, HEC.energy(), label="Euler-Cromer")
pylab.plot(HVVerlet.t, HVVerlet.energy(), label="velocity Verlet")

pylab.xlabel("t")
pylab.ylabel("E(t)")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.tight_layout()

pylab.savefig("pendulum-energy-100.png")

