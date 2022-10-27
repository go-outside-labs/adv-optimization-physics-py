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
theta0 = numpy.radians(10.0)

# period estimate
T = 2.0*math.pi*math.sqrt(L/g)*(1.0 + theta0**2/16.0)
print "period = ", T

tmax = 10.0*T
#dts = [0.5, 0.25, 0.125, 0.06125, 0.030625]
dts = [1, 0.1, 0.01, 0.001, 0.0001]

err = []
for dt in dts:
    HVVerlet = intVVerlet(theta0, dt, tmax, rhs)
    
    E = HVVerlet.energy()
    err.append(abs(E[-1]-E[0])/abs(E[0]))

    print dt, abs(E[-1]-E[0])/abs(E[-1])
    pylab.plot(HVVerlet.t, HVVerlet.theta, label="dt = %f" % (dt))


leg = pylab.legend(frameon=False)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')


pylab.savefig("Verlet-theta-dt.png")

pylab.clf()

pylab.scatter(dts, err)

pylab.plot(numpy.array(dts), err[0]*(dts[0]/numpy.array(dts))**-2, label="2nd order scaling")
pylab.plot(numpy.array(dts), err[0]*(dts[0]/numpy.array(dts))**-4, ls=":", label="4th order scaling")
pylab.plot(numpy.array(dts), err[0]*(dts[0]/numpy.array(dts))**-6, ls="--", label="6th order scaling")

leg = pylab.legend(loc=2, frameon=False)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')


ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.xlim(0.00001, 1.0)
pylab.ylim(1.e-12, 1.e-1)

pylab.savefig("Verlet-E-converge.png")





