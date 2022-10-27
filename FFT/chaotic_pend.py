# a damped, driven pendulum, from Pang, Ch. 4

import math
import numpy
import pylab

class ddpend:

    def __init__(self, theta0, omega0, q, b, omega_d):
        self.theta0 = theta0  # initial angular displacement
        self.omega0 = omega0  # initial angular displacement
        
        self.q = q             # damping parameter
        self.b = b             # forcing amplitude
        self.omega_d = omega_d # driving frequency

        self.t = None
        self.theta = None
        self.omega = None

    def rhs(self, t, theta, omega):
        """ return the RHS (thetadot(t), omegadot(t)) """

        thetadot = omega
        omegadot = -self.q*omega - math.sin(theta) \
            + self.b*math.cos(self.omega_d*t)

        return thetadot, omegadot

    def intRK4(self, dt, tmax):
        """ integrate the system using 4th-order Runge-Kutta """

        # initial condition
        t = 0.0
        theta = self.theta0
        omega = self.omega0

        # store the solution
        tHist = [t]
        thetaHist = [theta]
        omegaHist = [omega]

        # integrate
        while (t < tmax):

            thetadot1, omegadot1 = self.rhs(t, theta, omega)

            thetadot2, omegadot2 = self.rhs(t+0.5*dt, 
                                            theta+0.5*dt*thetadot1, 
                                            omega+0.5*dt*omegadot1)

            thetadot3, omegadot3 = self.rhs(t+0.5*dt, 
                                            theta+0.5*dt*thetadot2, 
                                            omega+0.5*dt*omegadot2)
            
            thetadot4, omegadot4 = self.rhs(t+dt, 
                                            theta+dt*thetadot3, 
                                            omega+dt*omegadot3)

            theta += (dt/6.0)*(thetadot1 + 2.0*thetadot2 + 2.0*thetadot3 + thetadot4)
            omega += (dt/6.0)*(omegadot1 + 2.0*omegadot2 + 2.0*omegadot3 + omegadot4)

            t += dt
    
            tHist.append(t)
            thetaHist.append(theta)
            omegaHist.append(omega)


        self.t = numpy.array(tHist)
        self.theta = numpy.array(thetaHist)
        self.omega = numpy.array(omegaHist)

    def restrictTheta(self):
        """ convert theta in place to be restricted to lie between -pi
            and pi.  This is done in a periodic fashion, with theta' =
            theta +/- 2n pi """
        
        # shift everything by pi, then restrict to lie between [0,
        # 2pi], then shift back by pi

        self.theta += math.pi

        n = 0
        while (n < len(self.theta)):
            self.theta[n] += - 2.0*math.pi*math.floor(self.theta[n]/(2.0*math.pi)) 
            n += 1
            
        self.theta -= math.pi


    def powerSpectrum(self):
        """ return the power spectrum of theta.  For the frequency
            component, return it in terms of omega """

        # power spectrum
        N = len(self.t)

        F = (2.0/N)*numpy.fft.rfft(self.theta)

        k = numpy.fft.fftfreq(N)[range(0,N/2+1)]
        if N % 2 == 0:
            k[-1] *= -1
            
        kfreq = 2.0*math.pi*k*N/max(self.t)

        return kfreq, F


#-----------------------------------------------------------------------------
# normal (non-damped, non-driven) pendulum

# Note, without damping or driving, all the power should be at the
# natural oscillation period of the pendulum.  For a small amplitude,
# with L = g, then the period is T = 2 pi, and the frequency is nu_k =
# 1/(2 pi).  We plot things in terms of the angular frequency, omega_k
# = 2 pi nu_k, so all the power will be at omega_k = 1

# For a large amplitude perturbation, the period will be longer, so
# the power will be at an omega_k < 1

q = 0.0
b = 0.0
omega_d = 2./3.

T_d = 2.0*math.pi/omega_d
dt = T_d/200.0

# these conditons give a large amplitude perturbation
#theta0 = 0.0
#omega0 = 2.0

# these conditions give a small amplitude, so the power for the undamped,
# non-driven pendulum should be at omega_k = 1
theta0 = 0.1
omega0 = 0.0

p0 = ddpend(theta0, omega0, q, b, omega_d)
p0.intRK4(dt, 100.0*T_d)

pylab.subplot(211)

pylab.plot(p0.theta, p0.omega)

pylab.xlabel(r"$\theta$")
pylab.ylabel(r"$\omega$")

# power spectrum
omega_k, F = p0.powerSpectrum()

pylab.subplot(212)

pylab.plot(omega_k, numpy.abs(F)**2)

pylab.xlim(0.,2.)
#pylab.ylim(1.e-4,1.0)

ax = pylab.gca()
#ax.set_yscale('log')

pylab.xlabel(r"$\omega_k$")
pylab.ylabel(r"power spectrum")

pylab.tight_layout()

pylab.savefig("pend_nodamping.png")


#-----------------------------------------------------------------------------
# non-chaotic pendulum
q = 0.5
b = 0.9
omega_d = 2./3.

T_d = 2.0*math.pi/omega_d
dt = T_d/200.0

theta0 = 0.0
omega0 = 2.0

p1 = ddpend(theta0, omega0, q, b, omega_d)
p1.intRK4(dt, 100.0*T_d)

pylab.clf()

pylab.subplot(211)

pylab.plot(p1.theta, p1.omega)

pylab.xlabel(r"$\theta$")
pylab.ylabel(r"$\omega$")

# power spectrum
omega_k, F = p1.powerSpectrum()

pylab.subplot(212)

pylab.plot(omega_k, numpy.abs(F)**2)
pylab.plot([omega_d, omega_d], [1.e-10,2.0*max(numpy.abs(F)**2)], ls=":")

pylab.xlim(0.,1.)
pylab.ylim(1.e-4,1.0)

ax = pylab.gca()
ax.set_yscale('log')

pylab.xlabel(r"$\omega_k$")
pylab.ylabel(r"power spectrum")

pylab.tight_layout()

pylab.savefig("pend_q0.5_b0.9_om0.666.png")


#-----------------------------------------------------------------------------
# Chaotic pendulum
q = 0.5
bmin = 0.9
db = 0.05
N = 20

B = numpy.arange(N)*db + bmin

omega_d = 2./3.

T_d = 2.0*math.pi/omega_d
dt = T_d/200.0

theta0 = 0.0
omega0 = 2.0


for b in B:

    p2 = ddpend(theta0, omega0, q, b, omega_d)
    p2.intRK4(dt, 500.0*T_d)

    p2.restrictTheta()

    pylab.clf()

    pylab.subplot(211)

    pylab.plot(p2.theta, p2.omega)

    pylab.title(r"$q = %3.2f, \, \omega_d = %4.3f, \, b = %3.2f$" % (q, omega_d, b))

    pylab.xlabel(r"$\theta$")
    pylab.ylabel(r"$\omega$")

    # power spectrum
    omega_k, F = p2.powerSpectrum()

    pylab.subplot(212)

    pylab.plot(omega_k, numpy.abs(F)**2)
    pylab.plot([omega_d, omega_d], [1.e-10,2.0*max(numpy.abs(F)**2)], ls=":")

    pylab.xlim(0.,6.)
    pylab.ylim(1.e-4,1.0)

    ax = pylab.gca()
    ax.set_yscale('log')

    pylab.xlabel(r"$\omega_k$")
    pylab.ylabel(r"power spectrum")

    pylab.tight_layout()

    pylab.savefig("pend_q0.5_b%3.2f_om0.666.png" % (b))



