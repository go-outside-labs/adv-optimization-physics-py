# solve the stiff linear system
#
# dY/dt = A Y
#
#           / -alpha   beta \
# with  A = |               |
#           \  alpha  -beta /
#
# this has the analytic solution for Y_B(0) = 0 of
#
# Y_B/Y_A = [ exp{(alpha+beta) t} - 1 ] / [ (beta/alpha) exp{(alpha+beta) t} + 1]
#
# this shows that a characteristic timescale is t ~ 1/(alpha + beta)
#
# and in the limit t -> infinity, Y_B/Y_A ~ alpha/beta
#
# This example came from Frank Timmes lecture notes
#
# With Backward-Euler, we have the linear system: 
#
# (I - dt A) Y^{n+1} = Y^n
#
# also do Runge-Kutta 4th order and VODE

import numpy
import pylab
from scipy.integrate import ode

alpha = 1.e9
beta = 1.e-5


def mkMatrix(dt):
    """ return I - dt A """
    return numpy.matrix([[1.0 + dt*alpha, -dt*beta     ],
                         [-dt*alpha,      1.0 + dt*beta]])

def rhs(t, Y):
    """ RHS of the system for explicit methods """
    return numpy.array([-alpha*Y[0] + beta*Y[1], alpha*Y[0] - beta*Y[1]])


def jac(t, Y):
    """ the Jacobian of the system (just the matrix A) """
    return numpy.array([ [-alpha, beta], [alpha, -beta] ])


def integrate(YA_0, YB_0, dt, tmax):
    """ perform a backward-Euler integration """

    YA = YA_0
    YB = YB_0

    tout = [0.0]
    yaout = [YA]
    ybout = [YB]

    t = 0.0
    while (t < tmax):
        
        # create the matrix
        J = mkMatrix(dt)

        b = numpy.array([YA, YB])

        # solve the linear system J x = b
        x = numpy.linalg.solve(J, b)

        YA = x[0]; YB = x[1]

        t += dt
        
        tout.append(t)
        yaout.append(YA)
        ybout.append(YB)

    return numpy.array(tout), numpy.array(yaout), numpy.array(ybout)


def trapIntegrate(YA_0, YB_0, dt, tmax):
    """ perform an implicit trapezoid integration:
        y^{n+1} = y^n + dt/2 [ f(y^n) + f(y^{n+1}) ] """

    YA = YA_0
    YB = YB_0

    tout = [0.0]
    yaout = [YA]
    ybout = [YB]

    t = 0.0
    while (t < tmax):
        
        # create the matrix -- using dt/2
        J = mkMatrix(dt/2)

        # create the explicit source term
        b = numpy.array([(1.0 - 0.5*dt*alpha)*YA + 0.5*dt*beta*YB, 
                         0.5*dt*alpha*YA + (1.0 - 0.5*dt*beta)*YB])

        # solve the linear system J x = b
        x = numpy.linalg.solve(J, b)

        YA = x[0]; YB = x[1]

        t += dt
        
        tout.append(t)
        yaout.append(YA)
        ybout.append(YB)

    return numpy.array(tout), numpy.array(yaout), numpy.array(ybout)


def RK4Integrate(YA_0, YB_0, dt, tmax):
    """ 4th-order explicit Runge-Kutta for comparison """

    YA = YA_0
    YB = YB_0

    tout = [0.0]
    yaout = [YA]
    ybout = [YB]

    t = 0.0
    while (t < tmax):
        
        k1 = rhs(t, numpy.array([YA, YB]) )
        k2 = rhs(t+0.5*dt, numpy.array([YA+0.5*dt*k1[0], YB+0.5*dt*k1[1]]) )
        k3 = rhs(t+0.5*dt, numpy.array([YA+0.5*dt*k2[0], YB+0.5*dt*k2[1]]) )
        k4 = rhs(t+dt, numpy.array([YA+dt*k3[0], YB+dt*k3[1]]) )

        YA += (dt/6.0)*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])
        YB += (dt/6.0)*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])

        t += dt
        
        tout.append(t)
        yaout.append(YA)
        ybout.append(YB)

    return numpy.array(tout), numpy.array(yaout), numpy.array(ybout)



def VODEIntegrate(YA_0, YB_0, dt, tmax):
    """ integrate using the VODE method, don't take a step smaller than dt """

    r = ode(rhs, jac).set_integrator("vode", method="bdf", 
                                     with_jacobian=True,
                                     atol=1.e-16, rtol=1.e-13,
                                     nsteps = 15000, order=5) #, min_step=dt)

    t = 0.0
    r.set_initial_value([YA_0, YB_0], t)

    tout = [t]
    yaout = [YA_0]
    ybout = [YB_0]

    while r.successful() and r.t < tmax:
        r.integrate(r.t+dt)

        tout.append(r.t)
        yaout.append(r.y[0])
        ybout.append(r.y[1])

    return numpy.array(tout), numpy.array(yaout), numpy.array(ybout)


#--------------------------------------------------------------------------
tref = 1.0/(alpha + beta)  # characteristic timescale
yequil = alpha/beta        # equilibrium YB/YA ratio

# R-K works well with a timestep < tref.  For dt > tref, it does horrible.
# Backward-Euler always remains stable.

#--------------------------------------------------------------------------
# first run, dt = 0.5 tref
dt = 0.5*tref
tmax = 100*tref

t, YA, YB = integrate(1.0, 0.0, dt, tmax)
t2, YA2, YB2 = trapIntegrate(1.0, 0.0, dt, tmax)
trk, YArk, YBrk = RK4Integrate(1.0, 0.0, dt, tmax)

pylab.scatter(t/tref, YB/YA/yequil, 
              marker="x", s=25, color="b", label="backward-Euler")

pylab.scatter(t2/tref, YB2/YA2/yequil, 
              marker="o", color="g", label="implicit trapezoid")

pylab.scatter(trk/tref, YBrk/YArk/yequil, 
              marker="*", color="r", label="explicit R-K 4")


pylab.plot(t/tref, (numpy.exp( (alpha+beta)*t) - 1.0)/ \
               ( (beta/alpha)*numpy.exp( (alpha+beta)*t) + 1.0)/yequil, 
           ls="-", color="k", label="analytic")

pylab.legend(loc=2, frameon=False, fontsize="x-small")


pylab.xlim(0.0, tmax/tref)

pylab.xlabel(r"$t \cdot (\alpha + \beta)$")
pylab.ylabel(r"$(Y_B/Y_A) / (\alpha/\beta)$")
pylab.title(r"$\Delta t = 0.5 / (\alpha + \beta)$")

ax = pylab.gca()
ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

pylab.savefig("stiff-linear-dt0.5tref.png")



#--------------------------------------------------------------------------
# first run, dt = 5 tref
pylab.clf()
dt = 5*tref
tmax = 500*tref

t, YA, YB = integrate(1.0, 0.0, dt, tmax)
t2, YA2, YB2 = trapIntegrate(1.0, 0.0, dt, tmax)
trk, YArk, YBrk = RK4Integrate(1.0, 0.0, dt, tmax)

tvode, YAvode, YBvode = VODEIntegrate(1.0, 0.0, dt, tmax)


pylab.scatter(t/tref, YB/YA/yequil, 
              marker="x", s=25, color="b", label="backward-Euler")

pylab.scatter(t2/tref, YB2/YA2/yequil, 
              marker="o", color="g", label="implicit trapezoid")

pylab.scatter(trk/tref, YBrk/YArk/yequil, 
              marker="*", color="r", label="explicit R-K 4")

pylab.plot(t/tref, (numpy.exp( (alpha+beta)*t) - 1.0)/ \
               ( (beta/alpha)*numpy.exp( (alpha+beta)*t) + 1.0)/yequil, 
           ls="-", color="k", label="analytic")

pylab.legend(loc=2, frameon=False, fontsize="x-small")


pylab.xlim(0.0, tmax/tref)

pylab.xlabel(r"$t \cdot (\alpha + \beta)$")
pylab.ylabel(r"$(Y_B/Y_A) / (\alpha/\beta)$")
pylab.title(r"$\Delta t = 5 / (\alpha + \beta)$")

ax = pylab.gca()
ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

pylab.savefig("stiff-linear-dt5tref.png")


pylab.scatter(tvode/tref, YBvode/YAvode/yequil, 
              marker="d", color="c", label="VODE (implicit B-D)")

pylab.legend(loc=2, frameon=False, fontsize="x-small")

pylab.savefig("stiff-linear-dt5tref-VODE.png")
