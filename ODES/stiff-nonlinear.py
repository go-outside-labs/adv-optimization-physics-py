# example of a stiff non-linear system -- this is the example defined in
# the VODE source code
#
# y1' = -0.04 y1 + 1.e4 y2 y3
# y2' =  0.04 y1 - 1.e4 y2 y3 - 3.e7 y2**2
# y3' =                         3.e7 y2**2
#
# y1(0) = 1
# y2(0) = 0
# y3(0) = 0
#
# long term behavior: y1 -> 0, y2 -> 0, y3 -> 1
#
# here we solve the system twice: once taking a small step initially
# and then increasing by 10x until we reach tmax -- this lets us see
# the intermediate behavior.  The second call has VODE handle the
# timestep throughout the entire integration and does 0 to tmax in one
# call.
#
# M. Zingale (2013-02-23)

import numpy
import pylab
import scipy
from scipy.integrate import ode

def rhs(t, Y):
    """ RHS of the system -- using 0-based indexing """
    y1 = Y[0]
    y2 = Y[1]
    y3 = Y[2]

    dy1dt = -0.04*y1 + 1.e4*y2*y3
    dy2dt =  0.04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
    dy3dt =                         3.e7*y2**2

    return numpy.array([dy1dt, dy2dt, dy3dt])


def jac(t, Y):
    """ J_{i,j} = df_i/dy_j """

    y1 = Y[0]
    y2 = Y[1]
    y3 = Y[2]
    
    df1dy1 = -0.04
    df1dy2 = 1.e4*y3
    df1dy3 = 1.e4*y2

    df2dy1 = 0.04
    df2dy2 = -1.e4*y3 - 6.e7*y2
    df2dy3 = -1.e4*y2

    df3dy1 = 0.0
    df3dy2 = 6.e7*y2
    df3dy3 = 0.0

    return numpy.array([ [ df1dy1, df1dy2, df1dy3 ],
                         [ df2dy1, df2dy2, df2dy3 ],
                         [ df3dy1, df3dy2, df3dy3 ] ])


def VODEIntegrate(Y0, dt, tmax):
    """ integrate using the VODE method, start with a timestep of dt and
        increase by 10x after each call until we reach tmax.  This is 
        the behavior used in the DVODE Fortran source. """

    r = ode(rhs, jac).set_integrator("vode", method="bdf", 
                                     with_jacobian=True,
                                     atol=1.e-10, rtol=1.e-10,
                                     nsteps = 15000, order=5) #, min_step=dt)

    t = 0.0
    r.set_initial_value(Y0, t)

    tout = [t]
    y1out = [Y0[0]]
    y2out = [Y0[1]]
    y3out = [Y0[2]]

    while r.successful() and r.t < tmax:
        r.integrate(r.t+dt)

        tout.append(r.t)
        y1out.append(r.y[0])
        y2out.append(r.y[1])
        y3out.append(r.y[2])

        dt = 10.0*dt

    return numpy.array(tout), \
        numpy.array(y1out), numpy.array(y2out), numpy.array(y3out)



Y0 = numpy.array([1.0, 0.0, 0.0])
tmax = 4.e7

# first call -- give an initial small timestep and have our integrate
# routine increase it by 10x each call to VODE to get to tmax
dt_init = 4.e-6

t, y1, y2, y3 = VODEIntegrate(Y0, dt_init, tmax)


# second call -- give an initial timestep of tmax -- VODE will handle
# everything
t_single, y1_single, y2_single, y3_single = VODEIntegrate(Y0, tmax, tmax)

ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.plot(t, y1, label=r"$y_1$", color="r")
pylab.plot(t, y2, label=r"$y_2$", color="g")
pylab.plot(t, y3, label=r"$y_3$", color="b")

pylab.scatter(t_single, y1_single, color="r", 
              label="r$y_1$ from single VODE call")

pylab.scatter(t_single, y2_single, color="g",
              label="r$y_2$ from single VODE call")

pylab.scatter(t_single, y3_single, color="b",
              label="r$y_3$ from single VODE call")

pylab.legend(loc=8, frameon=False, fontsize="small")

pylab.xlabel("$t$")
pylab.ylabel(r"$Y(t)$")
pylab.title("VODE reference stiff nonlinear system")

pylab.savefig("stiff-nonlinear.png")


