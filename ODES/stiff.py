# simple stiff ODE equation with analytic solution
#
# y' = -1.e3 (y - exp(-t) ) - exp(-t) 
# y(0) = 0
#
# analytic solution is
#
# y(t) = exp(-t) - exp(-1e3 t)
#
# this example comes from Byrne & Hindemarsh (1986)
#
# M. Zingale (2013-02-22)

import math
import numpy
import pylab

def analytic(t):
    return numpy.exp(-t) - numpy.exp(-1.e3*t)


def rhs(t, y):
    return -1.e3*(y - numpy.exp(-t)) - numpy.exp(-t)


def rk4(y0, dt, tmax):
    
    tsol = [0.0]
    ysol = [y0]

    t = 0.0
    y = y0

    while (t < tmax):
        ydot1 = rhs(t, y)
        ydot2 = rhs(t+0.5*dt, y+0.5*dt*ydot1)
        ydot3 = rhs(t+0.5*dt, y+0.5*dt*ydot2)
        ydot4 = rhs(t+dt, y+dt*ydot3)

        y += (dt/6.0)*(ydot1 + 2.0*ydot2 + 2.0*ydot3 + ydot4)

        t += dt

        tsol.append(t)
        ysol.append(y)

    return numpy.array(tsol), numpy.array(ysol)


def backwardEuler(y0, dt, tmax):

    tsol = [0.0]
    ysol = [y0]

    t = 0.0
    y = y0

    while (t < tmax):
        # an implicit discretication: y^{n+1} - y^n = dt ydot^{n+1}
        # and then solve analytically for y^{n+1}:
        ynew = (y + 1.e3*dt*numpy.exp(-t) - dt*numpy.exp(-t))/ \
               (1.0 + 1.e3*dt)

        y = ynew
        t += dt

        tsol.append(t)
        ysol.append(y)

    return numpy.array(tsol), numpy.array(ysol)

    


y0 = 0.0

# plot analytic solution
tt = numpy.linspace(0.0, 1.0, 1000)

pylab.plot(tt, analytic(tt), label="analytic solution")


# plot RK4 solution -- we need a step of 1.e-3 or smaller to get this
# right, because that is the fastest timescale.  If you go above 2.e-3,
# it blows up severely
tRK4, yRK4 = rk4(y0, 1.e-3, 1.0)

pylab.plot(tRK4, yRK4, label=r"R-K 4, $\tau = 10^{-3}$")

pylab.xlim(0.0, 0.1)

pylab.xlabel("t")
pylab.ylabel("y")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.savefig("stiff-rk4-dt-1e-3.png")

# now do dt = 2.5e-3
pylab.clf()

pylab.plot(tt, analytic(tt), label="analytic solution")

tRK4, yRK4 = rk4(y0, 2.5e-3, 1.0)

pylab.plot(tRK4, yRK4, label=r"R-K 4, $\tau = 2.5\times 10^{-3}$")

pylab.xlim(0.0, 0.1)

pylab.xlabel("t")
pylab.ylabel("y")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.savefig("stiff-rk4-dt-2.5e-3.png")


# now do dt = 5.e-3
pylab.clf()

pylab.plot(tt, analytic(tt), label="analytic solution")

tRK4, yRK4 = rk4(y0, 5.0e-3, 1.0)

pylab.plot(tRK4, yRK4, label=r"R-K 4, $\tau = 5\times 10^{-3}$")

pylab.xlim(0.0, 0.1)

pylab.xlabel("t")
pylab.ylabel("y")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

ax = pylab.gca()
ax.xaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))
ax.yaxis.set_major_formatter(pylab.ScalarFormatter(useMathText=True))

pylab.ylim(-1.e22, 0.0)

#for i in range(len(tRK4)):
#    print tRK4[i], yRK4[i]


pylab.savefig("stiff-rk4-dt-5e-3.png")


# implicit solution
pylab.clf()

pylab.plot(tt, analytic(tt), label="analytic solution")

tBE, yBE = backwardEuler(y0, 5.e-3, 1.0)

pylab.plot(tBE, yBE, label=r"backward-Euler, $\tau = 5\times 10^{-3}$")

tBE, yBE = backwardEuler(y0, 1.e-2, 1.0)

pylab.plot(tBE, yBE, label=r"backward-Euler, $\tau = 10^{-2}$")

pylab.xlim(0.0, 0.1)

pylab.xlabel("t")
pylab.ylabel("y")

leg = pylab.legend(loc=4)
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.savefig("stiff-be.png")

