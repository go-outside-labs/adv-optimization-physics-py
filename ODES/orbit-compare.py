# compare different ODE methods on the circular orbit problem


import numpy
import pylab
import math
from orbit import *

# circular orbit
o = orbit(1.0, 0.0)   # eccentricity = 0

# orbital period
P = o.keplerPeriod()

histEuler = o.intEuler(0.05, P)
histEulerCromer = o.intEulerCromer(0.05, P)
histRK2 = o.intRK2(0.05, P)
histRK4 = o.intRK4(0.05, P)

pylab.plot(histEuler.x, histEuler.y, label="Euler's method", color="k")

# mark the Sun
pylab.scatter([0],[0],s=250,marker=(5,1),color="k")
pylab.scatter([0],[0],s=200,marker=(5,1),color="y")


pylab.plot(histEulerCromer.x, histEulerCromer.y, label="Euler-Cromer method",
           color="r")
pylab.plot(histRK2.x, histRK2.y, label="R-K 2nd order",
           color="b")
pylab.plot(histRK4.x, histRK4.y, label="R-K 4th order",
           color="g")

leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.xlim(-2,2)
pylab.ylim(-2,2)

ax = pylab.gca()
ax.set_aspect("equal", "datalim")


pylab.savefig("orbit-compare.png")
