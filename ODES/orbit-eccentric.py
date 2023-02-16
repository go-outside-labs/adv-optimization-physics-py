# compare Euler's method and the Euler-Cromer method for eccentric orbits


import numpy
import pylab
import math
from orbit import *

# circular orbit
o = orbit(1.0, 0.6)  # eccentricity = 0.6

# period
P = o.keplerPeriod()

histEuler = o.intEuler(0.0125, P)
histEC = o.intEulerCromer(0.0125, P)
#histRK2 = o.intRK2(0.0125, P)


pylab.plot(histEuler.x, histEuler.y, label="Euler's method", color="k")

# mark the Sun
pylab.scatter([0],[0],s=250,marker=(20,1),color="k")
pylab.scatter([0],[0],s=200,marker=(20,1),color="y")

# draw a vertical line that the semi-major axis should fall on
yy = numpy.linspace(-2.0, 2.0, 100)
pylab.plot(0.0*yy, yy, ls=":", color="0.5")


pylab.plot(histEC.x, histEC.y, label="Euler-Cromer method", color="b")



leg = pylab.legend()
ltext = leg.get_texts()
pylab.setp(ltext, fontsize='small')
leg.draw_frame(0)

pylab.xlim(-2,2)
pylab.ylim(-2,2)

ax = pylab.gca()
ax.set_aspect("equal", "datalim")


pylab.savefig("orbit-eccentric.png")
