# integrate an elliptical orbit with fixed timestep 4th-order Runge-Kutta

import numpy
import pylab
import math
from orbit_adaptive import *

# circular orbit
o = orbit(1, 0.95)    # a, e

# period
P = o.keplerPeriod()

err = -1  # no adaptive stepping -- use the input dt always
dt = 0.0005
histRK4 = o.intRK4(dt, err, P)

print "number of steps = ", len(histRK4.t)-1

# plot the orbit
pylab.plot(histRK4.x, histRK4.y, label="4th order RK", color="b")
#pylab.scatter(histRK4.x, histRK4.y, marker="x", color="b")

# mark the Sun
pylab.scatter([0],[0],s=250,marker=(20,1),color="k")
pylab.scatter([0],[0],s=200,marker=(20,1),color="y")

pylab.xlim(-1.5,1.5)
pylab.ylim(-2.5,0.5)

ax = pylab.gca()
ax.set_aspect("equal", "datalim")

pylab.savefig("orbit-rk4-nonadaptive.png")


# plot the energy
pylab.clf()

E = histRK4.energy()
pylab.plot(histRK4.t, E/E[0])

pylab.ylabel(r"$E(t)/E_0$")
pylab.xlabel(r"t")

ax = pylab.gca()
ax.ticklabel_format(useOffset=False)

pylab.savefig("energy-rk4-nonadaptive.png")


# plot the timesteps
pylab.clf()

dt = histRK4.t[1:-1] - histRK4.t[0:-2]
tmid = 0.5*(histRK4.t[1:-1] + histRK4.t[0:-2])
pylab.scatter(tmid, dt)


ax = pylab.gca()
ax.set_yscale('log')

pylab.ylim(1.e-5,1.e-1)

pylab.savefig("dt-rk4-nonadaptive.png")
 

