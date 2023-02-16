"""
solve a scalar diffusion-reaction equation:

phi_t = kappa phi_{xx} + (1/tau) R(phi)

using operator splitting, with implicit diffusion

"""

import numpy
from scipy import linalg
import sys
import pylab

def frhs(phi, tau):
    """ reaction ODE righthand side """
    return 0.25*phi*(1.0 - phi)/tau


def react(gr, phi, tau, dt):
    """ react phi through timestep dt """

    phinew = gr.scratchArray()

    # integrate with a single R-K 4th order step
    k1 = frhs(phi[gr.ilo:gr.ihi+1], tau)
    k2 = frhs(phi[gr.ilo:gr.ihi+1]+k1*0.5*dt, tau)
    k3 = frhs(phi[gr.ilo:gr.ihi+1]+k2*0.5*dt, tau)
    k4 = frhs(phi[gr.ilo:gr.ihi+1]+k3*dt, tau)

    phinew[gr.ilo:gr.ihi+1] = phi[gr.ilo:gr.ihi+1] + \
        (dt/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)

    return phinew


def diffuse(gr, phi, kappa, dt):
    """ diffuse phi implicitly (C-N) through timestep dt """

    phinew = gr.scratchArray()
    
    alpha = kappa*dt/gr.dx**2

    # create the RHS of the matrix
    R = phi[gr.ilo:gr.ihi+1] + \
        0.5*alpha*(    phi[gr.ilo-1:gr.ihi] - 
                   2.0*phi[gr.ilo  :gr.ihi+1] + 
                       phi[gr.ilo+1:gr.ihi+2])

    
    # create the diagonal, d+1 and d-1 parts of the matrix
    d = (1.0 + alpha)*numpy.ones(gr.nx)
    u = -0.5*alpha*numpy.ones(gr.nx)
    u[0] = 0.0

    l = -0.5*alpha*numpy.ones(gr.nx)
    l[gr.nx-1] = 0.0

    # set the boundary conditions by changing the matrix elements

    # homogeneous neumann
    d[0] = 1.0 + 0.5*alpha
    d[gr.nx-1] = 1.0 + 0.5*alpha

    # dirichlet
    #d[0] = 1.0 + 1.5*alpha
    #R[0] += alpha*0.0

    #d[gr.nx-1] = 1.0 + 1.5*alpha
    #R[gr.nx-1] += alpha*0.0

    # solve
    A = numpy.matrix([u,d,l])
    phinew[gr.ilo:gr.ihi+1] = linalg.solve_banded((1,1), A, R)

    return phinew


def estDt(gr, kappa, tau):
    """ estimate the timestep """

    # use the proported flame speed
    s = numpy.sqrt(kappa/tau)
    dt = gr.dx/s
    return dt


class grid:

    def __init__(self, nx, ng=1, xmin=0.0, xmax=1.0, vars=None):
        """ grid class initialization """
        
        self.nx = nx
        self.ng = ng

        self.xmin = xmin
        self.xmax = xmax

        self.dx = (xmax - xmin)/nx
        self.x = (numpy.arange(nx+2*ng) + 0.5 - ng)*self.dx + xmin

        self.ilo = ng
        self.ihi = ng+nx-1

        self.data = {}

        for v in vars:
            self.data[v] = numpy.zeros((2*ng+nx), dtype=numpy.float64)


    def fillBC(self, var):

        if not var in self.data.keys():
            sys.exit("invalid variable")

        vp = self.data[var]

        # Neumann BCs
        vp[0:self.ilo+1] = vp[self.ilo]
        vp[self.ihi+1:] = vp[self.ihi]


    def scratchArray(self):
        return numpy.zeros((2*self.ng+self.nx), dtype=numpy.float64)


    def initialize(self, L, delta):
        """ initial conditions """

        # this looks like a smooth tophat delta is the thickness of the
        # transition and L is the length of the region set to 1

        phi = self.data["phi"]
        xc = 0.5*(self.xmin + self.xmax)
        x1 = xc - L/2
        x2 = xc + L/2

        phi[:] = 0.0
        phi[:] = 0.5*(1.0 + numpy.tanh((self.x-x1)/delta)) * \
                 0.5*(1.0 - numpy.tanh((self.x-x2)/delta))


def interpolate(x, phi, phipt):
    """ find the x position corresponding to phipt """

    idx = (numpy.where(phi >= 0.5))[0][0]
    xs   = numpy.array([x[idx-1],   x[idx],   x[idx+1]])
    phis = numpy.array([phi[idx-1], phi[idx], phi[idx+1]])

    xpos = 0.0

    m = 0
    while (m < len(phis)):
        # create Lagrange basis polynomial for point m
        l = None
        n = 0
        while (n < len(phis)):
            if n == m:
                n += 1
                continue
                
            if l == None:
                l = (phipt - phis[n])/(phis[m] - phis[n])
            else:
                l *= (phipt - phis[n])/(phis[m] - phis[n])

            n += 1
            
        xpos += xs[m]*l
        m += 1

    return xpos


def evolve(nx, kappa, tau, tmax, dovis=0, returnInit=0):
    """ 
    the main evolution loop.  Evolve 
  
     phi_t = kappa phi_{xx} + (1/tau) R(phi)

    from t = 0 to tmax
    """

    # create the grid
    gr = grid(nx, ng=1, xmin = 0.0, xmax=100.0,
              vars=["phi", "phi1", "phi2"])

    # pointers to the data at various stages
    phi  = gr.data["phi"]
    phi1 = gr.data["phi1"]
    phi2 = gr.data["phi2"]

    # initialize
    gr.initialize(20, 1)

    phiInit = phi.copy()

    # runtime plotting
    if dovis == 1: pylab.ion()
    
    t = 0.0
    while (t < tmax):

        dt = estDt(gr, kappa, tau)

        if (t + dt > tmax):
            dt = tmax - t

        # react for dt/2
        phi1[:] = react(gr, phi, tau, dt/2)
        gr.fillBC("phi1")

        # diffuse for dt
        phi2[:] = diffuse(gr, phi1, kappa, dt)
        gr.fillBC("phi2")

        # react for dt/2 -- this is the updated solution
        phi[:] = react(gr, phi2, tau, dt/2)
        gr.fillBC("phi")

        t += dt

        if dovis == 1:
            pylab.clf()
            pylab.plot(gr.x, phi)
            pylab.xlim(gr.xmin,gr.xmax)
            pylab.ylim(0.0,1.0)
            pylab.draw()

    print t

    if returnInit==1:
        return phi, gr.x, phiInit
    else:
        return phi, gr.x
    

# phi is a reaction progress variable, so phi lies between 0 and 1

kappa = 0.1
tau = 1.0

tmax1 = 60.0

nx = 256

phi1, x1 = evolve(nx, kappa, tau, tmax1)

tmax2 = 80.0

phi2, x2 = evolve(nx, kappa, tau, tmax2)

pylab.plot(x1, phi1)
pylab.plot(x2, phi2, ls=":")
pylab.savefig("flame.png")


# estimate the speed -- interpolate to x corresponding to where phi > 0.2
xpos1 = interpolate(x1, phi1, 0.2)
xpos2 = interpolate(x2, phi2, 0.2)

print (xpos1 - xpos2)/(tmax1 - tmax2), numpy.sqrt(kappa/tau)

# estimate the speed -- interpolate to x corresponding to where phi > 0.5
xpos1 = interpolate(x1, phi1, 0.5)
xpos2 = interpolate(x2, phi2, 0.5)

print (xpos1 - xpos2)/(tmax1 - tmax2), numpy.sqrt(kappa/tau)

# estimate the speed -- interpolate to x corresponding to where phi > 0.8
xpos1 = interpolate(x1, phi1, 0.8)
xpos2 = interpolate(x2, phi2, 0.8)

print (xpos1 - xpos2)/(tmax1 - tmax2), numpy.sqrt(kappa/tau)


# make a pretty plot
pylab.clf()

dt = 8.0
for i in range(0,10):
    tend = (i+1)*dt
    p, x, phi0 = evolve(nx, kappa, tau, tend, returnInit=1)

    c = 1.0 - (0.1 + i*0.1)
    pylab.plot(x, p, color=`c`)


pylab.plot(x, phi0, ls=":", color="0.5")

pylab.xlabel("$x$")
pylab.ylabel("$\phi$")
pylab.title(r"Diffusion-Reaction, $nx = %d, \, \kappa = %3.2f, \, \tau = %3.2f, \, \delta t = %3.2f$ (between lines)" % (nx ,kappa, tau, dt), fontsize=12)

pylab.tight_layout()

pylab.xlim(0,100)
pylab.savefig("flame_seq.png")

