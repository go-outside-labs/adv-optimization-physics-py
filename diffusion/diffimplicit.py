"""
solve the diffusion equation:

 phi_t = k phi_{xx} 

with a first-order (in time) implicit discretization

M. Zingale (2013-04-03)
"""

import numpy
from scipy import linalg
import sys
import pylab


def diffuseImplicit(gr, phi, k, dt):
    """ diffuse phi implicitly through timestep dt """

    phinew = gr.scratchArray()
    
    alpha = k*dt/gr.dx**2

    # create the RHS of the matrix
    R = phi[gr.ilo:gr.ihi+1] 
    
    # create the diagonal, d+1 and d-1 parts of the matrix
    d = (1.0 + 2.0*alpha)*numpy.ones(gr.nx)
    u = -alpha*numpy.ones(gr.nx)
    u[0] = 0.0

    l = -alpha*numpy.ones(gr.nx)
    l[gr.nx-1] = 0.0

    # set the boundary conditions by changing the matrix elements

    # homogeneous neumann
    d[0] = 1.0 + alpha
    d[gr.nx-1] = 1.0 + alpha


    # solve
    A = numpy.matrix([u,d,l])
    phinew[gr.ilo:gr.ihi+1] = linalg.solve_banded((1,1), A, R)

    return phinew


class grid:

    def __init__(self, nx, ng=1, xmin=0.0, xmax=1.0):
        """ grid class initialization """
        
        self.nx = nx
        self.ng = ng

        self.xmin = xmin
        self.xmax = xmax

        self.dx = (xmax - xmin)/nx
        self.x = (numpy.arange(nx+2*ng) + 0.5 - ng)*self.dx + xmin

        self.ilo = ng
        self.ihi = ng+nx-1

        # storage for the solution
        self.phi = numpy.zeros((nx+2*ng), dtype=numpy.float64)

    def fillBC(self):
        """ fill the Neumann BCs """

        # Neumann BCs
        self.phi[0:self.ilo]  = self.phi[self.ilo]
        self.phi[self.ihi+1:] = self.phi[self.ihi]


    def scratchArray(self):
        return numpy.zeros((2*self.ng+self.nx), dtype=numpy.float64)


    def phi_a(self, t, k, t0, phi1, phi2):
        """ analytic solution """

        xc = 0.5*(self.xmin + self.xmax)
        return (phi2 - phi1)*numpy.sqrt(t0/(t + t0)) * \
            numpy.exp(-0.25*(self.x-xc)**2/(k*(t + t0))) + phi1

    def norm(self, e):
        """ return the norm of quantity e which lives on the grid """
        if not len(e) == (2*self.ng + self.nx):
            return None

        return numpy.sqrt(self.dx*numpy.sum(e[self.ilo:self.ihi+1]**2))



def evolve(nx, k, t0, phi1, phi2, C, tmax):
    """ 
    the main evolution loop.  Evolve 
  
     phi_t = k phi_{xx} 

    from t = 0 to tmax
    """

    # create the grid
    gr = grid(nx, ng=1, xmax=1.0)

    # time info
    dt = C*0.5*gr.dx**2/k
    t = 0.0

    # initialize the data
    gr.phi[:] = gr.phi_a(0.0, k, t0, phi1, phi2)

    while (t < tmax):

        # make sure we end right at tmax
        if (t + dt > tmax):
            dt = tmax - t

        # diffuse for dt
        phinew = diffuseImplicit(gr, gr.phi, k, dt)

        gr.phi[:] = phinew[:]

        t += dt

    return gr


if __name__ == "__main__":

    # reference time
    t0 = 1.e-4

    # state coeffs
    phi1 = 1.0
    phi2 = 2.0

    k = 1.0


    #-------------------------------------------------------------------------
    # normal time

    tmax = 0.005

    nx = 128

    C = 0.8
    gr = evolve(nx, k, t0, phi1, phi2, C, tmax)

    pylab.plot(gr.x[gr.ilo:gr.ihi+1], gr.phi[gr.ilo:gr.ihi+1], color="r", label="C = 0.8")


    C = 2.0
    gr = evolve(nx, k, t0, phi1, phi2, C, tmax)

    pylab.plot(gr.x[gr.ilo:gr.ihi+1], gr.phi[gr.ilo:gr.ihi+1], color="g", label="C = 2.0")


    C = 10.0
    gr = evolve(nx, k, t0, phi1, phi2, C, tmax)

    pylab.plot(gr.x[gr.ilo:gr.ihi+1], gr.phi[gr.ilo:gr.ihi+1], color="b", label="C = 10.0")


    # analytic solution
    pylab.plot(gr.x[gr.ilo:gr.ihi+1], 
               gr.phi_a(tmax, k, t0, phi1, phi2)[gr.ilo:gr.ihi+1], 
               ls=":", color="0.5", label="analytic solution")


    pylab.legend(frameon=False)


    pylab.xlabel("$x$")
    pylab.ylabel(r"$\phi$")
    pylab.title("Backward-difference implicit diffusion, nx = %d, C = %3.2f, t = %5.2g" % (nx, C, tmax))

    pylab.savefig("diffimplicit.png")

    #-------------------------------------------------------------------------
    # early time

    pylab.clf()

    tmax = 0.0005

    nx = 128

    C = 0.8
    gr = evolve(nx, k, t0, phi1, phi2, C, tmax)

    pylab.plot(gr.x[gr.ilo:gr.ihi+1], gr.phi[gr.ilo:gr.ihi+1], color="r", label="C = 0.8")


    C = 2.0
    gr = evolve(nx, k, t0, phi1, phi2, C, tmax)

    pylab.plot(gr.x[gr.ilo:gr.ihi+1], gr.phi[gr.ilo:gr.ihi+1], color="g", label="C = 2.0")


    C = 10.0
    gr = evolve(nx, k, t0, phi1, phi2, C, tmax)

    pylab.plot(gr.x[gr.ilo:gr.ihi+1], gr.phi[gr.ilo:gr.ihi+1], color="b", label="C = 10.0")


    # analytic solution
    pylab.plot(gr.x[gr.ilo:gr.ihi+1], 
               gr.phi_a(tmax, k, t0, phi1, phi2)[gr.ilo:gr.ihi+1], 
               ls=":", color="0.5", label="analytic solution")


    pylab.legend(frameon=False)

    pylab.xlim(0.3,0.7)

    pylab.xlabel("$x$")
    pylab.ylabel(r"$\phi$")
    pylab.title("Backward-difference implicit diffusion, nx = %d, C = %3.2f, t = %5.2g" % (nx, C, tmax))

    pylab.savefig("diffimplicit-early.png")


