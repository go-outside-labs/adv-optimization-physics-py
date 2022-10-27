#!/usr/bin/env python

"""

Solve the diffusion equation with Crank-Nicolson time-discretization
using MG.

M. Zingale (2013-04-07)

"""
#from io import *
import numpy
import patch1d
import multigrid
import pylab


# global parameters for the analytic solution
phi1 = 1.0
phi2 = 2.0

# t0 is the initial time in the initial conditions -- the smaller you
# make this, the more the initial conditions approximate a delta
# function.  Careful though, too small, and you will not resolve this
# on your grid.
t0 = 1.e-4

# diffusion coefficient
k = 1.0

# the analytic solution
def phi_a(myg, t):
    xc = 0.5*(myg.xmin + myg.xmax)
    
    return (phi2 - phi1)*numpy.sqrt(t0/(t + t0)) * \
            numpy.exp(-0.25*(myg.x-xc)**2/(k*(t + t0))) + phi1


# the L2 error norm
def error(myg, r):

    # L2 norm of elements in r, multiplied by dx to
    # normalize
    return numpy.sqrt(myg.dx*numpy.sum((r[myg.ilo:myg.ihi+1]**2)))



def lap(gr, phi):
    """ compute the Laplacian of phi """

    lapphi = gr.scratchArray()

    ib = gr.ilo
    ie = gr.ihi

    lapphi[ib:ie+1] = (phi[ib-1:ie] - 2.0*phi[ib:ie+1] + phi[ib+1:ie+2])/gr.dx**2

    return lapphi



def evolve(nx, C, tmax, xmax=None):

    xmin = 0.0
    if xmax == None: xmax = 1.0

    # create a dummy patch to store some info in the same way the MG
    # solver will
    myGrid = patch1d.grid1d(nx, ng=1,
                            xmin=xmin, xmax=xmax)


    
    # initialize the data
    phi = myGrid.scratchArray()

    # initial solution -- this fills the GC too
    phi[:] = phi_a(myGrid, 0.0)
    
    # time info
    dt = C*0.5*myGrid.dx**2/k            
    t = 0.0

    # evolve
    while (t < tmax):

        if (t + dt > tmax):
            dt = tmax - t


        # create the multigrid object
        a = multigrid.ccMG1d(nx, xmin=xmin, xmax=xmax, 
                             alpha = 1.0, beta = 0.5*dt*k,
                             xlBCtype="neumann", xrBCtype="neumann",
                             verbose=0)
        
        # initialize the RHS
        a.initRHS(phi + 0.5*dt*k*lap(a.solnGrid, phi))
    
        # initialize the solution to 0
        a.initZeros()

        # solve to a relative tolerance of 1.e-11
        a.solve(rtol=1.e-11)

        # get the solution 
        v = a.getSolution()

        # store the new solution
        phi[:] = v[:]

        t += dt

    return a.solnGrid, phi


#-----------------------------------------------------------------------------
# test with various CFL

nx = 128
tmax = 5.e-3
xmax = 1.0

pylab.clf()

C = 0.8
grA, phiA = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grA.x[grA.ilo:grA.ihi+1], phiA[grA.ilo:grA.ihi+1], label = "C = 0.8")

C = 2.0
grB, phiB = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grB.x[grB.ilo:grB.ihi+1], phiB[grB.ilo:grB.ihi+1], label = "C = 2.0")

C = 10.0
grC, phiC = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grC.x[grC.ilo:grC.ihi+1], phiC[grC.ilo:grC.ihi+1], label = "C = 10.0")


pylab.plot(grA.x[grA.ilo:grA.ihi+1], 
           phi_a(grA, tmax)[grA.ilo:grA.ihi+1], 
           ls=":", color="0.5", label="analytic solution", lw=2)

pylab.legend(frameon=False)

pylab.xlim(0.0,1.0)

pylab.xlabel("$x$")
pylab.ylabel(r"$\phi$")
pylab.title("C-N implicit diffusion via MG, nx = %d, t = %5.2g" % (nx, tmax))

pylab.xlabel("x")
pylab.ylabel("$\phi$")

pylab.savefig("diffMG.png")


#-----------------------------------------------------------------------------
# test with various CFL -- earlier in time

nx = 128
tmax = 5.e-4
xmax = 1.0

pylab.clf()

C = 0.8
grA, phiA = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grA.x[grA.ilo:grA.ihi+1], phiA[grA.ilo:grA.ihi+1], label = "C = 0.8")

C = 2.0
grB, phiB = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grB.x[grB.ilo:grB.ihi+1], phiB[grB.ilo:grB.ihi+1], label = "C = 2.0")

C = 10.0
grC, phiC = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grC.x[grC.ilo:grC.ihi+1], phiC[grC.ilo:grC.ihi+1], label = "C = 10.0")


pylab.plot(grA.x[grA.ilo:grA.ihi+1], 
           phi_a(grA, tmax)[grA.ilo:grA.ihi+1], 
           ls=":", color="0.5", label="analytic solution", lw=2)

pylab.legend(frameon=False)

pylab.xlim(0.3, 0.7)

pylab.xlabel("$x$")
pylab.ylabel(r"$\phi$")
pylab.title("C-N implicit diffusion via MG, nx = %d, t = %5.2g" % (nx, tmax))

pylab.xlabel("x")
pylab.ylabel("$\phi$")

pylab.savefig("diffMG-early.png")


#-----------------------------------------------------------------------------
# test used for convergence

nx = 128
tmax = 5.e-3
xmax = 1.0

pylab.clf()

C = 0.8
grA, phiA = evolve(nx, C, tmax, xmax=xmax)

pylab.plot(grA.x[grA.ilo:grA.ihi+1], phiA[grA.ilo:grA.ihi+1], label = "C = 0.8", color="r")

pylab.plot(grA.x[grA.ilo:grA.ihi+1], 
           phi_a(grA, tmax)[grA.ilo:grA.ihi+1], 
           ls=":", color="0.5", label="analytic solution", lw=2)

pylab.xlim(0.0,1.0)

pylab.xlabel("$x$")
pylab.ylabel(r"$\phi$")
pylab.title("C-N implicit diffusion via MG, nx = %d, t = %5.2g" % (nx, tmax))

pylab.xlabel("x")
pylab.ylabel("$\phi$")

pylab.savefig("diffMGtest.png")


#-----------------------------------------------------------------------------
# convergence

pylab.clf()

# a characteristic timescale for diffusion if L^2/k
tmax = 5.e-3


N = [16, 32, 64, 128, 256, 512]

# CFL number
C = 0.8

# explicit dt for finest resolution
dt_finest = (1.0/N[-1])**2/k
N_finest = N[-1]

err = []

for nx in N:

    print nx

    # the present C-N discretization
    g, phi = evolve(nx, C, tmax)

    phi_analytic = phi_a(g, tmax)

    err.append(error(g, phi - phi_analytic))


pylab.clf()

N = numpy.array(N, dtype=numpy.float64)
err = numpy.array(err)

print "err = ", err

pylab.scatter(N, err, color="r", label="C-N implicit diffusion")
pylab.plot(N, err[len(N)-1]*(N[len(N)-1]/N)**2, color="k", label="$\mathcal{O}(\Delta x^2)$")

ax = pylab.gca()
ax.set_xscale('log')
ax.set_yscale('log')

pylab.xlabel(r"$N$")
pylab.ylabel(r"L2 norm of absolute error")
pylab.title("Convergence of C-N MG Diffusion, C = %3.2f" % (C))

pylab.legend(frameon=False, fontsize="small")

pylab.savefig("diffMG-converge.png")


