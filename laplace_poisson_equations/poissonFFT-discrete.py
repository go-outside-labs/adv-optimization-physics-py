# solve a 2-d Poisson equation by differentiating the discretized Poisson
# equation and then substituting in the inverse Fourier transform and solving
# for the amplitudes in Fourier space.
#
# This is the way that Garcia and NR do it.
#
# Note: we need a periodic problem for an FFT
#
# M. Zingale (2013-04-02)

import pylab
import numpy
import math


# the analytic solution
def true(x,y):
    pi = math.pi
    return numpy.sin(2.0*pi*x)**2*numpy.cos(4.0*pi*y) + \
        numpy.sin(4.0*pi*x)*numpy.cos(2.0*pi*y)**2


# the righthand side
def frhs(x,y):
    pi = math.pi
    return 8.0*pi**2*numpy.cos(4.0*pi*y)*(numpy.cos(4.0*pi*x) - 
                                          numpy.sin(4.0*pi*x)) - \
           16.0*pi**2*(numpy.sin(4.0*pi*x)*numpy.cos(2.0*pi*y)**2 + 
                       numpy.sin(2.0*pi*x)**2 * numpy.cos(4.0*pi*y))




# create the domain
xmin = 0.0
xmax = 1.0

ymin = 0.0
ymax = 1.0

Nx = 64
Ny = 64

dx = (xmax - xmin)/Nx
dy = (ymax - ymin)/Ny

x = (numpy.arange(Nx) + 0.5)*dx
y = (numpy.arange(Ny) + 0.5)*dy

x2d = numpy.repeat(x, Ny)
x2d.shape = (Nx, Ny)

y2d = numpy.repeat(y, Nx)
y2d.shape = (Ny, Nx)
y2d = numpy.transpose(y2d)


# create the RHS
f = frhs(x2d, y2d)
print f

# compatibility conditions require that the RHS sum to zero
print "sum of RHS: ", numpy.sum(f)

print x2d.shape

# FFT of RHS
F = numpy.fft.fft2(f)

# get the wavenumbers -- we need these to be physical, so multiply by Nx
kx = Nx*numpy.fft.fftfreq(Nx)
ky = Ny*numpy.fft.fftfreq(Ny)

# make 2-d arrays for the wavenumbers
kx2d = numpy.repeat(kx, Ny)
kx2d.shape = (Nx, Ny)

ky2d = numpy.repeat(ky, Nx)
ky2d.shape = (Ny, Nx)
ky2d = numpy.transpose(ky2d)

# here the FFT frequencies are in the order 0 ... N/2-1, -N/2, ...
# the 0 component is not a physical frequency, but rather it is the DC
# signal.  Don't mess with it, since we'll divide by zero
oldDC = F[0,0]
F = 0.5*F*dx*dx/(numpy.cos(2.0*math.pi*kx2d/Nx) + 
                 numpy.cos(2.0*math.pi*ky2d/Ny) - 2.0)# + 1.e-20)

F[0,0] = oldDC

# transform back to real space
fsolution = numpy.real(numpy.fft.ifft2(F))


#pylab.subplot(121)
#pylab.imshow(fsolution, origin="lower", interpolation="nearest")
#pylab.colorbar()

#pylab.subplot(122)
pylab.imshow(true(x2d,y2d))
pylab.colorbar()

#pylab.tight_layout()

pylab.savefig("poissonFFT.png")


# error
print Nx, Ny, numpy.sqrt(dx*dx*numpy.sum( ( (fsolution - true(x2d,y2d))**2).flat))



