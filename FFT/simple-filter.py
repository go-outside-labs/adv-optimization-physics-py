import pylab
import numpy
import math

# simple example of removing a frequency component from a
# two-frequency sine wave signal

# Since our input data is real, the negative frequency components
# don't include any new information, and are not interesting to us.
# The rfft routines understand this, and rfft takes n real points and
# returns n/2+1 complex output points.  The corresponding inverse
# knows this, and acts accordingly.
#
# these are the routines we want for real valued data

# note that the scipy version of rfft returns that data differently

# M. Zingale (2013-02-28)


def twoFreqSine(npts):

    # a pure sine with no phase shift will result in pure imaginary
    # signal
    f_0 = 0.2
    f_1 = 0.5

    xmax = 10.0/f_0

    # we call with endpoint=False -- if we include the endpoint, then for
    # a periodic function, the first and last point are identical -- this
    # shows up as a signal in the FFT.
    xx = numpy.linspace(0.0, xmax, npts, endpoint=False)

    # input frequency
    f_0 = 0.2

    f = 0.5*(numpy.sin(2.0*math.pi*f_0*xx) + numpy.sin(2.0*math.pi*f_1*xx))
    
    return xx, f




npts = 256

xx, f = twoFreqSine(npts)

# Forward transform: f(x) -> F(k)

# normalization factor: the 2 here comes from the fact that we neglect
# the negative portion of frequency space because our input function
# is real
norm = 2.0/npts   
fk = norm*numpy.fft.rfft(f)

ofk_r = fk.real.copy()
ofk_i = fk.imag.copy()


# the fftfreq returns the postive and negative (and 0) frequencies
# the newer versions of numpy (>=1.8) have an rfftfreq() function
# that really does what we want.
k = numpy.fft.fftfreq(len(xx))[range(0,npts/2+1)]

# the last element is negative, because of the symmetry, but should
# be positive (see http://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.fft.rfftfreq.html)
k[-1] *= -1


# since we don't include the endpoint in xx, to normalize things, we need
# max(xx) + dx to get the true range
kfreq = k*npts/(max(xx) + xx[1])


# filter all frequencies > 0.4
fk[kfreq > 0.4] = 0.0


# element 0 of fk is the DC component
fk_r = fk.real
fk_i = fk.imag


# Inverse transform: F(k) -> f(x)    
fkinv = numpy.fft.irfft(fk/norm)

# PLOT
pylab.rc("font", size=9)


pylab.subplot(411)

pylab.plot(xx, f)
pylab.xlabel("x")
pylab.ylabel("f(x)")


pylab.subplot(412)

pylab.plot(kfreq, ofk_r)
pylab.plot(kfreq, ofk_i, ls=":")
pylab.xlabel(r"$\nu_k$")
pylab.ylabel("original F(k)")

#pylab.xlim(0,1)

pylab.subplot(413)

pylab.plot(kfreq, fk_r)
pylab.plot(kfreq, fk_i, ls=":")
pylab.xlabel(r"$\nu_k$")
pylab.ylabel("filtered F(k)")


pylab.subplot(414)

pylab.plot(xx, fkinv.real)
pylab.xlabel("x")
pylab.ylabel(r"inverse filtered F(k)")

pylab.tight_layout()

f = pylab.gcf()
f.set_size_inches(9.0,8.0)

pylab.savefig("simple-filter.png", bbox_inches='tight', pad_inches=0.33)


