# two dimensional FFT example

# following http://matlabgeeks.com/tips-tutorials/how-to-do-a-2-d-fourier-transform-in-matlab/

import matplotlib.cm
import pylab
import numpy
import math

# get the figure
f = pylab.imread("luna_bw.png")

# two dimension FFT -- F is complex
F = numpy.fft.fft2(f)

print f.shape
print F.shape

# find the mag and phase -- shift to put 0 wavenumber at the center
F_mag = numpy.abs(numpy.fft.fftshift(F))
F_phase = numpy.angle(numpy.fft.fftshift(F))


pylab.rc("font", size=9)

pylab.subplot(131)
pylab.imshow(f, cmap=matplotlib.cm.Greys_r)
pylab.title("original image")

pylab.subplot(132)
pylab.imshow(numpy.log(F_mag))
pylab.title("|F(k)|")

pylab.subplot(133)
pylab.imshow(F_phase)
pylab.title("phase of F(k)")

f = pylab.gcf()
f.set_size_inches(10.0,6.0)

pylab.savefig("fft2d.png", bbox_inches="tight")

#-------------------------------------------------------------------------------
# scramble phase

pylab.clf()

Fnew_phase = 2.0*math.pi*numpy.random.rand(F_phase.shape[0], F_phase.shape[1])

# back to the complex representation
Fnew = F_mag*numpy.exp(1j*Fnew_phase)

fnew = numpy.fft.ifft2(numpy.fft.ifftshift(Fnew))

pylab.imshow(numpy.real(fnew), cmap=matplotlib.cm.Greys_r)
pylab.title(r"F$^{-1}$(F(k)) with scrampled phases")
pylab.savefig("fft2d_phasescamble.png", bbox_inches="tight")


#-------------------------------------------------------------------------------
# scramble amplitude

pylab.clf()

Fnew_mag = numpy.max(F_mag)*numpy.random.rand(F_mag.shape[0], F_mag.shape[1])

# back to the complex representation
Fnew = Fnew_mag*numpy.exp(1j*F_phase)

fnew = numpy.fft.ifft2(numpy.fft.ifftshift(Fnew))

pylab.imshow(numpy.real(fnew), cmap=matplotlib.cm.Greys_r)
pylab.title(r"F$^{-1}$(F(k)) with scrampled amplitudes")
pylab.savefig("fft2d_magscamble.png", bbox_inches="tight")


#-------------------------------------------------------------------------------
# filter out high frequencies

pylab.clf()

# http://glowingpython.blogspot.com/2011/08/fourier-transforms-and-image-filtering.html

F_orig = numpy.fft.fftshift(F)

P = numpy.zeros(F.shape, dtype=numpy.complex128)

frac = 0.25
rad = frac*int(min(F.shape)/2)

ic = F.shape[0]/2
jc = F.shape[1]/2

for i in range(F.shape[0]):
    for j in range(F.shape[1]):

        if math.sqrt( (i-ic)**2 + (j-jc)**2) < rad:
            P[i,j] = F_orig[i,j]


f_filtered = numpy.real(numpy.fft.ifft2(numpy.fft.ifftshift(P)))

pylab.subplot(131)
pylab.imshow(numpy.log(numpy.abs(F_orig)))
pylab.title("original |F(k)|")

pylab.subplot(132)
pylab.imshow(numpy.log(numpy.abs(P)))
pylab.title("filtered |F(k)|")

pylab.subplot(133)
pylab.imshow(f_filtered, cmap=matplotlib.cm.Greys_r)
pylab.title(r"filtered F$^{-1}$(F(k))")

f = pylab.gcf()
f.set_size_inches(10.0,6.0)

pylab.savefig("fft2d_filtered.png", bbox_inches="tight")


