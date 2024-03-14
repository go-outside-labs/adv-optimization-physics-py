"""

    This module produces the outputs/plots.
    MIA STEIN, spring/2013
    
"""

import pylab 
import numpy
from matplotlib import pyplot

def plotPhaseSpace( b, aTheta, aOmega, t, power, k):
    pylab.clf()
    pylab.cla()

    label = str(b)
    
    pylab.subplot(221)  
    pylab.plot( aTheta, aOmega, color="m", lw = 2)
    pylab.xlabel(r"$\theta$ (radians) ", fontsize=10)
    pylab.ylabel('$\omega$ (radians/seconds)', fontsize=10)
    pylab.grid(True)
    
    pylab.subplot(222)
    pylab.plot(  t, aTheta, color="g", lw = 2)
    pylab.ylabel(r"$\theta$ (radians)", fontsize=10)
    pylab.xlabel('t (seconds)', fontsize=10)
    pylab.grid(True)

    pylab.subplot(223)
    pyplot.grid(True)
    pyplot.plot(k, power, color="c", lw = 2)
    pyplot.ylabel("|F(k)$|^{2}$", fontsize=10)
    pyplot.xlabel(r"$\nu_k$ ($s^{-1}$)", fontsize=10)
    
    pylab.subplot(224)
    pyplot.yscale('log')
    pyplot.plot(2.0*numpy.pi*k, power, color="b", lw = 1)
    pylab.xlim(0,6)
    pyplot.grid(True)
    pyplot.xlabel(r"$\nu_k$ ($s^{-1}$)", fontsize=10)
    pyplot.ylabel("log |F(k)$|^{2}$", fontsize=10)
    

    pylab.savefig("plots/b-%s_phase_space.png" % (label) )
    
    return 0



def plotDFT( b, bDFT, k, bDFT_inv, aTheta, t):
    pylab.clf()
    pylab.cla()

    label = str(b)

    ymax = max( bDFT.real )
    imax = numpy.where(bDFT.real == ymax )[0]
    xmax = k[imax]

    pylab.subplot(221)
    pylab.plot(t,  aTheta.real, color="g", lw = 2)
    pylab.ylabel( r"$\theta$ (t)", fontsize=10)
    pylab.xlabel('t (seconds)', fontsize=10)
    pylab.grid(True)
    
    pylab.subplot(222)
    pylab.annotate("Frequency has the most power", xy=(xmax, ymax), xycoords='data', xytext=(+10, +30), textcoords='offset points', fontsize=10, arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    pylab.plot(k, bDFT.real, color="r", lw = 2, label="Real F(k)") 
    pylab.plot(k, bDFT.imag, color="b", lw = 1, label="Imaginary F(k)") 
    leg = pylab.legend(loc=4,labelspacing=0.0005)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(0)
    pylab.ylabel("F(k)", fontsize=10)
    pylab.xlabel(r"$\nu_k$ ($s^{-1}$)", fontsize=10)
    pylab.grid(True)
    
    pylab.subplot(223)
    pylab.plot(k, abs(bDFT.real), color="r", lw = 2) 
    pylab.xlabel(r"$\nu_k$ ($s^{-1}$)", fontsize=10)
    pylab.ylabel("|F(k)|", fontsize=10)
    pylab.grid(True)
    
    pylab.subplot(224)
    pylab.plot(t, bDFT_inv.real,  color="y", lw = 2)
    pylab.ylabel("Inverse F(k)", fontsize=10)
    pylab.xlabel('t (seconds)', fontsize=10)
    pylab.grid(True)

    pylab.savefig("plots/b-%s_dft.png" % (label) )
    
    return 0
    
