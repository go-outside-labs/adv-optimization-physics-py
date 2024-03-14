""""

    This program performs time series analysis for the damped driven pendulum 
    MIA STEIN, spring/2013 (based on Mike Zingale's codes)

"""


from pendulumIntegrator import rhs, intRK4
from constants import b0, dt, tmax, bstep, path, theta0, q, omegad
from plotter import plotPhaseSpace, plotDFT
import numpy
import os



def main():
    
    # create folder for plots
    try: 
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    numpy.seterr(divide='ignore')

    
    # generate results for variations of b
    for i in range(10):

        """
            vary b forward
        """
        b = b0 + bstep*i
    
        # integrate omega and theta for this b
        integral = intRK4( b, dt, tmax, rhs, theta0, q, omegad)
        aOmega = integral.omega
        aTheta = integral.theta
        aT = integral.t
        
        # calculate the discrete fourier transform and the inverse
        bDFT = numpy.fft.fft(aTheta)
        bDFT_inv = numpy.fft.ifft(bDFT)
        
        # calculate the frequency space
        #k = numpy.fft.fftfreq(len(aTheta), dt)
        k = 1./aT
        
        # generate plots
        plotPhaseSpace( b, aTheta, aOmega, aT, calculatePowerSpectrum(aTheta), k)
        plotDFT(b, bDFT, k, bDFT_inv, aTheta, aT) 
    
 
        """
            vary b backwards
        """
        b = b0 - bstep*i
        
        # integrate omega and theta for this b
        integral = intRK4( b, dt, tmax, rhs, theta0, q, omegad)
        aOmega = integral.omega
        aTheta = integral.theta
        aT = integral.t
        
        # calculate the discrete fourier transform and the inverse
        bDFT = numpy.fft.fft(aTheta)
        bDFT_inv = numpy.fft.ifft(bDFT)
        
        # calculate the frequency space
        #k = numpy.fft.fftfreq(len(aTheta), dt)
        k = 1./aT
        
        # generate plots
        plotPhaseSpace( b, aTheta, aOmega, aT, calculatePowerSpectrum(aTheta), k)
        plotDFT(b, bDFT, k, bDFT_inv, aTheta, aT) 
   
  
    print "\nDone!"




def calculatePowerSpectrum(theta):
    s = numpy.fft.fft(theta)
    return numpy.real(s*numpy.conjugate(s))/len(theta)**2



if __name__ == "__main__":
    main()