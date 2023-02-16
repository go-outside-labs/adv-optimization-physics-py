# simple example of linear regression.  We make up some "experimental"
# data by calling a random number generator to perturb a line and then
# fit to that line

import numpy
import pylab

def yExperiment(a1, a2, sigma, x):
    """ return the experimental data in a linear + random fashion a1
        is the intercept, a2 is the slope, and sigma is the error """

    N = len(x)

    # randn gives samples from the "standard normal" distribution
    r = numpy.random.randn(N)
    
    y = a1 + a2*x + sigma*r

    return y


def yExperiment2(a1, a2, a3, sigma, x):
    """ return the experimental data in a quadratic + random fashion,
        with a1, a2, a3 the coefficients of the quadratic and sigma is
        the error.  This will be poorly matched to a linear fit for 
        a3 != 0 """

    N = len(x)

    # randn gives samples from the "standard normal" distribution
    r = numpy.random.randn(N)
    
    y = a1 + a2*x + a3*x*x + sigma*r

    return y


def linearRegression(x, y, sigma, sigmaa=0):

    N = len(x)

    S = numpy.sum(1.0/sigma**2)
    
    xi_1 = numpy.sum(x/sigma**2)
    xi_2 = numpy.sum(x*x/sigma**2)

    eta = numpy.sum(y/sigma**2)
    mu = numpy.sum(x*y/sigma**2)

    a2 = (S*mu - xi_1*eta)/(xi_2*S - xi_1**2)
    a1 = (eta*xi_2 - mu*xi_1)/(xi_2*S - xi_1**2)

    chisq = numpy.sum( (a1 + a2*x - y)**2/sigma**2)
    chisq /= N-2

    if not sigmaa:
        return a1, a2, chisq

    else:
        sigmaFit = numpy.array([xi_2/(S*xi_2 - xi_1**2), 
                                S/(S*xi_2 - xi_1**2)])
        sigmaFit = numpy.sqrt(sigmaFit)

        return a1, a2, chisq, sigmaFit


#-----------------------------------------------------------------------------
N = 40
x = numpy.linspace(0.0, 100.0, N)

#-----------------------------------------------------------------------------
# test 1 -- linear data

print "Linear fit to linear data\n"

# make up the experimental data with errors
sigma = 25.0*numpy.ones(N)

y = yExperiment(10.0, 3.0, sigma, x)

pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the linear regression
a1, a2, chisq, sigmaFit = linearRegression(x, y, sigma, sigmaa=1)

pylab.plot(x, a1 + a2*x)

print "reduced chisq = ", chisq
print " a1 = %f +/- %f\n a2 = %f +/- %f\n" % (a1, sigmaFit[0], a2, sigmaFit[1])

pylab.savefig("linear-regression.png")


#-----------------------------------------------------------------------------
# test 2 -- quadratic data
    
pylab.clf()

print "Linear fit to quadratic data\n"

# make up the experimental data with errors
sigma = 5.0*numpy.ones(N)

y = yExperiment2(2.0, 1.50, -0.02, sigma, x)

pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the linear regression
a1, a2, chisq = linearRegression(x, y, sigma)

pylab.plot(x, a1 + a2*x)

print "reduced chisq = ", chisq

pylab.savefig("linear-regression-quad.png")
