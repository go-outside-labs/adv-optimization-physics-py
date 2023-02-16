# An example of general linear least squares fitting.  Here our basis
# functions can be non-linear, but we fit to a linear combination of
# the basis functions.  We consider simple polynomials (1, x, x**2, ...)
# and the Legendre polynomials


import numpy
import pylab


def polyBasis(M, x):
    """ the basis functions for the fit -- here x**n """

    j = numpy.arange(M)

    # simple polynomials
    return x**j


def legBasis(M, x):
    """ the basis functions for the fit -- here they are the Legendra
    polynomials """


    # note that the legendre polynomials are orthogonal in the interval 
    # [-1, 1] -- we need to convert our x to that range
    
    basis = []

    m = 0
    while (m < M):
        c = numpy.zeros(M)
        c[m] = 1.0

        basis.append(numpy.polynomial.legendre.legval(x, c))
        
        m += 1

    return numpy.array(basis)


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


def generalRegression(x, y, sigma, M, basis):
    """ here, M is the number of fitting parameters.  We will fit to
        a function that is linear in the a's, using the basis functions
        from basis() """


    N = len(x)

    # construct the design matrix -- A_{ij} = Y_j(x_i)/sigma_i -- this is 
    # N x M.  Each row corresponds to a single data point, x_i, y_i
    A = numpy.zeros((N, M), dtype=numpy.float64)

    i = 0
    while (i < N):
        A[i,:] = basis(M, x[i])/sigma[i]
        i += 1

    # construct the MxM matrix for the linear system, A^T A:
    ATA = numpy.dot(numpy.transpose(A), A)
    print "size of A^T A:", ATA.shape
    print "condition number of A^T A:", numpy.linalg.cond(ATA)

    # construct the RHS
    b = numpy.dot(numpy.transpose(A), y/sigma)

    # solve the system
    a = numpy.linalg.solve(ATA, b)

    # return the chisq
    chisq = 0
    i = 0
    while (i < N):
        chisq += (numpy.sum(a*basis(M, x[i])) - y[i])**2/sigma[i]**2
        i += 1

    chisq /= N-M

    return a, chisq


#-----------------------------------------------------------------------------
N = 40
x = numpy.linspace(0, 100.0, N)


#-----------------------------------------------------------------------------
# test 1 -- quadratic data with M = 3
    
print "\ntest 1: quadratic data with M = 3, simple poly\n "

pylab.clf()

# make up the experimental data with errors

sigma = 5.0*numpy.ones(N)

y = yExperiment2(2.0, 1.50, -0.02, sigma, x)

pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the regression with M = 3 (1, x, x^2)
M = 3
a, chisq = generalRegression(x, y, sigma, M, polyBasis)

print "a = ", a

pylab.plot(x, a[0] + a[1]*x + a[2]*x*x )

print "reduced chisq = ", chisq

pylab.savefig("general-regression-M3.png")
    

    

#-----------------------------------------------------------------------------
# test 2 -- quadratic data with M = 10
    
print "\ntest 2: quadratic data with M = 10, simple poly\n "

pylab.clf()

# same data as above
pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the regression with M = 10 (1, x, x^2, ...)
M = 10
a, chisq = generalRegression(x, y, sigma, M, polyBasis)

print "a = ", a

yfit = numpy.zeros((N), dtype=x.dtype)
i = 0
while (i < N):
    base = polyBasis(M, x[i])
    yfit[i] = numpy.sum(a*base)
    i += 1

pylab.plot(x, yfit)

print "reduced chisq = ", chisq

pylab.savefig("general-regression-M10.png")
    

    
#-----------------------------------------------------------------------------
# test 3 -- quadratic data with M = 3 with Legendre polynomial basis

print "\ntest 3: quadratic data with M = 10, Legendre poly\n "
    
pylab.clf()

# same data as above
pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the regression with M = 3 Legendra polynomials
M = 10
a, chisq = generalRegression(x, y, sigma, M, legBasis)

print "a = ", a

yfit = numpy.zeros((N), dtype=x.dtype)
i = 0
while (i < N):
    base = legBasis(M, x[i])
    yfit[i] = numpy.sum(a*base)
    i += 1

pylab.plot(x, yfit)

print "reduced chisq = ", chisq

pylab.savefig("general-regression-M10-leg.png")
    

#-----------------------------------------------------------------------------
# test 4 -- quadratic data in [-1, 1] with simple polynomials

print "\ntest 4: quadratic data on [-1,1] with M = 10, simple poly\n "

pylab.clf()

N = 40
x = numpy.linspace(-1.0, 1.0, N)
    
# make up the experimental data with errors

sigma = 1.0*numpy.ones(N)

y = yExperiment2(2.0, 1.50, -0.02, sigma, x)

pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the regression with M = 10 (1, x, x^2)
M = 10
a, chisq = generalRegression(x, y, sigma, M, polyBasis)

print "a = ", a

yfit = numpy.zeros((N), dtype=x.dtype)
i = 0
while (i < N):
    base = polyBasis(M, x[i])
    yfit[i] = numpy.sum(a*base)
    i += 1

pylab.plot(x, yfit)


print "reduced chisq = ", chisq

pylab.savefig("general-regression-M10-m1p1.png")



#-----------------------------------------------------------------------------
# test 5 -- quadratic data in [-1, 1] with Legendre polynomials

print "\ntest 5: quadratic data on [-1,1] with M = 10, Legendre poly\n "

pylab.clf()

pylab.scatter(x,y)
pylab.errorbar(x, y, yerr=sigma, fmt=None)


# do the regression with M = 10 (1, x, x^2)
M = 10
a, chisq = generalRegression(x, y, sigma, M, legBasis)

print "a = ", a

yfit = numpy.zeros((N), dtype=x.dtype)
i = 0
while (i < N):
    base = legBasis(M, x[i])
    yfit[i] = numpy.sum(a*base)
    i += 1

pylab.plot(x, yfit)


print "reduced chisq = ", chisq

pylab.savefig("general-regression-M10-m1p1-leg.png")

