# compute the inverse of a matrix via Newton iteration.  This is based
# on problem 4.12 in Garcia.  See also
# http://amca01.wordpress.com/2010/08/18/the-nsh-method-for-matrix-inversion/
# (which cites Gunther Schultz 1933 and Harold Hotelling 1943)

# also Wikipedia for invertible matrix notes:
#
#   Newton inverse is useful for inverting similar matrixes (use
#   previous inverse as starting guess).
#
# Note: this is very sensitive to the initial guess.  It will diverge
# quickly if the initial guess gives products that are all > 1


import numpy

tol = 1.e-12

def iter_inverse(A, Ainv0):
    """ here A is the matrix and Ainv0 is an initial guess to the
        inverse """

    Ainv = Ainv0.copy()

    err = 1.e10

    iter = 0
    while (err > tol):
        Ainv_new = 2.0*Ainv - numpy.dot(Ainv, numpy.dot(A, Ainv))

        err = numpy.max(numpy.abs(Ainv - Ainv_new))

        Ainv = Ainv_new.copy()

        iter += 1


    print "number of iterations = ", iter

    return Ainv


# some attempts
A = numpy.array([ [ 4,  3, 4, 10], 
                  [ 2, -7, 3,  0], 
                  [-2, 11, 1,  3],
                  [ 3, -4, 0,  2] ], dtype=numpy.float64)


# identity
print "calling with Ainv0 = I"
Ainv = iter_inverse(A, numpy.eye(4))

print numpy.dot(A, Ainv)

print  " "


# transpose scaled by maximum element **2
print "calling with Ainv0 = A^T/max(A)**2"
Ainv = iter_inverse(A, numpy.transpose(A)/numpy.max(numpy.abs(A))**2)

print numpy.dot(A, Ainv)

print  " "


# diagonal of 1/ A's diagonal
Ainv = numpy.diagflat(1.0/numpy.diag(A))

print "calling with diag(Ainv0) = 1.0/diag(A)"
Ainv = iter_inverse(A, numpy.transpose(A)/numpy.max(numpy.abs(A))**2)
print numpy.dot(A, Ainv)

print  " "


# matrix with all elements = 1/max(A)
Ainv = numpy.ones(A.shape)/numpy.max(numpy.abs(A))

print "calling with Ainv_ij = 1.0/max(|A|)"
Ainv = iter_inverse(A, numpy.transpose(A)/numpy.max(numpy.abs(A))**2)
print numpy.dot(A, Ainv)

print  " "


