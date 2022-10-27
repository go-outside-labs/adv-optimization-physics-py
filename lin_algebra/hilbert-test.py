import numpy
from gauss import *
from matmul import *
import sys

# tests of gaussian elimination using a Hilbert matrix as input.
# see "The Condition Number for A Matrix" by James Keesling


def Hilbert(n):
    """ return a Hilbert matrix, H_ij = (i + j - 1)^{-1} """

    H = numpy.zeros((n,n), dtype=numpy.float64)

    i = 1
    while (i <= n):
        j = 1
        while (j <= n):
            H[i-1,j-1] = 1.0/(i + j - 1.0)

            j += 1
        i += 1

    return H
    

for N in range(2,16):

    A = Hilbert(N)
    xorig = numpy.arange(N)
    b = mult_Ax(A, xorig)

    # this is our version from class
    # gaussElim changes A in place -- send a copy
    #x = gaussElim(A.copy(), b.copy())  

    # alternately, use the built-in solver 
    x = numpy.linalg.solve(A, b)

    err = numpy.max(numpy.abs(x - xorig))

    if N == 2:
        print "N, absolute error, condition number"

    print N, err, numpy.linalg.cond(A,p=1)

