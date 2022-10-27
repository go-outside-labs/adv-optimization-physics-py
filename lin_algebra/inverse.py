# find the inverse of a matrix using Gaussian elimination

import numpy

def inverse(AInput):
    """ return the inverse of AInput """

    A = AInput.copy()

    N = len(A[:,0])

    # A is square, with each dimension of length N
    if not (A.shape[0] == A.shape[1]):
        print "ERROR: A should be square"
        return None

    # create an identity matrix
    I = numpy.identity(N)

    # allocation for the inverse
    Ainv = numpy.zeros((N,N), dtype=A.dtype)

    # find the scale factors for each row -- this is used when pivoting
    scales = numpy.max(numpy.abs(A), 1)

    # keep track of the number of times we swapped rows
    numRowSwap = 0

    # we are essentially doing Gaussian elimination, but with A Ainv = I
    # each column of I represents a separate righthand side to an Ax = b
    # linear system

    # main loop over rows
    for k in range(N):
        
        # find the pivot row based on the size of column k -- only consider
        # the rows beyond the current row
        rowMax = numpy.argmax(A[k:, k]/scales[k:]) 
        if (k > 0): rowMax += k  # we sliced A from k:, correct for total rows

        # swap the row with the largest scaled element in the current column
        # with the current row (pivot) -- do this with b too!
        if not rowMax == k:
            A[[k, rowMax],:] = A[[rowMax, k],:]
            I[[k, rowMax],:] = I[[rowMax, k],:]
            numRowSwap += 1

        # do the forward-elimination for all rows below the current
        for i in range(k+1, N):
            coeff = A[i,k]/A[k,k]

            for j in range(k+1, N):
                A[i,j] += -A[k,j]*coeff

            A[i,k] = 0.0
            I[i,:] += -coeff*I[k,:]
    
    

    # back-substitution -- once for each column in the I matrix
    
    for c in range(N):

        # last solution is easy
        Ainv[N-1,c] = I[N-1,c]/A[N-1,N-1]

        for i in reversed(range(N-1)):
            sum = I[i,c]
            for j in range(i+1,N):
                sum += -A[i,j]*Ainv[j,c]
            Ainv[i,c] = sum/A[i,i]


    # determinant
    det = numpy.prod(numpy.diagonal(A))*(-1.0)**numRowSwap
    
    return Ainv


# output: numpy.savetxt("test.out", a, fmt="%5.2f", delimiter="  ")
# convert -font Courier-New-Regular -pointsize 20 text:test.out test.png


A = numpy.array([ [4, 3, 4, 10], [2, -7, 3, 0], [-2, 11, 1, 3], [3, -4, 0, 2] ], dtype=numpy.float64)
Ainv = inverse(A)
print "A . Ainv = \n", numpy.dot(A, Ainv)
print" Ainv = \n", Ainv

print " "

A = numpy.array([ [0, 1, 1], [1, 1, 0], [1, 0, 1] ], dtype=numpy.float64)
Ainv = inverse(A)
print "A . Ainv = \n", numpy.dot(A, Ainv)

print " "

A = numpy.array([ [0, 0, 0, 4], 
                  [0, 0, 3, 0], 
                  [5, 6, 7, 8],
                  [0, 4, 3, 2] ], dtype=numpy.float64)
Ainv = inverse(A)
print "A . Ainv = \n", numpy.dot(A, Ainv)



