"""

    This module finds the inverse of a matrix using Newton's method.
    MIA STEIN spring 2013.

"""


import numpy as npy


def inverseNewton(AInput):
    
    CONST_TOLERANCE = 10**-11
    A = AInput.copy()
    Asize = len(A[:,:])
    
    # verifies if A is square
    if not (A.shape[0] == A.shape[1]):
        print "ERROR: A should be square"
        return None


    # calculates initial guess and tolerance
    I = npy.identity(Asize)
    
    """ based on  Victor Pan and Robert Schreiber """
    X = npy.transpose(A)/(npy.linalg.norm(A, ord=1)*npy.linalg.norm(A, ord=npy.inf))       


    # iteration 
    dx = 1.0
    i = 0
    while (abs(dx) > abs(CONST_TOLERANCE)):
        X = npy.dot(X,2*I)-npy.dot(X,npy.dot(A,X))
        A_dx = npy.linalg.inv(A)
        dx = npy.linalg.norm(X-A_dx)
        i += 1
        
        
    return X, dx,i
