""""

    This program analysis the condition number of Hilbert Matrices.
    MIA STEIN, spring/2013 (based on Mike Zingale's codes)

"""


import numpy as npy
from scipy import linalg 
from gaussElimination import gaussElim
from createHilbert import createHilbert
from calculatePrecision import calculatePrecision





def main():
    
    aCond = []
    aOlimit = []
    epsilon = calculatePrecision()
    
    
    for N in range(2,16):     
        H, x = createHilbert(N)
    
        # calculate the error
        b = npy.dot(H,x)
        xtilde = gaussElim(H, b) 
        error_calculated = abs(npy.linalg.norm(x - xtilde))
        
        # verifies when it reaches O(1)
        if (error_calculated >= 1.0):
            aOlimit.append([N, x, xtilde])
        
        # calculates the conditional number use the numpy API to get the conditional number
        cond_npy = npy.linalg.cond(H, p="fro")
        cond_formal = npy.linalg.norm(npy.dot(abs(H), abs(linalg.inv(H))))
        cond_comp = (npy.linalg.norm(x -xtilde)/npy.linalg.norm(x))/epsilon
        m = npy.log10(cond_npy)
        aCond.append([error_calculated,cond_formal, cond_comp, cond_npy, m])



    
    for N in range(2,16): 
        print "\n*** Hilbert Matrix N x N =", N,"x",N, ": ***" 
        print "Calculated error =",aCond[N-2][0]
        print "Condition number:\n    Formal definition (cond(A)=|A||A^-1|) =", aCond[N-2][1] 
        print "    Calculated with the machine precision =",  aCond[N-2][2] 
        print "    From the Numpy API =", aCond[N-2][3] 
        print "    Number of digit of accuracy lost in solving Hx=b, giving by log(cond(H)) is ", aCond[N-2][4]
        
        
        
        
    print "\n\n\nThe error for the N x N Hilbert matrix becomes O(1) when N =", aOlimit[0][0]
    print "as we can see when comparing \nx =", aOlimit[0][1], " to \nxtilde =", aOlimit[0][2]
    print "\n\nDone!"


if __name__ == "__main__":
    main()