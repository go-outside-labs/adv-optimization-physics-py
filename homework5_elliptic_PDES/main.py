"""
    Any vector field U can be decomposed into a divergence free term
    Ud and the gradient of a scalar, phi:
        U = Ud + phi
    This program recovers a divergemce-free filed on a 2-d grid.
    
    MIA STEIN, based on M. Zingale's codes, spring 2013
"""


from part_a import doPartA
from part_b import doPartB, error
from part_c import doPartC
import os
import numpy


def main():

    """
        Do you want to make plots???
    """
    DO_PLOTS = 1


    """ create folder for plots """
    try: 
        os.makedirs("plots/")
    except OSError:
        if not os.path.isdir("plots/"):
            raise
    numpy.seterr(divide='ignore')
    
    
    """ [0,1]x[0,1]"""
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0
    
    
    """ setting the number of cells """
    nx = [32,64]
    ny = [32,64]
    ng = 1


    for i in range(len(nx)):
        print "Calculating for [%d,%d] cells..." %(nx[i], ny[i])  
        
        """ set the limits, grid limits   """
        ilo = ng
        ihix = ng + nx[i] - 1
        ihiy = ng + ny[i] - 1
    
        """ coordinates of centers """
        dx = (xmax - xmin)/nx[i]
        dy = (ymax - ymin)/ny[i]
        
        """
            a) Constructing the vector field
        """
        
        Ustar, Ustar_final, phi_true, Ud = doPartA(nx[i], ny[i], ng, dx, dy, DO_PLOTS)
        
        
        """
            b) Solving discrete Poisson equation using Gauss-Seidel relaxation
        """
        
        phi_num  = doPartB(Ustar, phi_true, nx[i], ny[i],  ng, dx, dy, DO_PLOTS)
        
        
        """
            c) Recovering the original and calculating the error
        """
        
        
        doPartC(Ustar_final, phi_num, Ud, nx[i], ny[i], xmin, xmax, ymin, ymax, DO_PLOTS)
        phi_true_final = phi_true[ilo:ihiy+1]
        print "The error is %.8f.\n" %(error(ilo, ihix, ihiy, dx, phi_num - phi_true_final))    
    
    
    print "\nDone!"




if __name__ == "__main__":
    main()
