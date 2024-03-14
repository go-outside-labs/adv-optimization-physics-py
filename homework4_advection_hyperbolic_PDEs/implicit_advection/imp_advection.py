"""
    implicit * upwind *finite-difference * implementation the * linear advection *
    
    MIA STEIN, (totally) based on Mike Zingale's code, Spring/2013
    
        The linear advection equation provides a simple problem to 
        explore methods for hyperbolic problems
            at + uax = 0
        where u represents the speed at which the information propagates.

"""


import numpy
import pylab


class FDgrid:

    def __init__(self, nx, ng, xmin=0.0, xmax=1.0):

        self.xmin = xmin
        self.xmax = xmax
        self.ng = ng
        self.nx = nx
        
        # defome dim grid
        self.ilo = ng
        self.ihi = ng+nx-1

        # physical coords
        self.dx = (xmax - xmin)/(nx-1.0)
        self.x = xmin + (numpy.arange(nx+2*ng)-ng)*self.dx

        # storage for the solution
        self.a = numpy.zeros((nx+2*ng), dtype=numpy.float64)


    def scratchArray(self):
        """ return a scratch array dimensioned for our grid """
        return numpy.zeros((self.nx+2*self.ng), dtype=numpy.float64)
    

    def fillBCs(self):
        """ fill the a single ghostcell with periodic boundary conditions """
        self.a[self.ilo-1] = self.a[self.ihi-1]
        self.a[self.ihi+1] = self.a[self.ilo+1]




def tridiag(a, b, d):
    """ solve the linear system Ax = d where A has the form:

          a_i x_{i-1} + b_i x_i  = d_i """
    N = len(a)
    if not (len(b)  == len(d) == N):
        print "ERROR: vectors not the right size"
        return None

    # forward elimination
    dprime = numpy.zeros((N), dtype=a.dtype)
    dprime[0] = d[0]/b[0]

    for i in range(1,N-1):
        dprime[i] = (d[i] - dprime[i-1]*a[i])/(b[i])
        
    dprime[N-1] = (d[N-1] - dprime[N-2]*a[N-1])/(b[N-1])
    # back substitution
    x = numpy.zeros((N), dtype=a.dtype)
    x[N-1] = dprime[N-1]
    for i in reversed(range(0,N-1)):
        x[i] = dprime[i] 

    return x


"""
    define constants and i.c.
"""
# grid
nx = [65, 257]
ng = 1

# define the CFL and speed
C = [0.05, 0.1, 0.5, 1.0, 10.0]
u = 1.0 


"""
    loop in each value of C and nx (many different solutions
"""
for cc in range(len(C)):
    for nxx in range(len(nx)):
        
        g = FDgrid(nx[nxx], ng)
        g_im = FDgrid(nx[nxx], ng)
        
        # time info
        dt = C[cc]*g.dx/u
        t = 0.0
        tmax = 1.0*(g.xmax - g.xmin)/u      # On a domain [0,1], one period is simply: 1/u
        if ( (C[cc] >= 10.0) & (nx[nxx] == 65) ):
            tmax = 0.1*tmax
        if ( (C[cc] >= 10.0) & (nx[nxx] == 257) ):
            tmax = 0.0001*tmax
            
        # initialize the data -- tophat
        g.a[numpy.logical_and(g.x >= 0.333, g.x <= 0.666)] = 1.0
        g_im.a[numpy.logical_and(g_im.x >= 0.333, g_im.x <= 0.666)] = 1.0

        # create a copy of the initial tophat for plotting
        ainit = g.a.copy()

        # define arrays
        anew = g.scratchArray()
        anew_im = g_im.scratchArray()
        a = g_im.scratchArray()
        b = g_im.scratchArray()
        d = g_im.scratchArray()

        """
            loop in many periods (evolution loop)
        """
        
        """ I use your triagdia, so the solution is not as robust as if I had used
        some method such as numpy.linealg.solve(B,b) and the matrix is not triangular
        but the way the solver see it works"""
        i = g.ilo 
        while (i <= g.ihi+2):   # the two is because of the extra padding from the gosth cell, since im
            #compare with the explicit system, I decide to leave like thos
                a[i-1] =  -C[cc]
                b[i-1] = (C[cc] + 1)
                i += 1
            
        
        while (t < tmax):

            # fill the boundary conditions
            g.fillBCs()
            g_im.fillBCs()

            """
                 loop over zones: note since we are periodic and both endpoints
                 are on the computational domain boundary, we don't have to
                 update both g.ilo and g.ihi -- we could set them equal instead.
                 But this is more general
             """
            # upwind discretization explicit
            i = g.ilo    
            while (i <= g.ihi):
                anew[i] = g.a[i] - C[cc]*(g.a[i] - g.a[i-1])  
                i += 1

            
            # upwind discretization implicit               
            d[:] = g_im.a[:]
            anew_im= tridiag(a, b, d)
              
            
            # store the updated solution
            g.a[:] = anew[:]
            g_im.a[:] = anew_im[:]
        

            t += dt
            

        pylab.clf()
        pylab.cla()

        pylab.title("Upwind Adv. Eq.: C=%.2f, %d grid points" %(C[cc], nx[nxx]-1))
        
        pylab.plot(g.x[g.ilo:g.ihi+1], ainit[g.ilo:g.ihi+1], ls=":")
        pylab.plot(g.x[g.ilo:g.ihi+1], g.a[g.ilo:g.ihi+1], label="Explicit Solution")
        pylab.plot(g_im.x[g_im.ilo:g_im.ihi+1], g_im.a[g_im.ilo:g_im.ihi+1], label="Implicit Solution")
        
        pylab.xlabel('Grid', fontsize=10)
        pylab.ylabel('Solution', fontsize=10)
        
        if (C[cc] == 1):
            pylab.text(0.03, 0.95, "No diffusion for C=1 for Explicit!", color="r", fontsize=9,verticalalignment='top')
            pylab.text(0.03, 0.9, "Solution perfectly advected.", color="r", fontsize=9,verticalalignment='top')
        
        if (C[cc] <= 0.1):
            pylab.text(0.03, 0.95, "The most accurate solutions", color="r", fontsize=9,verticalalignment='top')
            pylab.text(0.03, 0.9, "still need a small C.", color="r", fontsize=9,verticalalignment='top')
        
        if (C[cc] > 1):
            pylab.text(0.03, 1.8, "Only implicit solutions", color="r", fontsize=9,verticalalignment='top')
            pylab.text(0.03, 1.3, "converge for any C!", color="r", fontsize=9,verticalalignment='top')
        
        leg = pylab.legend(loc=1,labelspacing=0.0005)
        ltext = leg.get_texts()
        pylab.setp(ltext, fontsize='small')
        leg.draw_frame(0)
        
        outfile = "implicit_advect_C-%0.2f_dx-%d.png" %(C[cc], nx[nxx]-1)
        pylab.savefig(outfile)


print("Done!")


    
