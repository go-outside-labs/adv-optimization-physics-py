"""
The patch module allows for a grid to be created and for data to be
defined on that grid.

Typical usage:

create the grid

   grid = grid2d(nx, ny)


create the data that lives on that grid

   data = ccData2d(grid)

   bcObj = bcObject(xlb="reflect", xrb="reflect", 
                    ylb="outflow", yrb="outflow")
   data.registerVar("density", bcObj)
   ...

   data.create()


initialize some data

   dens = data.getVarPtr("density")
   dens[:,:] = ...


fill the ghost cells

   data.fillBC("density")

M. Zingale (2013-03-28)

"""

import numpy
import sys



class bcObject:
    """
    boundary condition container -- hold the BCs on each boundary
    for a single variable
    """

    def __init__ (self, 
                  xlb="outflow", xrb="outflow", 
                  oddReflectDir=""):

        valid = ["outflow", "periodic", 
                 "reflect", "reflect-even", "reflect-odd",
                 "dirichlet", "neumann"]

        # note: "reflect" is ambiguous and will be converted into
        # either reflect-even (the default) or reflect-odd if
        # oddReflectDir specifies the corresponding direction ("x",
        # "y")

        # -x boundary
        if (xlb in valid):
            self.xlb = xlb
            if (self.xlb == "reflect"):
                if (oddReflectDir == "x"):
                    self.xlb = "reflect-odd"
                else:
                    self.xlb = "reflect-even"
            
        else:
            sys.exit("ERROR: xlb = %s invalid BC" % (xlb))

        # +x boundary
        if (xrb in valid):
            self.xrb = xrb
            if (self.xrb == "reflect"):
                if (oddReflectDir == "x"):
                    self.xrb = "reflect-odd"
                else:
                    self.xrb = "reflect-even"

        else:
            sys.exit("ERROR: xrb = %s invalid BC" % (xrb))


        # periodic checks
        if ((xlb == "periodic" and not xrb == "periodic") or
            (xrb == "periodic" and not xlb == "periodic")):
            sys.exit("ERROR: both xlb and xrb must be periodic")

    def __str__(self):
        """ print out some basic information about the BC object """

        string = "BCs: -x: %s  +x: %s " % \
            (self.xlb, self.xrb)

        return string

    
class grid1d:
    """
    the 1-d grid class.  The grid object will contain the coordinate
    information (at various centerings).

    A basic (1-d) representation of the layout is:

    |     |      |     X     |     |      |     |     X     |      |     |
    +--*--+- // -+--*--X--*--+--*--+- // -+--*--+--*--X--*--+- // -+--*--+
       0          ng-1    ng   ng+1         ... ng+nx-1 ng+nx      2ng+nx-1

                         ilo                      ihi
    
    |<- ng ghostcells->|<---- nx interior zones ----->|<- ng ghostcells->|

    The '*' marks the data locations.
    """

    def __init__ (self, nx, ng=1, xmin=0.0, xmax=1.0):

        """
        The class constructor function.

        The only data that we require is the number of points that
        make up the mesh.

        We optionally take the extrema of the domain, number of ghost
        cells (assume 1)
        """

        # size of grid
        self.nx = nx
        self.ng = ng

        self.qx = 2*ng+nx

        # domain extrema
        self.xmin = xmin
        self.xmax = xmax

        # compute the indices of the block interior (excluding guardcells)
        self.ilo = ng
        self.ihi = ng+nx-1

        # define the coordinate information at the left, center, and right
        # zone coordinates
        self.dx = (xmax - xmin)/nx

        self.xl = (numpy.arange(nx+2*ng) - ng)*self.dx + xmin
        self.xr = (numpy.arange(nx+2*ng) + 1.0 - ng)*self.dx + xmin
        self.x = 0.5*(self.xl + self.xr)

    def scratchArray(self, dtype=numpy.float64):
        """ 
        return a standard numpy array dimensioned to have the size
        and number of ghostcells as the parent grid
        """

        return numpy.zeros((2*self.ng+self.nx), dtype=dtype)

    def __str__(self):
        """ print out some basic information about the grid object """

        return "1-d grid: nx = " + `self.nx` + ", ng = " + `self.ng`

    def __eq__(self, other):
        """ are two grids equivalent? """
        result = (self.nx == other.nx) and (self.ng == other.ng) and \
                 (self.xmin == other.xmin) and (self.xmax == other.xmax)

        return result


class ccData1d:
    """
    the cell-centered data that lives on a grid.

    a ccData1d object is built in a multi-step process before it can
    be used.  We pass in a grid object to describe where the data
    lives:

        myData = patch.ccData1d(myGrid)

    register any variables that we expect to live on this patch.  Here
    bcObject describes the boundary conditions for that variable.

        myData.registerVar('density', bcObject)
        myData.registerVar('x-momentum', bcObject)
        ...

    finally, finish the initialization of the patch

        myPatch.create()

    This last step actually allocates the storage for the state
    variables.  Once this is done, the patch is considered to be
    locked.  New variables cannot be added.
    
    """

    def __init__ (self, grid, dtype=numpy.float64):
        
        self.grid = grid

        self.dtype = dtype
        self.data = None

        self.vars = []
        self.nvar = 0

        self.BCs = {}

        self.initialized = 0


    def registerVar(self, name, bcObject):
        """ 
        register a variable with ccData1d object.  Here we pass in a
        bcObject that describes the boundary conditions for that
        variable.
        """

        if (self.initialized == 1):
            sys.exit("ERROR: grid already initialized")

        self.vars.append(name)
        self.nvar += 1

        self.BCs[name] = bcObject


    def create(self):
        """
        called after all the variables are registered and allocates
        the storage for the state data
        """

        if (self.initialized) == 1:
            sys.exit("ERROR: grid already initialized")

        self.data = numpy.zeros((self.nvar,
                                 2*self.grid.ng+self.grid.nx), 
                                dtype=self.dtype)
        self.initialized = 1

        
    def __str__(self):
        """ print out some basic information about the ccData2d object """

        if (self.initialized == 0):
            myStr = "ccData1d object not yet initialized"
            return myStr

        myStr = "cc data: nx = " + `self.grid.nx` + \
                       ", ng = " + `self.grid.ng` + "\n" + \
                 "   nvars = " + `self.nvar` + "\n" + \
                 "   variables: \n" 
                 
        ilo = self.grid.ilo
        ihi = self.grid.ihi

        n = 0
        while (n < self.nvar):
            myStr += "%16s: min: %15.10f    max: %15.10f\n" % \
                (self.vars[n],
                 numpy.min(self.data[n,ilo:ihi+1]), 
                 numpy.max(self.data[n,ilo:ihi+1]) )
            myStr += "%16s  BCs: -x: %-12s +x: %-12s \n" %\
                (" " , self.BCs[self.vars[n]].xlb, 
                       self.BCs[self.vars[n]].xrb)
            n += 1
 
        return myStr
    
    def getVarPtr(self, name):
        """
        return a pointer to the data array for the variable described
        by name.  Any changes made to this are automatically reflected
        in the ccData2d object
        """
        n = self.vars.index(name)
        return self.data[n,:]

    def zero(self, name):
        n = self.vars.index(name)
        self.data[n,:] = 0.0
        
    def fillBCAll(self):
        """
        fill boundary conditions on all variables
        """
        for name in self.vars:
            self.fillBC(name)

                
    def fillBC(self, name):
        """ 
        fill the boundary conditions.  This operates on a single state
        variable at a time, to allow for maximum flexibility

        we do periodic, reflect-even, reflect-odd, and outflow

        each variable name has a corresponding bcObject stored in the
        ccData2d object -- we refer to this to figure out the action
        to take at each boundary.
        """
    
        # there is only a single grid, so every boundary is on
        # a physical boundary (except if we are periodic)

        # Note: we piggy-back on outflow and reflect-odd for
        # Neumann and Dirichlet homogeneous BCs respectively, but
        # this only works for a single ghost cell

    
        n = self.vars.index(name)

        # -x boundary
        if (self.BCs[name].xlb == "outflow" or 
            self.BCs[name].xlb == "neumann"):

            i = 0
            while i < self.grid.ilo:
                self.data[n,i] = self.data[n,self.grid.ilo]
                i += 1                

        elif (self.BCs[name].xlb == "reflect-even"):
        
            i = 0
            while i < self.grid.ilo:
                self.data[n,i] = self.data[n,2*self.grid.ng-i-1]
                i += 1

        elif (self.BCs[name].xlb == "reflect-odd" or
              self.BCs[name].xlb == "dirichlet"):
        
            i = 0
            while i < self.grid.ilo:
                self.data[n,i] = -self.data[n,2*self.grid.ng-i-1]
                i += 1

        elif (self.BCs[name].xlb == "periodic"):

            i = 0
            while i < self.grid.ilo:
                self.data[n,i] = self.data[n,self.grid.ihi-self.grid.ng+i+1]
                i += 1
            

        # +x boundary
        if (self.BCs[name].xrb == "outflow" or
            self.BCs[name].xrb == "neumann"):

            i = self.grid.ihi+1
            while i < self.grid.nx+2*self.grid.ng:
                self.data[n,i] = self.data[n,self.grid.ihi]
                i += 1
                
        elif (self.BCs[name].xrb == "reflect-even"):

            i = 0
            while i < self.grid.ng:
                i_bnd = self.grid.ihi+1+i
                i_src = self.grid.ihi-i

                self.data[n,i_bnd] = self.data[n,i_src]
                i += 1

        elif (self.BCs[name].xrb == "reflect-odd" or
              self.BCs[name].xrb == "dirichlet"):

            i = 0
            while i < self.grid.ng:
                i_bnd = self.grid.ihi+1+i
                i_src = self.grid.ihi-i
                
                self.data[n,i_bnd] = -self.data[n,i_src]
                i += 1

        elif (self.BCs[name].xrb == "periodic"):

            i = self.grid.ihi+1
            while i < 2*self.grid.ng + self.grid.nx:
                self.data[n,i] = self.data[n,i-self.grid.ihi-1+self.grid.ng]
                i += 1

    def restrict(self, varname):
        """
        restrict the variable varname to a coarser grid (factor of 2
        coarser) and return an array with the resulting data (and same
        number of ghostcells)            
        """

        fG = self.grid
        fData = self.getVarPtr(varname)

        # allocate an array for the coarsely gridded data
        ng_c = fG.ng
        nx_c = fG.nx/2

        cData = numpy.zeros((2*ng_c+nx_c), dtype=self.dtype)

        ilo_c = ng_c
        ihi_c = ng_c+nx_c-1

        # fill the coarse array with the restricted data -- just
        # average the 2 fine cells into the corresponding coarse cell
        # that encompasses them.

        # This is done by shifting our view into the fData array and
        # using a stride of 2 in the indexing.
        cData[ilo_c:ihi_c+1] = \
            0.5*(fData[fG.ilo  :fG.ihi+1:2] + fData[fG.ilo+1:fG.ihi+1:2])
                  
        return cData


    def prolong(self, varname):
        """
        prolong the data in the current (coarse) grid to a finer
        (factor of 2 finer) grid.  Return an array with the resulting
        data (and same number of ghostcells).

        We will reconstruct the data in the zone from the
        zone-averaged variables using the centered-difference slopes

                  (x)      
        f(x,y) = m    x/dx + <f> 

        When averaged over the parent cell, this reproduces <f>.

        Each zone's reconstrution will be averaged over 2 children.  

        |           |     |     |     |
        |    <f>    | --> |     |     | 
        |           |     |  1  |  2  |  
        +-----------+     +-----+-----+ 

        We will fill each of the finer resolution zones by filling all
        the 1's together, using a stride 2 into the fine array.  Then
        the 2's, this allows us to operate in a vector
        fashion.  All operations will use the same slopes for their
        respective parents.

        """

        cG = self.grid
        cData = self.getVarPtr(varname)

        # allocate an array for the coarsely gridded data
        ng_f = cG.ng
        nx_f = cG.nx*2

        fData = numpy.zeros((2*ng_f+nx_f), dtype=self.dtype)

        ilo_f = ng_f
        ihi_f = ng_f+nx_f-1

        # slopes for the coarse data
        m_x = cG.scratchArray()
        m_x[cG.ilo:cG.ihi+1] = \
            0.5*(cData[cG.ilo+1:cG.ihi+2] - cData[cG.ilo-1:cG.ihi])


        # fill the '1' children
        fData[ilo_f:ihi_f+1:2] = \
            cData[cG.ilo:cG.ihi+1] - 0.25*m_x[cG.ilo:cG.ihi+1] 

        # fill the '2' children
        fData[ilo_f+1:ihi_f+1:2] = \
            cData[cG.ilo:cG.ihi+1] + 0.25*m_x[cG.ilo:cG.ihi+1] 
                  
        return fData
        

if __name__== "__main__":

    # illustrate basic mesh operations

    myg = grid1d(16, xmax=1.0)

    mydata = ccData1d(myg)

    bc = bcObject()

    mydata.registerVar("a", bc)
    mydata.create()


    a = mydata.getVarPtr("a")
    a[:] = numpy.exp(-(myg.x - 0.5)**2/0.1**2)

    print mydata
