"""
    We construct the vector field 
        Ustar = Ud + G\phi
    where G is the discrete gradient.
    
    MIA STEIN, spring/2013
    
"""


import numpy
import math
import pylab
import operator



def scalarField(x,y):
    
    """ this is the analytical function, it can be padded, doenst matter,
    there results are in function of spacing"""
    phi = 0.1*numpy.cos(2.0*math.pi*y)*numpy.cos(2.0*math.pi*x)
    return phi


def calculateDivergenceFreeTerm(x,y):
    
    """ this is the analytical function, it can be padded, doenst matter,
    there results are in function of spacing"""
    Udx = - numpy.sin(math.pi*x)**2*numpy.sin(2.0*math.pi*y) 
    Udy = numpy.sin(math.pi*y)**2*numpy.sin(2.0*math.pi*x)
    return [Udx, Udy]


def do_Plot_item_a(Ustar, U, grad_phi):
    """ plotting stuff: here we plot the first row with Ustar and
        the second rows with Ud and grad phi, to see how they sum up """
    pylab.clf()
    pylab.cla()
    
    f = pylab.figure() 
    f.text(.5, .95, r"Top: $U = U_d + \nabla \phi$, where "  r"$U_d$ =  - $\sin^2 ( \pi x)\sin(2 \pi y)\hat x$ + $\sin^2(\pi y)\sin(2\pi x)\hat y$ and $\phi = \frac{1}{10} \cos(2\pi y) \cos(2\pi x$) ",  horizontalalignment='center', size=9)

    pylab.subplot(221) 
    pylab.imshow(Ustar[0])
    pylab.xlabel("Vector "r"$U_d$ in $\hat x$:" ,size=8)
    pylab.ylabel("U in terms of # of cells in " r"$\hat x$", size=8)

    pylab.subplot(222)
    pylab.imshow(Ustar[1])
    pylab.xlabel(r"Vector $\nabla \phi$ in $ \hat x $:", size=8)
    pylab.ylabel("U in terms of # of cells in " r"$\hat y$",size=8)

    pylab.subplot(223)
    pylab.imshow(U [0])
    pylab.ylabel("# of cells",size=8)
    pylab.xlabel("# of cells",size=8)
    
    
    pylab.subplot(224)
    pylab.imshow(grad_phi [0])
    pylab.xlabel("# of cells",size=8)
    pylab.ylabel("# of cells",size=8)
    pylab.savefig("plots/item_a_Ustar.png")
    
    return 0


def doPartA(nx, ny, ng, dx, dy, DO_PLOTS):
    """ we test if does any difference to
        start place the first point before 0 and the last after 1.
        So here is the part  didnt understand, it does change your
        solution and the speed of converegence, and I havent figure it
        out yet..."""
    #x = numpy.linspace(-dx, 1.0 + dx, nx + 2.0*ng)
    #y = numpy.linspace(-dx, 1.0 + dy, ny + 2.0*ng)
    
    x = (numpy.arange(nx + 2.0*ng) -ng + 0.5) * dx
    y = (numpy.arange(ny + 2.0*ng) -ng + 0.5) * dy

    print x
    print y

    """ the bug in my code is here, but I cant figure it yet, I know
    it will disturb the gradient etc, but what the correct way? both
    ways are wrong and give different answers!!!"""
    
    

    """ We generate the 2d from here, that's the shape of the final vector
        array"""
    x2d = numpy.repeat(x, ny + 2.0*ng)
    x2d.shape = (nx+2.0*ng, ny + 2.0*ng)

    y2d = numpy.repeat(y, nx + 2.0*ng)
    y2d.shape = (ny+2.0*ng, nx + 2.0*ng)
    y2d = numpy.transpose(y2d)
    
    print x2d[:,0]
    print y2d[0,:]
   
    """ calculate compoents """
    phi = scalarField(x2d ,y2d)
    U = calculateDivergenceFreeTerm(x2d,y2d)

    """ 
        The gradient is computed using central differences in the interior and 
        first differences at the boundaries. 
        The returned gradient hence has the same shape as the input array.
    """
    grad_phi = numpy.gradient(phi, dx, dy)
    Ustar = map(operator.add, grad_phi, U)
    
    
    """ unpadding """
    Ustar_finalx  = Ustar[0]
    Ustar_finaly = Ustar[1]
    aux_Ux = Ustar_finalx[ng:ng + nx ]
    aux_Uy = Ustar_finaly[ng:ng + nx ]
    Ustar_final = [aux_Ux, aux_Uy]

    if (DO_PLOTS == 1):
        do_Plot_item_a(Ustar, U, grad_phi)
                         
    return Ustar, Ustar_final, phi, U
