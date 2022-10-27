# example of error estimation with numerical derivatives and using
# Richardson extrapolation to reduce the leading order error.

import math


# function we are differentiating
def fun(x):
    return math.sin(x)


# analytic derivative (for comparison)
def fprime(x):
    return math.cos(x)


# difference equation
def diff(x, h, fun):
    return (fun(x+h) - fun(x-h))/(2*h)


# desired tolerance -- be careful not to go too close to machine
# epsilon, or else roundoff error will rule
tol = 1.e-7


# starting h for differencing
h = 0.125

# point where we want the derivative
x0 = 1.0

err = 100.0

# initial derivative
d0 = diff(x0, h, fun)

print "h, d, rel err, analytic rel err"

while (err > tol):
    d1 = diff(x0, h/2, fun)
    
    # relative error between the h and h/2 estimates
    err = abs(d1 - d0)/d1

    # combination of h and h/2 estimates to eliminate leading error
    # term
    d = (4*d1-d0)/3.0

    print h, d, err, abs(d - fprime(x0))/fprime(x0)

    d0 = d1
    h = h/2



    
