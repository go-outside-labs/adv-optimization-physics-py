"""

    This module calculates the machine precision
    MIA STEIN (based on Mike Zingale's code), spring 2013.

"""

def calculatePrecision():

    x = 1.0
    epsilon = 1.0

    while (not x + epsilon == x):
        epsilon = epsilon/2.0

    return 2*epsilon