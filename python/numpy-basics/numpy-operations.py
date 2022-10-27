import numpy as np
import math

a = np.arange(12).reshape(3,4)
print a

# sum over the first axis (rows), so this gives the sum in each column
print a.sum(axis=0)

# sum over all elements
print a.sum()

# find the minimum and maximum element
print a.min(), a.max()

print " "


# universal functions
b = a*math.pi/12.0
print b

c = np.cos(b)
print c

print " "

# add two arrays
d = b + c
print d



