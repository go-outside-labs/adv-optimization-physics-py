# some examples on views and copies from the NumPy tutorial

import numpy as np


a = np.arange(12)


b = a   # no new object is created -- b and a point to the same object
print b is a   

b.shape = 3,4   # changes the shape of a too
print a.shape


print " "
print a
print b
print " "
b[1,1] = -1    # changes a[1,1] too
print a

print " "


# views / shallow copy
print "views"

c = a[:]   # or c = a.view()

print c is a
print c.base is a

c.shape = 12

print a
print c
print " "

c[2] = 100.0
print a
print c

print " "

# deep copy
print "deep copying"

d = a.copy()

print a
print d
print " "

d[:,:] = 0.0

print a
print d

