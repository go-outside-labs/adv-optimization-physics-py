import numpy as np

# basic array creation
a = np.arange(15).reshape(3,5)
print a

print a.ndim
print a.shape
print a.size
print a.dtype
print a.itemsize

print " "

# create an array from a list
b = np.array( [1.0,2.0,3.0,4.0])
print b
print b.dtype

print " "

# create an empty array of specified size
c = np.zeros((10,7), dtype=np.float64)
print c

# create an array with 10 elements ranging from 0 to 1
d = np.linspace(0, 1, 10)
print d


