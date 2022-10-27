import numpy as np

# example from scipy.org NumPy tutorial
def f(x,y):
    return 10*x+y

a = np.fromfunction(f, (5,4),dtype=int)

print a

print " "
print a[0:2,0:2]

print " "
print a[:,1]

print " "
print a.flatten()


print " "
for row in a:   # iteration is done over fiest axis
    print row


print " "
for element in a.flat:
    print element,





