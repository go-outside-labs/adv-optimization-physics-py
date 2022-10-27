#---------------------------------------------------------------------------
# integers
a = 1  
b = 2

print a/b
print a+b


#---------------------------------------------------------------------------
# floats
r = 2.0  
pi = 3.14159

A = pi*r**2

print A
print a  # note case sensitive


#---------------------------------------------------------------------------
# strings
s = "str"  
t = "ing"

print s, t

u = s + t
print u

# slicing
print u[0:2]
print u[-1]


#---------------------------------------------------------------------------
# lists
a = [1, 2.0, "three", 4]
print a
print a[0]
print len(a)

# copy -- this is just a pointer
b = a
print a
print b
a[0] = "one"
print b

# shallow copy
c = a[:]
print a
print c
a[1] = "two"
print a
print c

