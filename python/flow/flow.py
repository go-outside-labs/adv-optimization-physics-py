# control flow -- here white space is important

# while loop
n = 0
while n < 10:
    print n
    n += 1



# if statements
x = 10
if x < 0:
    print "negative"
elif x == 0:
    print "zero"
else:
    print "positive"



# note: there is no case statement like in C/Fortran


# looping over items in a list
alist = [1, 2.0, "three", 4]
for a in alist:
    print a



# break
alist = [1, 2.0, "three", 4]
n = 0
for a in alist:
    if a == "three":
        break
    else:
        n += 1

print n

# simpler way
print alist.index("three")








