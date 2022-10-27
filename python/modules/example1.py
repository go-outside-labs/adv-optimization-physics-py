import fibo

fibo.fib(1000)

fibList = fibo.fib2(100)
print " "

print fibList

# can assign the function name to a variable to make it easier to call
f = fibo.fib2
print f(100)


# find out the names that a module defines
print dir(fibo)
