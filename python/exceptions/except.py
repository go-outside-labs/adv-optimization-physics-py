# try openning a file for reading
try: f = open("file.txt", "r")
except IOError:
    print "I/O Error"


# undefined variable
x = 1.0
try: x + y
except NameError:
    print "undefined variable"



# example from tutorial
def this_fails():
    x = 1/0

try:
    this_fails()
except ZeroDivisionError as detail:
    print 'Handling run-time error:', detail

