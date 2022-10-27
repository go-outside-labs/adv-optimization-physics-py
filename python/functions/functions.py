# functions

globalVariable = 5


def printFunc(numTimes, msg="my function"):
    
    i = 0
    while (i < numTimes and i < globalVariable):
        print msg
        i += 1

    a = ["this", "is", "my", "list"]

    return a, msg



list, printedMsg = printFunc(2, "message")

print list
print printedMsg

# note here that the global variable comes into play
printFunc(10)






