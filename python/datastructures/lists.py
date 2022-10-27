# more fun with lists


# create an empty list
myList = []

myList.append(1)
myList.append(2)

print myList

myList.append(20.0)
myList.append(-1.0)

myList.sort()

print myList
print " "

# get a specific element
print myList[1]
print " "

# remove (and return) the last item
end = myList.pop()
print end
print myList
print " "

# multiple instances of the same number
myList.append(1)
myList.append(1)

print myList
print myList.count(1)



