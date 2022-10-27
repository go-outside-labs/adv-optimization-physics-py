# create a simple dictionary
myDict = {"key1":1, "key2":2, "key3":3}

print myDict
print myDict["key1"]
print " "


# add a key:pair -- notice that the values can be any data type
myDict["newkey"] = "new"
print myDict

# and so can the keys
myDict[5] = 5
print myDict
print " "


# loop over the elements
for k, v in myDict.iteritems():
    print "key = %s, value = %s" % (`k`, `v`)

print " "

# just get the keys
keys = myDict.keys()
print keys

print " "

# check whether a key exists
print "key1" in keys

print "dummykey" not in keys



