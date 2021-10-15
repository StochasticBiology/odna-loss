import sys

# open files for IO
residuestatsfile = sys.argv[1]
codonstatsfile = sys.argv[2] 
extension = sys.argv[3]

outstr = extension

f = open(residuestatsfile, "r")
lineset = f.readline().rstrip("\n").split(",")
outstr = outstr+str(lineset[1:]).replace(", ", ",").replace("[","").replace("]","").replace("'","") + ","
f.close()

f = open(codonstatsfile, "r")
lineset = f.readline().rstrip("\n").split(",")
outstr = outstr+str(lineset[1:]).replace(", ", ",").replace("[","").replace("]","").replace("'","") 
f.close()
print(outstr)
