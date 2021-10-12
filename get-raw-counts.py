import sys

srcfilename = sys.argv[1]
rawoccurrencefilename = sys.argv[2]

# open and read source data file
infile = open(srcfilename, "r")
lines = infile.readlines()

speciesrawgenedict = {}
for line in lines:
    if line[0] == '>':
        lineset = line.rstrip("\n").split(",")
        genelabel = lineset[1]
        if genelabel in speciesrawgenedict:
            speciesrawgenedict[genelabel].append(lineset[0])
        else:
            speciesrawgenedict[genelabel] = [lineset[0]]

rawgenedict = {}
for gene in speciesrawgenedict:
    rawgenedict[gene] = len(set(speciesrawgenedict[gene]))

# write occurrence counts to file
fp = open(rawoccurrencefilename, "w")
fp.write("GeneLabel,Occurrence\n")
for gene in rawgenedict:
    fp.write(gene+","+str(rawgenedict[gene])+"\n")

fp.close()
