import sys

threshold = float(sys.argv[1])
srcfilename = sys.argv[2]
repfilename = sys.argv[3]
statsfilename = sys.argv[4]
barcodefilename = sys.argv[5]
speciesfilename = sys.argv[6]
occurrencefilename = sys.argv[7]

#threshold = 10
#srcfilename = "Data/mt-dna.fasta"
#repfilename = "Data/mt-blast-replace.csv"
#statsfilename = "tmp1"
#barcodefilename = "tmp2"
#speciesfilename = "tmp3"
#occurrencefilename = "tmp4"

if "mt" in srcfilename:
    compartment = "MT" 
elif "pt" in srcfilename:
    compartment = "PT" 
else:
    compartment = "??" 

# read in replacement dictionary
infile = open(repfilename, "r")
lines = infile.readlines()
replacements = {}
for line in lines:
    lineset = line.rstrip("\n").split(",")
    replacements[lineset[0]] = lineset[1]

infile.close()
# print(replacements)

# open and read source data file
infile = open(srcfilename, "r")
lines = infile.readlines()

# first, populate a dictionary of "true" gene labels (from replacement dictionary) with their occurrence frequency
# this will help prune rare genes later
speciesgenedict = {}
for line in lines:
    if line[0] == '>':
        lineset = line.rstrip("\n").split(",")
        genelabel = lineset[1]
        if genelabel in replacements:
            truelabel = replacements[genelabel]
        else:
            truelabel = genelabel
        if truelabel in speciesgenedict:
            speciesgenedict[truelabel].append(lineset[0])
        else:
            speciesgenedict[truelabel] = [lineset[0]]

genedict = {}
for gene in speciesgenedict:
    genedict[gene] = len(set(speciesgenedict[gene]))

# write occurrence counts to file
fp = open(occurrencefilename, "w")
fp.write("GeneLabel,Occurrence\n")
for gene in genedict:
    fp.write(gene+","+str(genedict[gene])+"\n")

fp.close()


# construct and sort the set of genes of interest
geneset = [gene for gene in genedict if genedict[gene] > threshold and "orf" not in gene and "ymf" not in gene and "ycf" not in gene and "_" not in gene and "anon" not in gene and "-i" not in gene and "oi" not in gene]
geneset.sort()
print("From file "+srcfilename+" we have "+str(len(geneset))+" genes of interest:\n"+str(geneset)+"\n")

# now, read through the source data again, record genes of interest that are present for each species
# and outputting their summary statistics to file
outfile = open(statsfilename, "a")
#outfile = open(statsfilename, "w")
#outfile.write("Species,Compartment,GeneLabel,Length,Hydro,Hydro_i,MolWeight,pKa1,pKa2,A_Glu,CW,GC,Uni1,Uni2,Robust,GC12,GC3\n")
speciesdict = {}
for line in lines:
    if line[0] == '>':
        lineset = line.rstrip("\n").split(",")
        species = lineset[0][1:]
        genelabel = lineset[1]
        if genelabel in replacements:
            truelabel = replacements[genelabel]
        else:
            truelabel = genelabel
        if truelabel in geneset:
            if species in speciesdict:
                speciesdict[species].append(truelabel)
            else:
                speciesdict[species] = [truelabel]
            if lineset[4] != "statserror":
                outfile.write(species.replace("_", " ")+","+compartment+","+truelabel+","+str(lineset[4:18]).replace(", ",",").replace("[","").replace("]","").replace("'","")+"\n")

outfile.close()

# construct summary barcodes of presence/absence for genes of interest
outfile = open(barcodefilename, "w")
specieslistfile = open(speciesfilename, "w")
outfile.write("Species,"+str(geneset).replace(", ", ",").replace("[", "").replace("]", "").replace("'","")+"\n")
barcodes = {}
for species in speciesdict:
    barcodes[species] = []
    for gene in geneset:
        if gene in speciesdict[species]:
            barcodes[species].append(1)
        else:
            barcodes[species].append(0)
    outfile.write(species.replace("_", " ")+","+str(barcodes[species]).replace(", ", ",").replace("[", "").replace("]", "")+"\n")
    specieslistfile.write(species.replace("_", " ")+"\n")

outfile.close()
specieslistfile.close()


