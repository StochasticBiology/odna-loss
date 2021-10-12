import re
import os.path
import sys

complexlabel = sys.argv[1]
interfacefilename = sys.argv[2]
pdbfilename = sys.argv[3]
extraspath = sys.argv[4]
mtoccurrencefilename = sys.argv[5]
ptoccurrencefilename = sys.argv[6]
outfilename = sys.argv[7]
plotscriptfilename = sys.argv[8]
plotfilename = sys.argv[9]

#complexlabel = "1oco"
#interfacefilename = "Prelims/1oco-lig.txt"
#pdbfilename = "Prelims/1oco-pdb.html"
#extraspath = "Prelims/"
#mtoccurrencefilename = "Data/mt-gene-occurrence.csv"
#ptoccurrencefilename = "Data/pt-gene-occurrence.csv"
#outfilename = "Data/1oco-data.csv"
#plotscriptfilename = "plot-1oco.py"
#plotfilename = "Data/plot-1oco.png"


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

fp = open(interfacefilename, "r")
# parse cut-and-pasted PISA output tables to get adjacency matrix
# awkwardly, these come in different formats for different complexes. we first examine the table header to figure out which format we're using, then count the number of fields in a given row to figure out which columns we should be extracting
formattype = 0
results = []
for line in fp.readlines():
    # detect the start of the results table, and set format type depending on what fields are present
    if "iNres" in line:
        if "Symmetry" in line:
            formattype = 1
        elif "Id" in line:
            formattype = 2
        else:
            formattype = 3
    # detect the end of the results table
    if "Copyright" in line or "Download" in line:
        formattype = 0
    row = line.split()
    rowdata = []
    # grab different fields depending on the recognised output format of the table
    if formattype == 1:
        if len(row) == 20:
            rowdata = [row[2], row[7], row[14]]
        if len(row) == 19:
            rowdata = [row[1], row[6], row[13]]
    if formattype == 2:
        if len(row) == 18:
            rowdata = [row[2], row[7], row[12]]
        if len(row) == 17:
            rowdata = [row[1], row[6], row[11]]
    if formattype == 3:
        rowdata = [row[1], row[6], row[11]]
    # if this line contained a record on an interface
    if len(rowdata) > 0:
        if is_number(rowdata[2]):
            print(str(rowdata)+"\n")
            # identify whether there is a ligand involved in this interface
            ligand = 0
            if "[" in rowdata[0] or "[" in rowdata[1]:
                ligand = 1
            # fold ligands out, replacing with interactions between the subunits they connect
            # the regex subs find [CCC]A:999, replaces with A, where CCC and 999 can be anything (i.e. replaces ligand bound to a chain with the chain ref)
            rowdata[0] = re.sub("\[[^\]]*\]", "", rowdata[0])
            rowdata[0] = re.sub(":.*", "", rowdata[0])
            rowdata[1] = re.sub("\[[^\]]*\]", "", rowdata[1])
            rowdata[1] = re.sub(":.*", "", rowdata[1])
            rowdata[2] = float(rowdata[2])
            # record energy and interface number (0 or 1) both counting and not counting ligand-mediated interfaces
            if rowdata[0] != rowdata[1]:
                if ligand == 1:
                    resultrow = [rowdata[0], rowdata[1], rowdata[2], 1, 0, 0]
                else:
                    resultrow = [rowdata[0], rowdata[1], rowdata[2], 1, rowdata[2], 1]
                results.append(resultrow)

fp.close()

# build dictionary of chains containing their energy and interface counts
dict = {}
for result in results:
    if result[0] in dict:
        dict[result[0]] = [sum(i) for i in zip(dict[result[0]], result[2:6])]
    else:
        dict[result[0]] = result[2:6]
    if result[1] in dict:
        dict[result[1]] = [sum(i) for i in zip(dict[result[1]], result[2:6])]
    else:
        dict[result[1]] = result[2:6]

# now read chain identities from the PDB html file
# this is finicky HTML parsing, subsetting out the elements of the messy HTML page that correspond to chain labels and gene labels
repdict = {}
fp = open(pdbfilename, "r")
uid = 1
for line in fp.readlines():
    if "Entity ID" in line:
        for newline in line.split("Entity ID"):
            cols = newline.split("<td")
            if len(cols) > 5:
                chains = cols[2].split(">")[2].split("<")[0].replace(" ", "").split(",")
                # specific fix here for one instance: the final entry in 5mlc
                if chains != ['']:
                    if "<strong>Gene Names</strong>:" in cols[5]:
                        name = cols[5].split("<strong>Gene Names</strong>:")[1].split("<br>")[0].split("&")[0].replace(" ", "").replace(",", "_")
                    else:
                        name = "uid-"+str(chains).replace(",","")
                        uid = uid+1
                    for chain in chains:
                        repdict[chain] = name

fp.close()

# specific fixes for some instances causing bugs
if "5mlc" in pdbfilename:
    repdict['A'] = "rna-A"
    repdict['B'] = "rna-B"
    repdict['C'] = "rna-C"
if "2h88" in pdbfilename:
    repdict['C'] = "sdhc"
    repdict['P'] = "sdhc"
if "2wsc" in pdbfilename:
    repdict['3'] = "lhca3"
if "5o31" in pdbfilename:
    repdict['P'] = "nad9"
    repdict['j'] = "nad2"
    repdict['k'] = "nad3"
    repdict['l'] = "nad8"
    
# compile a dictionary of interface details by gene label
resultsdict = {}
for key in dict:
    # this messy bit includes a stoichiometry count of 1 for each copy of a gene product found, so that these can be summed for overall stoichiometry
    # the count is added as a sublist, then the list of lists is flattened
    record = [[1], dict[key]]
    record = [item for sublist in record for item in sublist]
    if repdict[key] in resultsdict:
        resultsdict[repdict[key]] = [sum(i) for i in zip(resultsdict[repdict[key]], record)]
    else:
        resultsdict[repdict[key]] = record

# output results from this dictionary
fp = open(outfilename, "a")
for key in resultsdict:
    fp.write(complexlabel+","+key+","+str(resultsdict[key]).replace("[","").replace("]","").replace(" ","")+","+str(resultsdict[key][1]/resultsdict[key][0])+"\n")

fp.close()

# label aliases set contains list of organelle-encoded subunits
fp = open("Prelims/label-aliases.csv", "r")
aliasdict = {}
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    aliasdict[lineset[0]] = lineset[1]

fp.close()

occurdict = {}
fp = open(mtoccurrencefilename, "r")
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if lineset[0] != "GeneLabel":
        occurdict[lineset[0]] = float(lineset[1])

fp.close()
fp = open(ptoccurrencefilename, "r")
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if lineset[0] != "GeneLabel":
        occurdict[lineset[0]] = float(lineset[1])

fp.close()

# construct Python script for pymol plotting
fp = open(plotscriptfilename, "w")
#  -- preamble 
fp.write('from pymol import cmd\n')
#fp.write('cmd.load("'+complexlabel+'")\n')
fp.write('cmd.fetch("'+complexlabel+'")\n')
fp.write('cmd.bg_color("white")\n')
fp.write('cmd.set("transparency_mode", 1)\n')
fp.write('cmd.select("sele", "all")\n')
fp.write('cmd.hide("lines", "sele")\n')
fp.write('cmd.set("surface_color", "white", "sele")\n')
fp.write('cmd.set("transparency", 0.75, "sele")\n')

if complexlabel not in ("5o31", "2h88", "5xte", "1oco", "6cp3"):
    col1 = "0x0000FF"
    col2 = "0x4444FF"
    col3 = "0x8888FF"
    col4 = "0xBBBBFF"
else:
    col1 = "0xFF0000"
    col2 = "0xFF4444"
    col3 = "0xFF8888"
    col4 = "0xFFBBBB"

for chain in repdict:
    gene = repdict[chain]
    if gene in aliasdict:
        gene = aliasdict[gene]
    gene = gene.lower()
    if gene in occurdict:
        fp.write('cmd.select("sele", "chain '+chain+'")\n')
        if occurdict[gene] > 1000:
            fp.write('cmd.set("surface_color", "'+col1+'", "sele")\n')
        elif occurdict[gene] > 100:
            fp.write('cmd.set("surface_color", "'+col2+'", "sele")\n')
        elif occurdict[gene] > 10:
            fp.write('cmd.set("surface_color", "'+col3+'", "sele")\n')
        else:
            fp.write('cmd.set("surface_color", "'+col4+'", "sele")\n')
        fp.write('cmd.set("transparency", 0.5, "sele")\n')

fp.write('cmd.select("sele", "all")\n')
fp.write('cmd.show("surface", "sele")\n')

extrafilename = extraspath+"plot-"+complexlabel+"-extras.txt"
if os.path.isfile(extrafilename):
    extrafp = open(extrafilename)
    for line in extrafp.readlines():
        fp.write(line)
    extrafp.close()


# if present, append extras to pymol script (view, monomer, etc)
#more $1-extras.txt >> $1-vis.py

# outro to pymol script, showing and saving image
fp.write('cmd.hide("nonbonded")\n')
fp.write('cmd.select("none")\n')
# keep these for transparent backgrounds... makes the rendering much more computationally intensive though!
#fp.write('cmd.set("ray_opaque_background", 0)\n')
#fp.write('cmd.ray(600,600)\n')
fp.write('cmd.png("'+plotfilename+'")\n')
fp.write('cmd.quit()\n')
fp.close()



