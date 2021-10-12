from ete3 import Tree #, TreeStyle, NodeStyle, TextFace
import sys

barcodefilename = sys.argv[1]
treefilename = sys.argv[2]
transfilename = sys.argv[3]
sistersfilename = sys.argv[4]

#barcodefilename = "Data/mt-barcodes-manual.csv"
#treefilename = "Prelims/mt-tree-manual.phy"
#transfilename = "Data/mt-trans-manual.txt"
 
# first read barcodes for all species
fp = open(barcodefilename, "r")
species = []
barcodes = {}
genecount = {}
for line in fp.readlines():
    linesplit = line.rstrip("\n").split(",")
    species.append(linesplit[0])
    barcodes[linesplit[0]] = linesplit[1:]
    if linesplit[0] != "Species":
        genecount[linesplit[0]] = sum([float(element) for element in barcodes[linesplit[0]]])
    ntraits = len(barcodes[linesplit[0]])

fp.close()
genenames = barcodes['Species']

# now read Newick tree. because some NCBI entries awkwardly include parentheses and colons, we'll strip these out line-by-line, remembering if a line started with a parenthesis, before joining everything together
fp = open(treefilename, "r")
treeraw = []
for line in fp.readlines():
    if line[0] == "(":
        treeraw.append("(")
    if line[0] == ")":
        treeraw.append(")")
    treeraw.append(line.replace("(", "").replace(")", "").replace(": ", ""))

treeraw = "".join(treeraw)
tree = Tree(treeraw.lower().replace("'",""), format=1)

# go through our tree from leaves up
change = 1
print("Pruning\n")
while change == 1:
    change = 0
    toprune = []
    for node in tree.traverse("postorder"):
        if node.is_leaf() and node.name not in species:
            toprune.append(node.name)
    print(str(toprune)+"\n")
    for nodename in toprune:
        node = tree.search_nodes(name=nodename)[0]
        node.delete()
        change = 1

# go through our tree from leaves up
for node in tree.traverse("postorder"):
    # A verbose output
    if not node.is_leaf() and node.name not in barcodes: # that is, for all internal (ancestral) nodes
        # reference (end state) vector
        ref = ['0']*ntraits
        # go through children and assign 1s to ancestor if a descendant has them
        for childnode in node.children:
            for i, c in enumerate(barcodes[childnode.name]):
                if c == '1':
                    ref[i] = '1'
        genecount[node.name] = sum([float(element) for element in ref])
        barcodes[node.name] = ref

# output reconstructed nodes 
fp = open(transfilename, "w")
fp.write("Ancestor,Descendant,")
fp.write(str([f'a.{i}' for i in barcodes['Species']]).replace("[", "").replace("]", "").replace("'", "").replace(" ",""))
fp.write(",")
fp.write(str([f'd.{i}' for i in barcodes['Species']]).replace("[", "").replace("]", "").replace("'", "").replace(" ",""))
fp.write("\n")
for node in tree.traverse("postorder"):
    # A verbose output
    if not node.is_leaf(): # For all internal nodes
        for childnode in node.children:
            print("".join(barcodes[node.name]), "(", node.name, ") -> ", "".join(barcodes[childnode.name]), "(", childnode.name, ")")
            if barcodes[node.name] != barcodes[childnode.name]:
                fp.write(node.name+","+childnode.name+",")
                fp.write(str(barcodes[node.name]).replace("[","").replace("]","").replace("'","").replace(" ",""))
                fp.write(",")
                fp.write(str(barcodes[childnode.name]).replace("[","").replace("]","").replace("'","").replace(" ",""))
                fp.write("\n")
#                print(" ".join(str(x) for x in list(barcodes[node.name])), file=fp)
#                print(" ".join(str(x) for x in list(barcodes[childnode.name])), file=fp)

fp.close()


# first pass at identifying sets of sister nodes with different gene counts
fp = open(sistersfilename, "w")
fp.write("ParentNode,OffspringDifferences\n")
for node in tree.traverse("postorder"):
    if not node.is_leaf(): # For all internal nodes
        thisparentdict = {}
        allleaves = 1
        for childnode in node.children:
            if not childnode.is_leaf():
                allleaves = 0
            noderef = int(genecount[childnode.name])
            if noderef not in thisparentdict:
                thisparentdict[noderef] = [childnode.name]
            else:
                thisparentdict[noderef].append(childnode.name)
        if len(thisparentdict) > 1 and allleaves == 1:
            thislist = []
            for key in thisparentdict:
                thislist.append(str(key)+" "+str(thisparentdict[key]))
            fp.write(node.name+","+" vs ".join(thislist).replace(",", "_")+"\n")

fp.close()

# label nodes with barcodes for checking output
#ts = TreeStyle()
#ts.show_leaf_name = True
#for leaf in tree.iter_leaves():
#  thisleafcontent = TextFace(" ".join(str(x) for x in list(barcodes[leaf.name])))
#  leaf.add_face(thisleafcontent, 0, "aligned")

# output check tree to file
#fname = str(arg1)+"-check.png"
#tree.render(str(fname), w=800, tree_style=ts)

