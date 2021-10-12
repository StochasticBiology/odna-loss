from Bio.Seq import Seq

aadict = {}
fp = open("Prelims/stats-residue.csv", "r")
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if lineset[0] != "Residue":
        aadict[lineset[0]] = lineset[1:]
    else:
        aalabels = lineset[1:]

fp.close()

codondict = {}
fp = open("Prelims/stats-codon.csv", "r")
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if lineset[0] != "Codon":
        residue = str(Seq(lineset[0]).translate())
        doublelist = [residue, lineset[1:], aadict[residue]]
        codondict[lineset[0]] = [item for sublist in doublelist for item in sublist]
    else:
        codonlabels = lineset[1:]

fp.close()

fp = open("feature-corr.csv", "w")
fp.write("Codon,Residue,"+str(codonlabels).replace("[", "").replace("]", "").replace("'", "").replace(", ", ",")+","+str(aalabels).replace("[", "").replace("]", "").replace("'","").replace(", ", ",")+"\n")
for codon in codondict:
    fp.write(codon+","+str(codondict[codon]).replace("[", "").replace("]", "").replace("'","").replace(", ", ",")+"\n")

fp.close()
