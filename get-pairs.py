import sys

infilename = sys.argv[1]
outfilename = sys.argv[2]

fp = open(infilename, "r")
outfile = open(outfilename, "w")
linenum = 1
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if linenum > 1:
        accessionstr = lineset[2]
        accessionstr2 = lineset[4]
        labelstr = lineset[0]
        outfile.write("wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\&id="+accessionstr+"\&rettype=fasta_cds_na -O Downloads/"+labelstr+"-symbiont.fasta\n")
        outfile.write("wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\&id="+accessionstr2+"\&rettype=fasta_cds_na -O Downloads/"+labelstr+"-partner.fasta\n")
    linenum = linenum+1

