# takes FASTA file and extracts gene statistics information

# needs the BioPython toolbox
from Bio import SeqIO
import re
import sys

# takes nucleotide and protein sequence and returns gene statistics
def getstats(mydna, mypro, codondict, residuedict):

  # populate list of statistics
  nfeatures = len(residuedict['A'])+len(codondict['AAA'])
  nrfeatures = len(residuedict['A'])
  genestats = [0 for i in range(1,nfeatures+1)]

  # first go through protein sequence, looking up by character code
  for code in mypro:
      if code not in residuedict:
          print("Didn't find "+str(code)+"\n")
          return -999
      record = residuedict[code]
      for i in range(0, len(record)):
          genestats[i] = genestats[i] + record[i]

  # now nucleotide sequence, looking up by codon
  if len(mydna) % 3 != 0:
    print(mydna)
    print("Length not a multiple of 3, but ignoring")
#    return -999
  for i in range(0,len(mydna),3):
    if i+2 < len(mydna):
      codon = mydna[i]+mydna[i+1]+mydna[i+2]
      if codon not in codondict:
        print("Didn't find ",+codon+"\n")
        return -999
      record = codondict[codon]
      for i in range(nrfeatures,nrfeatures+len(record)):
          genestats[i] = genestats[i] + record[i-nrfeatures]

  return genestats  

# SeqIO.parse returns a SeqRecord iterator x
# list(x.annotations.keys())
# ['comment', 'source', 'taxonomy', 'keywords', 'references', 'accessions', 'molecule_type', 'data_file_division', 'date', 'organism', 'sequence_version', 'topology']

# this contains features which are type SeqFeature
# dir(x.features[1])
# ['__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__iter__', '__len__', '__module__', '__new__', '__nonzero__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_flip', '_get_location_operator', '_get_ref', '_get_ref_db', '_get_strand', '_set_location_operator', '_set_ref', '_set_ref_db', '_set_strand', '_shift', 'extract', 'id', 'location', 'location_operator', 'qualifiers', 'ref', 'ref_db', 'strand', 'translate', 'type']

# open files for IO
srcfile = sys.argv[1]
outfile = sys.argv[2]
codonstatsfile = sys.argv[3] 
residuestatsfile = sys.argv[4]
label = sys.argv[1].split(".")[0].replace("-",",")

codondict = {}
fp = open(codonstatsfile, "r")
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if lineset[0] == "Codon":
        codondict[lineset[0]] = [entry for entry in lineset[1:]]
    else:
        codondict[lineset[0]] = [float(entry) for entry in lineset[1:]]
        

fp.close()
residuedict = {}
fp = open(residuestatsfile, "r")
for line in fp.readlines():
    lineset = line.rstrip("\n").split(",")
    if lineset[0] == "Residue":
        residuedict[lineset[0]] = [entry for entry in lineset[1:]]
    else:
        residuedict[lineset[0]] = [float(entry) for entry in lineset[1:]]

fp.close()


f1 = open(outfile, "a")


# go through records in the genbank file
for rec in SeqIO.parse(srcfile, "fasta"):
  desc = rec.description

  dnastr = rec.seq
  prostr = rec.seq.translate()

  # get gene stats and produce string
  genestats = getstats(dnastr, prostr, codondict, residuedict)
  if genestats == -999:
    statstr = "statserror"
  else:
    statstr = str(genestats).replace(", ",",").replace("[","").replace("]","")
  # output to files
  f1.write(label+","+desc.replace(",","-")+","+statstr+"\n")
    



      
