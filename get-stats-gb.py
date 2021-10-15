# takes GenBank files and extracts annotation and gene statistics information

# needs the BioPython toolbox
from Bio import SeqIO
import sys

# open files for IO
srcfile = sys.argv[1]
outfile = sys.argv[2]
codonstatsfile = sys.argv[3] 
residuestatsfile = sys.argv[4]


# takes nucleotide and protein sequence and returns gene statistics
def getstats(mydna, mypro, codondict, residuedict, permissive):

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
  if permissive == False and len(mydna) % 3 != 0:
#    print(mydna)
#    print("Length not a multiple of 3")
    return -999
  for i in range(0,len(mydna),3):
    if i+2 < len(mydna):
      codon = mydna[i]+mydna[i+1]+mydna[i+2]
      if codon not in codondict:
        print("Didn't find "+codon+"\n")
        return -999
      record = codondict[codon]
      for i in range(nrfeatures,nfeatures+1):
          genestats[i] = genestats[i] + record[i-nrfeatures]

  return genestats  

# SeqIO.parse returns a SeqRecord iterator x
# list(x.annotations.keys())
# ['comment', 'source', 'taxonomy', 'keywords', 'references', 'accessions', 'molecule_type', 'data_file_division', 'date', 'organism', 'sequence_version', 'topology']

# this contains features which are type SeqFeature
# dir(x.features[1])
# ['__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__iter__', '__len__', '__module__', '__new__', '__nonzero__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_flip', '_get_location_operator', '_get_ref', '_get_ref_db', '_get_strand', '_set_location_operator', '_set_ref', '_set_ref_db', '_set_strand', '_shift', 'extract', 'id', 'location', 'location_operator', 'qualifiers', 'ref', 'ref_db', 'strand', 'translate', 'type']


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


f1 = open(outfile+"-dna.fasta", "w")
f2 = open(outfile+"-pro.fasta", "w")

counter = 0
# go through records in the genbank file
for rec in SeqIO.parse(srcfile, "genbank"):
  print(counter)
  counter = counter+1
  if rec.features:
    # we've found a suitable record -- grab the organism name
    if "organism" in rec.annotations:
      orgname = rec.annotations["organism"].replace(" ", "_").replace(",", "_")
    else:
      orgname = "anon"
    print(orgname)

    # loop through the genetic features in this record
    for feature in rec.features:
      # if this is a CDS, grab gene name and any other info if present
      if feature.type == "CDS":
        if "gene" in feature.qualifiers:
          genename = feature.qualifiers["gene"][0].replace(" ", "_").replace(",", "")
        else:
          genename = "anon"
        if "db_xref" in feature.qualifiers:
          genelabel1 = feature.qualifiers["db_xref"][0].replace(" ", "_").replace(",", "")
        else:
          genelabel1 = "anon"
        if "locus_tag" in feature.qualifiers:
          genelabel2 = feature.qualifiers["locus_tag"][0].replace(" ", "_").replace(",", "")
        else:
          genelabel2 = "anon"
        if "transl_table" in feature.qualifiers:
          ttable = int(feature.qualifiers["transl_table"][0])
        else:
          ttable = 1
        if "transl_except" in feature.qualifiers:
          permissive = True
        else:
          permissive = False
          
        try:
          # produce strings of nucleotide and translated protein sequence
          dnastr = str(feature.location.extract(rec).seq)
          if "translation" in feature.qualifiers:
            prostr = feature.qualifiers['translation'][0]
          else:
            prostr = str(feature.location.extract(rec).seq.translate(table=ttable))

          # get gene stats and produce string
          genestats = getstats(dnastr, prostr, codondict, residuedict, permissive)
          if genestats == -999:
            print("debugging: statserror"+str(feature))
            statstr = "statserror"
          else:
            statstr = str(genestats).replace(", ",",").replace("[","").replace("]","")

          # output to files
          f1.write(">"+str(orgname).lower()+","+"".join(genename).lower()+","+"".join(genelabel1)+","+"".join(genelabel2)+","+statstr+"\n")
          f2.write(">"+str(orgname).lower()+","+"".join(genename).lower()+","+"".join(genelabel1)+","+"".join(genelabel2)+","+statstr+"\n")
          f1.write(dnastr+"\n")
          f2.write(prostr+"\n")
        except:
          print("Sequence reading error")
