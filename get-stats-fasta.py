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
        print("Didn't find "+codon+"\n")
        return -999
      record = codondict[codon]
      for i in range(nrfeatures,nfeatures):
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
refsfile = sys.argv[3]
codonstatsfile = sys.argv[4] 
residuestatsfile = sys.argv[5]

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
f2 = open(refsfile, "a")

# we want to identify genes that are part of particular organelle complexes
# these regexs identify protein labels that describe these complexes
re_ci_nu = re.compile("NADH dehydrogenase|[Uu]biquinone oxidoreductase")
re_cii_nu = re.compile("[Ss]uccinate dehydrogenase|[cC]o[qQ] reductase")
re_ciii_nu = re.compile("[Cc]ytochrome [Bb]|[Cc]ytochrome [Cc] reductase")
re_civ_nu = re.compile("[Cc]ytochrome [cC] oxidase")
re_cv_nu = re.compile("[Aa][Tt][Pp] synthase|ATPase sub")
re_mr_nu = re.compile("[Rr]ibosomal.*[Mm]itochondri")
re_psi_nu = re.compile("[Pp]hotosystem I ")
re_psii_nu = re.compile("[Pp]hotosystem II ")
re_cb_nu = re.compile("[Cc]ytochrome [Bb]6|[Cc]ytochrome f|[Pp]lastocyanin reductase")
re_rb_nu = re.compile("bi.phosphate [Cc]arboxylase")
re_pr_nu = re.compile("[Rr]ibosomal.*[cC]hloroplast")

# list of regexs and labels for their corresponding complexes
re_prots = [re_ci_nu, re_cii_nu, re_ciii_nu, re_civ_nu, re_cv_nu, re_mr_nu, re_psi_nu, re_psii_nu, re_cb_nu, re_rb_nu, re_pr_nu]
prots_lab = ["mt-c1", "mt-c2", "mt-c3", "mt-c4", "mt-c5", "mt-ri", "pt-pa", "pt-pb", "pt-cb", "pt-rb", "pt-ri"]

# several genes match the regexs above but are not within the complex of interested -- perhaps assembly factors, "similar to"s, regulators or otherwise. these negative regexs catch such genes
re_ignore = re.compile("assembly|alternative|containing|dependent|chaperone|kinase|NADH-cytochrome|coupling|maturase|vacuolar|biogenesis|repair|LOW QUALITY PROTEIN|synthetase|activator|reticulum|activase|synthesis|lyase|like| non|transporting|lipid|autoinhibited|membrane|type|required|QUALITY|precursor|inhibitor|proteasomal|proteasome|E1|various|regulatory|Clp|calcium|vesicle|b-245|b5|WRNIP|AAA|Cation|family|remodelling")

print(srcfile)

# go through records in the genbank file
for rec in SeqIO.parse(srcfile, "fasta"):
  desc = rec.description
  for re_prot in re_prots:
    if(re.search(re_prot, desc) != None and re.search(re_ignore, desc) == None):

      dnastr = rec.seq
      prostr = rec.seq.translate()

      # get gene stats and produce string
      genestats = getstats(dnastr, prostr, codondict, residuedict)
      if genestats == -999:
        statstr = "statserror"
      else:
        statstr = str(genestats).replace(", ",",").replace("[","").replace("]","")
        species = rec.name.split("_")[0].lstrip(">")+" "+rec.name.split("_")[1]
        compartment = rec.name.split("_")[-1]
        # output to files
        f1.write(species+","+compartment+","+prots_lab[re_prots.index(re_prot)]+","+statstr+"\n")
    
      f2.write(desc+" | "+prots_lab[re_prots.index(re_prot)]+"\n")


      
