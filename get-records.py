import sys

infilename = sys.argv[1]
outfilename = sys.argv[2]

fp = open(infilename, "r")
linenum = 1
accessiondict = {}
for line in fp.readlines():
    # fields in this file include Organism Name,Organism Groups,Strain,BioSample,BioProject,Assembly,Level,Size(Mb),GC%,Replicons,WGS,Scaffolds,CDS,...
    # we grab Organism Name [0], number of CDS records [12], and replicon accessions [9]
    lineset = line.rstrip("\n").split(",")
    if linenum > 1 and float(lineset[12]) > 1000:
        # get list of accessions
        accessions = lineset[9].split(";")
        for accession in accessions:
            species = lineset[0].strip('"')
            # extract both the specific accession code and information about the cellular compartment it is from
            code = accession.split(":")[-1].split("/")[0].strip('"')
            if accession.find("plast") != -1 or accession.find("Plast") != -1:
                organelle = "PT"
            if accession.find("mitoch") != -1 or accession.find("Mitoch") != -1:
                organelle = "MT"
            if accession.find("chromoso") != -1 or accession.find("Chromoso") != -1 or accession.find("linkage gro") != -1 or accession.find("Linkage gro") != -1:
                organelle = "NU"
            print(species+","+organelle+","+code)
            label = species.replace(" ", "_")+"_"+organelle
            # add each accession to dictionary with key species-compartment
            if label in accessiondict:
                accessiondict[label].append(code)
            else:
                accessiondict[label] = [code]
    linenum = linenum + 1

# next create a bash script to download these accessions
# one download command -- and resulting file -- for each species-compartment combination
# the script will include species-compartment labels in the downloaded records via sed
scriptfile = open(outfilename, "w")
scriptfile.write("# automatically generated script to download accessions\n")
for label in accessiondict:
    accessionstr = str(accessiondict[label]).replace("'", "").replace(", ", ",").replace("[", "").replace("]", "")
    scriptfile.write("wget https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore\&id="+accessionstr+"\&rettype=fasta_cds_na -O tmp1\n")
    scriptfile.write("sed 's/>/>"+label+" /g' tmp1 > Downloads/cds-"+label+".txt\n\n")
