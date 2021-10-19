### master script for organelle genome analysis

# requirements: Python 3 with BioPython, ete3; PyMOL; local BLAST; R (several libraries needed)
# libraries needed: ape, arm, blme, caper, cowplot, e1071, geiger, GGally, ggnewscale, ggplot2, ggpubr, ggpval, ggrepel, ggtree*, ggtreeExtra*, glmnet, gridExtra, hexbin, igraph, lme4, logistf, mombf, nlme, phangorn, phytools, randomForest, stringdist, stringr, tree
# these libraries can be installed with a command-line option below
# * ggtree needs installing via Bioconductor

# takes a command-line argument -- a single string, concatenated with commas or other non-whitespace symbols, determining which aspects of the pipeline will be run
# the options are:

## Housekeeping
# installrpackages      -- install the R libraries required
# installothers         -- install other software

## Data curation and production
# downloadorganelles*           -- download organelle genome data
# processorganelles*            -- process downloaded organelle genome data to get statistics
# fullblast                     -- BLAST all-vs-all organelle gene records
# blastdictionary               -- construct a label replacement dictionary from BLAST analysis
# manuallabel*%                 -- label and process organelle gene records given a manually constructed replacement dictionary
# blastlabel                    -- label and process organelle gene records given a BLAST-constructed replacement dictionary
# processtreesmanual*^          -- process taxonomy trees for manual case
# processtreesblast             -- process taxonomy trees for BLAST case
# getindicessimple*             -- use simple barcode analysis to estimate retention indices
# downloadgenomes*              -- download whole genome data
# parsegenomes*                 -- parse downloaded whole genome data
# complexes*%                   -- analyse energetics of organelle protein complexes
# downloadotherorganelles*      -- download genomes for symbionts and partners
# processotherorganelles*       -- parse symbiont and partner data

## Statistics
# datavisualisation*         -- visualise barcodes
# indexregression*           -- regression analysis for various retention indices within organelles
# bindingenergy*             -- relationship between binding energy and retention
# nuclearvsorganelle*        -- nuclear/organelle comparison
# otherorganellepredictors*  -- analysis for other organelles
# supportingstatistics*      -- correlations, checks, etc

## Manuscript
# latextable             -- convert tabular output to LaTeX form

# if no argument is provided, a default pipeline involving the modules marked with a * is run. this assumes that organelle and whole genome data has been downloaded, and avoids the computationally intensive BLAST dictionary construction.

# a couple of steps in the "Data curation and production" pipeline require manual work. these are the construction of taxonomy trees using NCBI's Common Taxonomy Tree tool. the user should take the list of species produced by the pipeline and download a Phylip tree from the online tool.
# several steps in the "Statistics" pipeline also make use of these taxonomy trees -- by default these use the preconstructed versions.
# * -- part of default pathway
# % -- makes use of already-constructed content
# ^ -- makes optional use of already-constructed content

# the script will store downloaded and processed data in the Data/ subdirectory
# where necessary, use is made of the already-constructed content in the Prelims/ subdirectory
# some outputs from running code are of use in debugging (and science) and are sent to files in the Outputs/ directory

# the I/O for individual scripts are designed to be interpretable from this wrapper script, so that the user can track which files inform which parts of the pipeline. typically a script will take many command-line arguments describing input files and output filenames (or locations if many files are created). parameters are also provided at the command line (see individual notes below)

# process command-line arguments
if [ $# -eq 0 ]
then
    echo "No modules selected! Not running anything. See script preamble for options."
    exit 1
fi

commandstr=$1

if [[ $commandstr == *default* ]]; then
    commandstr="downloadorganelles,processorganelles,manuallabel,processtreesmanual,getindicessimple,downloadgenomes,parsegenomes,complexes,downloadotherorganelles,processotherorganelles,datavisualisation,indexregression,bindingenergy,nuclearvsorganelle,otherorganellepredictors,supportingstatistics,latextable"
fi

if [[ $commandstr == *reprocess* ]]; then
    commandstr="processorganelles,manuallabel,processtreesmanual,getindicessimple,parsegenomes,complexes,processotherorganelles,datavisualisation,indexregression,bindingenergy,nuclearvsorganelle,otherorganellepredictors,supportingstatistics,latextable"
fi

if [[ $commandstr == *statistics* ]]; then
    commandstr="datavisualisation,indexregression,bindingenergy,nuclearvsorganelle,otherorganellepredictors,supportingstatistics,latextable"
fi

if [[ $commandstr == *reduced* ]]; then
    commandstr="processorganelles,manuallabel,processtreesmanual,getindicessimple,parsegenomes,complexes"
fi

echo "Command string is $commandstr"

################
### housekeeping section

# create required directories if not already present
[ ! -d "Data/" ] && mkdir Data/
[ ! -d "Outputs/" ] && mkdir Outputs/
[ ! -d "Plots/" ] && mkdir Plots/
[ ! -d "Downloads/" ] && mkdir Downloads/

if [[ $commandstr == *installrpackages* ]]; then
    echo "Installing packages..."
    R -e "install.packages(c(\"ape\", \"arm\", \"blme\", \"BMA\", \"caper\", \"cowplot\", \"e1071\", \"geiger\", \"GGally\", \"ggnewscale\", \"ggplot2\", \"ggpubr\", \"ggpval\", \"ggrepel\", \"glmnet\", \"gridExtra\", \"hexbin\", \"igraph\", \"lme4\", \"logistf\", \"mombf\", \"nlme\", \"phangorn\", \"phytools\", \"randomForest\", \"stringdist\", \"stringr\", \"tree\"))"
    R -e "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) { install.packages(\"BiocManager\") } ; BiocManager::install(\"ggtree\") ; BiocManager::install(\"ggtreeExtra\")"
fi

if [[ $commandstr == *installothers* ]]; then
    sudo apt install ncbi-blast+
    sudo apt install pip3
    sudo apt install pymol
    pip3 install ete3
    pip3 install biopython
fi
    
################
### organelle genome section
if [[ $commandstr == *downloadorganelles* ]]; then
    echo "Downloading records..."
    # download, extract and concatenate gbff files
    # the download may take several minutes
    # the followup should take a few seconds (as of summer 2020 this is unzipping and concatenating maybe 700MB of data)
    cd Downloads
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.genomic.gbff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.genomic.gbff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.genomic.gbff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.genomic.gbff.gz
    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.3.genomic.gbff.gz

    echo "Unzipping records..."
    gunzip *gz
    cat mitochondrion*gbff > mito.genomic.gbff
    cat plastid*gbff > plastid.genomic.gbff
    cd ..
fi

if [[ $commandstr == *processorganelles* ]]; then
    # convert genomic.gbff files to nucleotide and protein FASTA files
    # the header line for each entry in these files contains the set of gene statistics for use later
    # this should take maybe half an hour
    # this calls a Python script which uses BioPython
    echo "Processing records..."
    python3 get-stats-gb.py Downloads/mito.genomic.gbff Data/mt Prelims/stats-codon.csv Prelims/stats-residue.csv > Outputs/mt-process.txt
    python3 get-stats-gb.py Downloads/plastid.genomic.gbff Data/pt Prelims/stats-codon.csv Prelims/stats-residue.csv > Outputs/pt-process.txt
    python3 get-raw-counts.py Data/mt-dna.fasta Data/mt-gene-raw-occurrence.csv
    python3 get-raw-counts.py Data/pt-dna.fasta Data/pt-gene-raw-occurrence.csv
fi

if [[ $commandstr == *fullblast* ]]; then
    echo "Performing full BLAST..."
    # BLAST these nucleotide files to get labelling clusters
    # this takes some time, even using multiple cores (default 6, set with -num_threads). MT set perhaps a day, PT set perhaps several. Output files are 30gb MT and 91gb PT. 
    # a good bit of free memory (several GB) is also required
    # requires local BLAST functionality (blastx)
    makeblastdb -in Data/mt-pro.fasta -dbtype prot
    blastx -query Data/mt-dna.fasta -db Data/mt-pro.fasta -num_threads 6 -outfmt "6 qseqid sseqid qlen slen length pident evalue bitscore" > Data/mt-pro-blast.txt
    makeblastdb -in Data/pt-pro.fasta -dbtype prot
    blastx -query Data/pt-dna.fasta -db Data/pt-pro.fasta -num_threads 6 -outfmt "6 qseqid sseqid qlen slen length pident evalue bitscore" > Data/pt-pro-blast.txt
fi

if [[ $commandstr == *blastdictionary* ]]; then
    echo "Labelling full BLAST output..."
    # process BLAST output to create replacement label list
    # this code takes two parameters. the first is a threshold applied to a score based on BLAST statistics, determining what is a "good hit" (0.333 is perfect)
    # the second is an "intersection" parameter. briefly, if a proportion x of appearances of gene1 occur as a "good hit" with gene2, treat them as identical
    gcc -o3 labeller.c -lm -o labeller.ce
    ./labeller.ce Data/mt-gene-raw-occurrence.csv Data/mt-pro-blast.txt Data/mt-blast-replace.csv 0.25 0.25 > Outputs/mt-blast-dictionary.txt
    ./labeller.ce Data/pt-gene-raw-occurrence.csv Data/pt-pro-blast.txt Data/pt-blast-replace.csv 0.25 0.25 > Outputs/pt-blast-dictionary.txt
fi

if [[ $commandstr == *manuallabel* ]]; then
    echo "Assigning manual labels..."
    # label individual gene data and produce species barcodes
    # using manual replacement list
    # 10 is a threshold for inclusion (a gene must be present in >=10 species): a parameter of the pipeline
    python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv  Species,Compartment,GeneLabel, > Data/mt-stats-manual.csv
    python3 process-labels.py 10 Data/mt-dna.fasta Prelims/mt-manual-replace.csv Data/mt-stats-manual.csv Data/mt-barcodes-manual.csv Data/mt-species-manual.txt Data/mt-gene-occurrence-manual.csv
    python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv  Species,Compartment,GeneLabel, > Data/pt-stats-manual.csv
    python3 process-labels.py 10 Data/pt-dna.fasta Prelims/pt-manual-replace.csv Data/pt-stats-manual.csv Data/pt-barcodes-manual.csv Data/pt-species-manual.txt Data/pt-gene-occurrence-manual.csv 
fi

if [[ $commandstr == *blastlabel* ]]; then
    echo "Assigning BLAST labels..."
    # as above, using BLAST-derived replacement list
    # some cleaning from the MT BLAST output goes on here. the seds replace "nd", the most common form for some "nad" genes, with "nad". the Python script ignores gene labels containing "-i" or "oi" -- some rare-ish isoforms have these.
    python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv Species,Compartment,GeneLabel, > Data/mt-stats-blast.csv
    python3 process-labels.py 10 Data/mt-dna.fasta Data/mt-blast-replace.csv Data/mt-stats-blast.csv Data/mt-barcodes-blast.csv Data/mt-species-blast.txt Data/mt-gene-occurrence-blast.csv
    sed -i 's/,nd/,nad/g' Data/mt-barcodes-blast.csv
    sed -i 's/,nd/,nad/g' Data/mt-stats-blast.csv
    sed -i 's/nd/nad/g' Data/mt-gene-occurrence-blast.csv
    python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv Species,Compartment,GeneLabel, > Data/pt-stats-blast.csv
    python3 process-labels.py 10 Data/pt-dna.fasta Data/pt-blast-replace.csv Data/pt-stats-blast.csv Data/pt-barcodes-blast.csv Data/pt-species-blast.txt Data/pt-gene-occurrence-blast.csv 
fi

## manual step!
# now use Taxonomy Common Tree to get phylip trees, using mt-species-[manual/blast].txt and pt-species-[manual/blast].txt
# put these in the Downloads/ subdirectory as e.g. mt-tree-manual.phy
# if this step isn't taken, preconstructed trees will be used

if [[ $commandstr == *processtreesmanual* ]]; then
    echo "Processing trees with manual data..."
    # embed species barcodes on taxonomy, and thus compute transitions between barcodes for use in evolutionary inference
    if [ -f Downloads/mt-tree-manual.phy ]; then
	python3 tree-parse.py Data/mt-barcodes-manual.csv Downloads/mt-tree-manual.phy Data/mt-trans-manual.csv Data/mt-sisters-manual.csv > Outputs/mt-manual-out.txt 
    else
	echo "Using preconstructed MT tree"
        python3 tree-parse.py Data/mt-barcodes-manual.csv Prelims/mt-tree-manual.phy Data/mt-trans-manual.csv Data/mt-sisters-manual.csv > Outputs/mt-manual-out.txt 
    fi
    if [ -f Downloads/pt-tree-manual.phy ]; then
	python3 tree-parse.py Data/pt-barcodes-manual.csv Downloads/pt-tree-manual.phy Data/pt-trans-manual.csv Data/pt-sisters-manual.csv > Outputs/pt-manual-out.txt 
    else
        echo "Using preconstructed PT tree"
        python3 tree-parse.py Data/pt-barcodes-manual.csv Prelims/pt-tree-manual.phy Data/pt-trans-manual.csv Data/pt-sisters-manual.csv > Outputs/pt-manual-out.txt 
    fi    
fi	   

if [[ $commandstr == *processtreesblast* ]]; then
    echo "Processing trees with BLAST data..."
    python3 tree-parse.py Data/mt-barcodes-blast.csv Prelims/mt-tree-manual.phy Data/mt-trans-blast.csv Data/mt-sisters-blast.csv > Outputs/mt-blast-out.txt
    python3 tree-parse.py Data/pt-barcodes-blast.csv Prelims/pt-tree-manual.phy Data/pt-trans-blast.csv Data/pt-sisters-blast.csv > Outputs/pt-blast-out.txt
fi

if [[ $commandstr == *getindicessimple* ]]; then
    echo "Computing retention indices..."
    Rscript --vanilla get-indices-transitions.R 0 Data/mt-trans-manual.csv Data/mt-simple-manual-indices.csv
    Rscript --vanilla get-indices-transitions.R 0 Data/pt-trans-manual.csv Data/pt-simple-manual-indices.csv
    # indices not computed for BLAST stats here because of the tight correlation between BLAST and manual protocols
#    Rscript --vanilla get-indices-transitions.R 0 Data/mt-trans-blast.csv Data/mt-simple-blast-indices.csv
#    Rscript --vanilla get-indices-transitions.R 0 Data/pt-trans-blast.csv Data/pt-simple-blast-indices.csv
    Rscript --vanilla get-indices-transitions.R 1 Data/pt-trans-manual.csv Prelims/pt-tree-manual.phy Data/pt-simple-manual-green-indices.csv Data/pt-simple-manual-red-indices.csv

    Rscript --vanilla get-indices-barcodes.R Data/mt-barcodes-manual.csv Data/mt-barcode-manual-indices.csv
    Rscript --vanilla get-indices-barcodes.R Data/pt-barcodes-manual.csv Data/pt-barcode-manual-indices.csv
fi

################
### whole genome section

if [[ $commandstr == *downloadgenomes* ]]; then
    echo "Downloading genome records..."
    python3 get-records.py Prelims/eukaryotes.csv ./download.sh
    # now execute these Entrez calls to get the data
    chmod +x ./download.sh
    ./download.sh
fi

if [[ $commandstr == *parsegenomes* ]]; then
    echo "Parsing genome records..."
    python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv Species,Compartment,GeneLabel, > Data/all-stats.csv
#    echo Species,Compartment,GeneLabel,Length,Hydro,Hydro_i,MolWeight,pKa1,pKa2,A_Glu,CW,GC,Uni1,Uni2,Robust,GC12,GC3 > Data/all-stats.csv
    echo > Data/all-refs.txt
    # parse the resulting datafiles to extract quantitative data
    for file in Downloads/cds*txt
    do
	python3 get-stats-fasta.py $file Data/all-stats.csv Data/all-refs.txt Prelims/stats-codon.csv Prelims/stats-residue.csv
    done
fi

################
### protein complex section

if [[ $commandstr == *complexes* ]]; then
    echo "Processing complex records..."
    COMPLEXES=('1oco' '2h88' '5iu0' '5mlc' '5xte' '6fkf' '1q90' '2wsc' '5mdx' '5o31' '6cp3')
    # respectively CIV, CII, rubisco, chlororibo, CIII, chloro-atp, cb6f, PSI, PSII, CI, CV
    
    echo PDBLabel,GeneLabel,Stoichiometry,EnergywLigand,InterfaceswLigand,EnergywoLigand,InterfaceswoLigand,MeanEnergywLigand > Data/complex-data.csv

    for complex in "${COMPLEXES[@]}"
    do
	echo $complex
	python3 parse-interfaces.py $complex Prelims/$complex-lig.txt Prelims/$complex-pdb.html Prelims/ Data/mt-gene-occurrence-manual.csv Data/pt-gene-occurrence-manual.csv Data/complex-data.csv plot-$complex.py Plots/plot-$complex.png > Outputs/$complex-tmp.txt
	pymol plot-$complex.py
    done
fi

if [[ $commandstr == *downloadotherorganelles* ]]; then
    echo "Downloading symbiont records..."
    python3 get-pairs.py Prelims/symbionts.csv ./download-pairs.sh
    chmod +x download-pairs.sh
    ./download-pairs.sh
fi

if [[ $commandstr == *processotherorganelles* ]]; then
    python3 get-feature-labels.py Prelims/stats-residue.csv Prelims/stats-codon.csv  SystemLabel,Partner,GeneLabel, > Data/symbiont-all-stats.csv
   
    # parse the resulting datafiles to extract quantitative data
    for file in Downloads/*-symbiont.fasta
    do
	python3 get-stats-symbiont.py $file Data/symbiont-all-stats.csv Prelims/stats-codon.csv Prelims/stats-residue.csv 
    done
    for file in Downloads/*-partner.fasta
    do
	python3 get-stats-symbiont.py $file Data/symbiont-all-stats.csv Prelims/stats-codon.csv Prelims/stats-residue.csv 
    done
fi

##############
### analysis and plotting code

if [[ $commandstr == *indexregression* ]]; then
    echo "Analysing -- index regression..."
    # first produce summaries of gene statistics according to different sampling protocols
    Rscript analysis-summarise-stats.R Data/mt-stats-manual.csv Data/pt-stats-manual.csv Prelims/mt-tree-manual.phy Prelims/pt-tree-manual.phy Data/mt-stats-means-manual.csv Data/pt-stats-means-manual.csv Plots/average-stats.png

    # now do Bayesian model selection with linear model for different sampling protocols
    Rscript analysis-model-fit.R Data/mt-stats-means-manual.csv Data/mt-gene-occurrence-manual.csv Data/mt-simple-manual-indices.csv Data/pt-stats-means-manual.csv Data/pt-gene-occurrence-manual.csv Data/pt-simple-manual-indices.csv Sampled Plots/ Data/ simple

    Rscript analysis-model-fit.R Data/mt-stats-means-manual.csv Data/mt-gene-occurrence-manual.csv Data/mt-simple-manual-indices.csv Data/pt-stats-means-manual.csv Data/pt-gene-occurrence-manual.csv Data/pt-simple-manual-indices.csv AltSampled Plots/ Data/ simple-alt

    Rscript analysis-model-fit.R Data/mt-stats-means-manual.csv Data/mt-gene-occurrence-manual.csv Data/mt-barcode-manual-indices.csv Data/pt-stats-means-manual.csv Data/pt-gene-occurrence-manual.csv Data/pt-barcode-manual-indices.csv Sampled Plots/ Data/ barcode

    # now do alternative models for retention index regression
    # the final two parameters here are the test (as opposed to training) proportion for the dataset split, and the number of different training/test sets to explore for each method
    # simple-manual-sampled
    Rscript analysis-index-regression-alt.R Data/mt-stats-means-manual.csv Data/mt-gene-occurrence-manual.csv Data/mt-simple-manual-indices.csv Data/pt-stats-means-manual.csv Data/pt-gene-occurrence-manual.csv Data/pt-simple-manual-indices.csv Sampled simple-manual-sampled Data/ Plots/ 0.5 100
    # simple-blast-sampled -- not done here
    # Rscript analysis-index-regression-alt.R Data/mt-stats-means-blast.csv Data/mt-gene-occurrence-blast.csv Data/mt-simple-blast-indices.csv Data/pt-stats-means-blast.csv Data/pt-gene-occurrence-blast.csv Data/pt-simple-blast-indices.csv Sampled simple-blast-sampled Data/ Plots/ 0.5 100
fi

if [[ $commandstr == *datavisualisation* ]]; then
    echo "Analysing -- data visualisation..."
    Rscript analysis-data-visualisation.R Data/mt-barcodes-manual.csv Prelims/mt-tree-manual.phy Data/mt-stats-means-manual.csv Data/mt-simple-manual-indices.csv Data/pt-barcodes-manual.csv Prelims/pt-tree-manual.phy Data/pt-stats-means-manual.csv Data/pt-simple-manual-indices.csv Plots/
fi

if [[ $commandstr == *bindingenergy* ]]; then
    echo "Analysing -- complex energies..."
    # the final parameter here is a threshold for occurrences, used when we binarise genes into "highly" vs "not highly" retained in oDNA
    Rscript analysis-complexes.R Data/complex-data.csv Data/mt-gene-occurrence-manual.csv Data/pt-gene-occurrence-manual.csv Prelims/label-aliases.csv Prelims/complex-labels.csv Plots/ Data/ 1000
    Rscript analysis-complexes.R Data/complex-data.csv Data/mt-gene-occurrence-manual.csv Data/pt-gene-occurrence-manual.csv Prelims/label-aliases.csv Prelims/complex-labels.csv Plots/ Data/ 0
    Rscript analysis-complexes-alt.R Data/complex-data.csv Data/mt-gene-occurrence-manual.csv Data/pt-gene-occurrence-manual.csv Prelims/label-aliases.csv Prelims/complex-labels.csv Plots/complexes-alt-0.png Data/complex-model-summaries-alt-0.txt 0
    Rscript analysis-complexes-alt.R Data/complex-data.csv Data/mt-gene-occurrence-manual.csv Data/pt-gene-occurrence-manual.csv Prelims/label-aliases.csv Prelims/complex-labels.csv Plots/complexes-alt-1000.png Data/complex-model-summaries-alt-1000.txt 1000
fi

if [[ $commandstr == *nuclearvsorganelle* ]]; then
    echo "Analysing -- encoding compartment..."
    Rscript analysis-nuc-org-glm.R Data/mt-stats-manual.csv Data/pt-stats-manual.csv Data/all-stats.csv manual Data/ Plots/ 0.5 10
    Rscript analysis-nuc-org-alt.R Data/mt-stats-manual.csv Data/pt-stats-manual.csv Data/all-stats.csv manual Data/ Plots/ 0.5 10
    Rscript analysis-phylo-plot.R Prelims/whole-genome-species.phy Data/all-stats.csv Data/mt-stats-manual.csv Data/pt-stats-manual.csv Plots/phy-plot
fi

if [[ $commandstr == *otherorganellepredictors* ]]; then
    echo "Analysing -- other symbionts..."
    Rscript analysis-other-symbionts.R Data/symbiont-all-stats.csv Plots/
fi

if [[ $commandstr == *supportingstatistics* ]]; then
    echo "Analysing -- supporting statistics..."
    # correlations between genetic features
    python3 features-corr.py
    R CMD BATCH analysis-feature-corr.R
    # links between features and energy centrality 
    R CMD BATCH analysis-energy-hydro.R
    # manual/BLAST correlation
    Rscript analysis-si-plots.R
fi

##############
### manuscript preparation
if [[ $commandstr == *latextable* ]]; then
    echo "Reformatting tables..."
    chmod +x latex-table.sh
    # these just convert text output from R into LaTeX-formatted tables
    ./latex-table.sh Data/model-sel-simple-bayeslm-stats.csv > Data/latex-modelsel-simple.tex
    ./latex-table.sh Data/model-sel-barcode-bayeslm-stats.csv > Data/latex-modelsel-barcode.tex
    ./latex-table.sh Data/index-regression-simple-manual-sampled.csv > Data/latex-regression-simple.tex
    ./latex-table.sh Data/nuc-org-manual-results.csv > Data/latex-nucorg.tex
    # this pulls statistics from the code output and embeds in LaTeX
    Rscript latex-interpret.R
fi
