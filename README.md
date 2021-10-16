# odna-loss
Bioinformatics pipeline for data acquisition and analysis on oDNA loss across eukaryotes

Requirements
----

Python 3 with `BioPython`, `ete3`; PyMOL; local BLAST; R with libraries:  `ape`, `arm`, `BMA`, `blme`, `caper`, `cowplot`, `e1071`, `geiger`, `GGally`, `ggnewscale`, `ggplot2`, `ggpubr`, `ggpval`, `ggrepel`, `ggtree*`, `ggtreeExtra*`, `glmnet`, `gridExtra`, `hexbin`, `igraph`, `lme4`, `logistf`, `mombf`, `nlme`, `phangorn`, `phytools`, `randomForest`, `stringdist`, `stringr`, `tree` (these R libraries can be installed with a single command-line option below).  (* `ggtree` needs installing via Bioconductor)

A modern machine with 35GB hard disk space should be sufficient for the pipeline without the full BLAST comparison (see below). If you are doing the full BLAST comparison, you'll need at least 160GB hard disk space and several GB of memory, and ideally several cores (default 6).

Outline 
=====

The pipeline consists of (Bash) calls to download data from NCBI, Python and R scripts, and one piece of C code. Broadly, each script solves a subproblem in the workflow -- parsing sequence data, model selection, etc. The Bash script `run-all.sh` serves both as a wrapper for the whole pipeline and to explain the workflow. Use it if you have a Bash environment. If not, open it in a text editor, where you will be able to see which individual scripts are invoked for each part of the pipeline.

The `Prelims/` directory contains some pre-existing scientific data needed for some parts of the analysis, including amino acid properties, eukaryotic accessions, etc.

The scientific aspects of the pipeline come in two parts: "Data curation and production" and "Statistics". There are also "Housekeeping" and "Manuscript" parts for installing software and manuscript preparation. Each of these parts has several submodules. Each submodule can be invoked by passing its name as part of a command-line argument to `run-all.sh`. 

For example,
`./run-all.sh downloadorganelles,processorganelles` would run the first two submodules in the "Data curation and production" part.

`./run-all.sh default` is interpreted as `./run-all.sh downloadorganelles,processorganelles,manuallabel,processtreesmanual,getindicessimple,downloadgenomes,parsegenomes,complexes,downloadotherorganelles,processotherorganelles,datavisualisation,indexregression,bindingenergy,nuclearvsorganelle,otherorganellepredictors,supportingstatistics` and runs the default pipeline without full BLAST analysis.

`./run-all.sh reprocess` is interpreted as `./run-all.sh processorganelles,manuallabel,processtreesmanual,getindicessimple,parsegenomes,complexes,processotherorganelles,datavisualisation,indexregression,bindingenergy,nuclearvsorganelle,otherorganellepredictors,supportingstatistics` and re-runs the default analysis pipeline without full BLAST analysis on pre-downloaded data.

The options are:

* Housekeeping
    * `installrpackages`      -- install the R libraries required
    * `installothers`      -- install other software 
* Data curation and production
    * `downloadorganelles`    -- download organelle genome data
    * `processorganelles`    -- process downloaded organelle genome data to get statistics
    * `fullblast`             -- BLAST all-vs-all organelle gene records
    * `blastdictionary`       -- construct a label replacement dictionary from BLAST analysis
    * `manuallabel`         -- label and process organelle gene records given a manually constructed replacement dictionary
    * `blastlabel`           -- label and process organelle gene records given a BLAST-constructed replacement dictionary
    * `processtreesmanual`  -- process taxonomy trees for manual case
    * `processtreesblast`     -- process taxonomy trees for BLAST case
    * `getindicessimple`     -- use simple barcode analysis to estimate retention indices
    * `downloadgenomes`       -- download whole genome data
    * `parsegenomes`         -- parse downloaded whole genome data
    * `complexes`           -- analyse energetics of organelle protein complexes
    * `downloadotherorganelles`      -- download genomes for symbionts and partners
    * `processotherorganelles`      -- parse symbiont and partner data
* Statistics
    * `datavisualisation`         -- visualise barcodes
    * `indexregression`           -- regression analysis for various retention indices within organelles
    * `bindingenergy`             -- relationship between binding energy and retention
    * `nuclearvsorganelle`        -- nuclear/organelle comparison
    * `otherorganellepredictors`  -- analysis for other organelles
    * `supportingstatistics`       -- correlations, checks, etc
* Manuscript
    * `latextable`             -- convert tabular output to LaTeX form

With a few exceptions, no aspects of the pipeline are very computationally intensive -- each should run in under an hour on a modern machine (most substantially faster). `fullblast` is an optional all-against-all BLAST comparison of all genes in the dataset. This takes some time, even using multiple cores (default 6, set with `-num_threads`). The MT set takes perhaps a day, the PT set perhaps several, and the followup `blastlabel` will likely take an hour or so. Output files are 30GB (MT) and 91GB (PT). Several GB of free memory is also required. `downloadgenomes` involves downloading quite a few whole-genome records from NCBI and can take perhaps half or a full day depending on download rates; `downloadorganelles` and `downloadotherorganelles` also involve getting (fewer) gene records and will take a few minutes at 200kB.

Extendability
----
This fork is attempting to make the addition of new residue and codon features into the analysis as automated as possible. Hopefully only two aspects of the "Data curation and production" part need changing to introduce a new feature: (a) including its values in `Prelims/stats-residues.csv` or `Prelims/stats-codons.csv`; and (b) editing `lengthNormalise.R` to impose its normalisation and give it labels for the "Statistics" part. The individual experiments in the "Statistics" part will need to be edited to ask scientific questions about any new features added; in particular, the model structure in `analysis-model-fit.R`. 
