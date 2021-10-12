# odna-loss
Bioinformatics pipeline for data acquisition and analysis on oDNA loss across eukaryotes

Requirements
----

Python 3 with `BioPython`, `ete3`; PyMOL; local BLAST; R with libraries:  `ape`, `arm`, `blme`, `caper`, `cowplot`, `e1071`, `geiger`, `GGally`, `ggnewscale`, `ggplot2`, `ggpubr`, `ggpval`, `ggrepel`, `ggtree*`, `ggtreeExtra*`, `glmnet`, `gridExtra`, `hexbin`, `igraph`, `lme4`, `logistf`, `mombf`, `nlme`, `phangorn`, `phytools`, `randomForest`, `stringdist`, `stringr`, `tree` (these R libraries can be installed with a single command-line option below).  (* `ggtree` needs installing via Bioconductor)

I haven't computed the exact amount of hard disk space needed, but probably best to have c.10GB free. If you're doing the full BLAST comparison (see below), you'll need at least 130GB hard disk space and several GB of memory, and ideally several cores (default 6).

Outline 
=====

The pipeline consists of Python and R scripts and one piece of C code. Broadly, each script solves a subproblem in the workflow -- parsing sequence data, model selection, etc. The Bash script `run-all.sh` serves both as a wrapper for the whole pipeline and to explain the workflow. Use it if you have a Bash environment. If not, open it in a text editor, where you will be able to see which individual scripts are invoked for each part of the pipeline.

The `Prelims/` directory contains some pre-existing scientific data needed for some parts of the analysis, including amino acid properties, eukaryotic accessions, etc.

The scientific aspects of the pipeline come in two parts: "Data curation and production" and "Statistics". There are also "Housekeeping" and "Manuscript" parts for installing software and manuscript preparation. Each of these parts has several submodules. Each submodule can be invoked by passing its name as part of a command-line argument to `run-all.sh`. 

The options are:

* Housekeeping
    * `installrpackages`      -- install the R libraries required
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
* Statistics
    * `datavisualisation`         -- visualise barcodes
    * `indexregression`           -- regression analysis for various retention indices within organelles
    * `bindingenergy`             -- relationship between binding energy and retention
    * `nuclearvsorganelle`        -- nuclear/organelle comparison
    * `otherorganellepredictors`  -- analysis for other organelles
    * `supportingstatistics`       -- correlations, checks, etc
* Manuscript
    * `latextable`             -- convert tabular output to LaTeX form

With one exception `fullblast`, no aspects of the pipeline are very computationally intensive -- each should run in under an hour on a modern machine (most substantially faster). `fullblast` is an optional all-against-all BLAST comparison of all genes in the dataset. This takes some time, even using multiple cores (default 6, set with `-num_threads`). The MT set takes perhaps a day, the PT set perhaps several, and the followup `blastlabel` will likely take an hour or so. Output files are 30GB (MT) and 91GB (PT). Several GB of free memory is also required.
