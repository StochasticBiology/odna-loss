#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3) {
  stop("Need mode (0, normal; 1, split Rhodophyta), and filenames")
}

# if mode = 0, arguments are: 0 ; infile ; outfile
# if mode = 1, arguments are: 1 ; infile ; treefile ; outfilegreen ; outfilered

df = read.csv(args[2], header=T, stringsAsFactor=F)
ngenes = (ncol(df)-2)/2

# construct summed gene counts of ancestor and descendant state for each transition
df$a.sum = rowSums(df[,3:(3+ngenes-1)])
df$d.sum = rowSums(df[,(3+ngenes):(3+2*ngenes-1)])

if(args[1] == 0) {
  gene.list = data.frame(GeneLabel=NULL, Index=NULL)
  # loop through genes
  for(i in 1:ngenes) {
    # find those transitions that involve the loss of gene i
    refs = which(df[,i+2] == 1 & df[,i+ngenes+2] == 0)
    # grab the average of the ancestor and descendant gene sums for these transitions
    # higher values mean we lose this gene while retaining more others
    gene.set = (df$a.sum[refs]+df$d.sum[refs])/2
    # invert and normalise
    gene.set = (ngenes-gene.set)/ngenes
    gene.list = rbind(gene.list, data.frame(GeneLabel = gsub("a[.]", "", colnames(df)[2+i]), Index = mean(gene.set)))
  }
  write.csv(gene.list, args[3], row.names = FALSE, quote= FALSE)
} else {
  library(phytools)

  treefilename = args[3]
  # read taxonomy tree
  treeString = tolower(paste(readLines(treefilename), collapse=""))
  treeString = gsub(" ", "_", treeString)
  pt.tree <- read.tree(text=treeString)

  all.labels = c(pt.tree$tip.label, pt.tree$node.label)
  pt.tree.labels = tolower(gsub("_", " ", all.labels))

  # numbers of leaves and internal nodes
  n.tips = length(pt.tree$tip.label)
  n.nodes = length(pt.tree$node.label)

  rhodophytes = pt.tree.labels[getDescendants(pt.tree, which(pt.tree.labels == "rhodophyta"))]
  
  a.rhodophyte.refs = which(sapply(df$Ancestor, function(y) y %in% rhodophytes))
  d.rhodophyte.refs = which(sapply(df$Descendant, function(y) y %in% rhodophytes))
  rhodophyte.refs = c(a.rhodophyte.refs, d.rhodophyte.refs)
  pt.green = df[-rhodophyte.refs,]
  pt.red = df[rhodophyte.refs,]

  # first green
  gene.list = data.frame(GeneLabel=NULL, Index=NULL)
  for(i in 1:ngenes) {
    refs = which(pt.green[,i+2] == 1 & pt.green[,i+ngenes+2] == 0)
    gene.set = (pt.green$a.sum[refs]+pt.green$d.sum[refs])/2
    gene.set = (ngenes-gene.set)/ngenes
    gene.list = rbind(gene.list, data.frame(GeneLabel = gsub("a[.]", "", colnames(pt.green)[2+i]), Index = mean(gene.set)))
  }
  write.csv(gene.list, args[4], row.names = FALSE, quote= FALSE)

  # now red
  gene.list = data.frame(GeneLabel=NULL, Index=NULL)
  for(i in 1:ngenes) {
    refs = which(pt.red[,i+2] == 1 & pt.red[,i+ngenes+2] == 0)
    gene.set = (pt.red$a.sum[refs]+pt.red$d.sum[refs])/2
    gene.set = (ngenes-gene.set)/ngenes
    gene.list = rbind(gene.list, data.frame(GeneLabel = gsub("a[.]", "", colnames(pt.red)[2+i]), Index = mean(gene.set)))
  }
  write.csv(gene.list, args[5], row.names = FALSE, quote= FALSE)
}

  
