#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 7) {
  stop("Need: MT stats file, PT stats file, MT tree, PT tree, MT output, PT output, output plot")
}

#args = c("Data/mt-stats-manual.csv", "Data/pt-stats-manual.csv", "Prelims/mt-tree-manual.phy", "Prelims/pt-tree-manual.phy", "Data/mt-stats-means-manual.csv", "Data/pt-stats-means-manual.csv", "Plots/average-stats.png")

message("Loading libraries...")

library(ggplot2)
library(gridExtra)

library(phytools)
set.seed(121)

source("lengthNormalise.R")

mt.stats.file = args[1]
pt.stats.file = args[2]
# tree files looped through below
mt.out.file = args[5]
pt.out.file = args[6]
out.plot = args[7]

message("Loading data...")

mt = read.csv(mt.stats.file, header=T, stringsAsFactor=F)
pt = read.csv(pt.stats.file, header=T, stringsAsFactor=F)

# implementation note -- galaxaura rugosa has been included rather than chondrus crispus here because, as of Oct 2021, the latter seems to be absent from the NCBI PT dataset

#sample.set.1 = c("homo sapiens", "reclinomonas americana", "arabidopsis thaliana", "plasmodium falciparum", "saccharomyces cerevisiae", "chondrus crispus")
sample.set.1 = c("homo sapiens", "reclinomonas americana", "arabidopsis thaliana", "plasmodium falciparum", "saccharomyces cerevisiae", "galaxaura rugosa")

#sample.set.2 = c("chondrus crispus", "reclinomonas americana")
sample.set.2 = c("galaxaura rugosa", "reclinomonas americana")
for(j in c(3,4)) {
  tree.file = args[j]
  treeString = tolower(paste(readLines(tree.file), collapse=""))
  treeString = gsub(" ", "_", treeString)
  tree <- read.tree(text=treeString)

  # list of all node labels
  all.labels = c(tree$tip.label, tree$node.label)

  # numbers of leaves and internal nodes
  n.tips = length(tree$tip.label)
  n.nodes = length(tree$node.label)

  # identify tree root
  luca.node = grep("eukaryota", all.labels)

  # identify those clades that are directly descended from root...
  clades.refs = tree$edge[which(tree$edge[,1]==luca.node),2]

  # ... and which are not leaves
  clades.refs = clades.refs[clades.refs > n.tips]
  for(clade.ref in clades.refs) {
    des = getDescendants(tree, clade.ref)
    leaves = des[des <= n.tips]
    species = sample(leaves, 1)
    sample.set.2 = c(sample.set.2, all.labels[species])
  }
}
sample.set.2 = gsub("_", " ", sample.set.2)


message("Computing summaries...")

mt.df = mt[1,]
mt.df = mt.df[-1,]
mt.df$Protocol = NULL
mt.df$Species = mt.df$Compartment = NULL
for(genelabel in unique(mt$GeneLabel)) {
  subset = mt[mt$GeneLabel == genelabel,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  genestats$Protocol = "All"
  mt.df = rbind(mt.df, as.data.frame(genestats, stringsAsFactors=F))

  subset = mt[mt$GeneLabel == genelabel & mt$Species %in% sample.set.1,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  genestats$Protocol = "Sampled"
  mt.df = rbind(mt.df, as.data.frame(genestats, stringsAsFactors=F))

  subset = mt[mt$GeneLabel == genelabel & mt$Species %in% sample.set.2,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  genestats$Protocol = "AltSampled"
  mt.df = rbind(mt.df, as.data.frame(genestats, stringsAsFactors=F))
}
mt.df = lengthNormalise(mt.df)


pt.df = pt[1,]
pt.df = pt.df[-1,]
pt.df$Protocol = NULL
pt.df$Species = pt.df$Compartment = NULL
for(genelabel in unique(pt$GeneLabel)) {
  subset = pt[pt$GeneLabel == genelabel,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  genestats$Protocol = "All"
  pt.df = rbind(pt.df, as.data.frame(genestats, stringsAsFactors=F))

  subset = pt[pt$GeneLabel == genelabel & pt$Species %in% sample.set.1,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  genestats$Protocol = "Sampled"
  pt.df = rbind(pt.df, as.data.frame(genestats, stringsAsFactors=F))

  subset = pt[pt$GeneLabel == genelabel & pt$Species %in% sample.set.2,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  genestats$Protocol = "AltSampled"
  pt.df = rbind(pt.df, as.data.frame(genestats, stringsAsFactors=F))
}
pt.df = lengthNormalise(pt.df)

message("Writing...")

write.csv(mt.df, mt.out.file, row.names=F, quote=F)
write.csv(pt.df, pt.out.file, row.names=F, quote=F)

n.features = 9
#png(out.plot, width=1200,height=800)
#par(mfrow=c(4,n.features))
plots = list()
for(i in 1:n.features) {
  to.plot = data.frame(GeneLabel=NULL, Sampled=NULL, AltSampled=NULL)
  for(j in 1:length(unique(mt.df$GeneLabel))) {
    label = mt.df$GeneLabel[j]
    ref.1 = which(mt.df$GeneLabel == label & mt.df$Protocol == "Sampled")
    ref.2 = which(mt.df$GeneLabel == label & mt.df$Protocol == "AltSampled")
    if(length(ref.1) > 0 & length(ref.2) > 0) {
      to.plot = rbind(to.plot, data.frame(GeneLabel=label, Sampled=mt.df[ref.1,i], AltSampled=mt.df[ref.2,i]))
    }
  }
  plots[[length(plots)+1]] = ggplot(to.plot, aes(x=Sampled, y=AltSampled)) +
                               geom_point() +
			       geom_abline(slope=1, intercept=0) +
			       ggtitle(paste("MT", colnames(mt.df)[i], "\nr=", round(cor(to.plot$Sampled, to.plot$AltSampled, use = "complete.obs"), digits=3))) +
			       theme_light()
}
for(i in 1:n.features) {
  to.plot = data.frame(GeneLabel=NULL, Sampled=NULL, AltSampled=NULL)
  for(j in 1:length(unique(pt.df$GeneLabel))) {
    label = pt.df$GeneLabel[j]
    ref.1 = which(pt.df$GeneLabel == label & pt.df$Protocol == "Sampled")
    ref.2 = which(pt.df$GeneLabel == label & pt.df$Protocol == "AltSampled")
    if(length(ref.1) > 0 & length(ref.2) > 0) {
      to.plot = rbind(to.plot, data.frame(GeneLabel=label, Sampled=pt.df[ref.1,i], AltSampled=pt.df[ref.2,i]))
    }
  }
  plots[[length(plots)+1]] = ggplot(to.plot, aes(x=Sampled, y=AltSampled)) +
                               geom_point() +
			       geom_abline(slope=1, intercept=0) +
			       ggtitle(paste("PT", colnames(mt.df)[i], "\nr=", round(cor(to.plot$Sampled, to.plot$AltSampled, use = "complete.obs"), digits=3))) +
			       theme_light()
}

for(i in 1:n.features) {
  to.plot = data.frame(GeneLabel=NULL, Sampled=NULL, All=NULL)
  for(j in 1:length(unique(mt.df$GeneLabel))) {
    label = mt.df$GeneLabel[j]
    ref.1 = which(mt.df$GeneLabel == label & mt.df$Protocol == "Sampled")
    ref.2 = which(mt.df$GeneLabel == label & mt.df$Protocol == "All")
    if(length(ref.1) > 0 & length(ref.2) > 0) {
      to.plot = rbind(to.plot, data.frame(GeneLabel=label, Sampled=mt.df[ref.1,i], All=mt.df[ref.2,i]))
    }
  }
  plots[[length(plots)+1]] = ggplot(to.plot, aes(x=Sampled, y=All)) +
                               geom_point() +
			       geom_abline(slope=1, intercept=0) +
			       ggtitle(paste("MT", colnames(mt.df)[i], "\nr=", round(cor(to.plot$Sampled, to.plot$All, use = "complete.obs"), digits=3))) +
			       theme_light()
}
for(i in 1:n.features) {
  to.plot = data.frame(GeneLabel=NULL, Sampled=NULL, All=NULL)
  for(j in 1:length(unique(pt.df$GeneLabel))) {
    label = pt.df$GeneLabel[j]
    ref.1 = which(pt.df$GeneLabel == label & pt.df$Protocol == "Sampled")
    ref.2 = which(pt.df$GeneLabel == label & pt.df$Protocol == "All")
    if(length(ref.1) > 0 & length(ref.2) > 0) {
      to.plot = rbind(to.plot, data.frame(GeneLabel=label, Sampled=pt.df[ref.1,i], All=pt.df[ref.2,i]))
    }
  }
  plots[[length(plots)+1]] = ggplot(to.plot, aes(x=Sampled, y=All)) +
                               geom_point() +
			       geom_abline(slope=1, intercept=0) +
			       ggtitle(paste("PT", colnames(mt.df)[i], "\nr=", round(cor(to.plot$Sampled, to.plot$All, use = "complete.obs"), digits=3))) +
			       theme_light()
}

res.factor = 3
png(out.plot, width=1200*res.factor,height=800*res.factor,res=72*res.factor)
grid.arrange(grobs=plots, nrow=4)
dev.off()
