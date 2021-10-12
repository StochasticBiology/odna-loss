#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 9) {
  stop("Need 9 argumnets -- MT barcodes, tree, stats, indices; PT barcodes, tree, stats, indices; output directory")
}

# args = c("Data/mt-barcodes-manual.csv", "Prelims/mt-tree-manual.phy", "Data/mt-stats-means-manual.csv", "Data/mt-simple-manual-indices.csv", "Data/pt-barcodes-manual.csv", "Prelims/pt-tree-manual.phy", "Data/pt-stats-means-manual.csv", "Data/pt-simple-manual-indices.csv", "Plots/")


message("Loading libraries...")

library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(gridExtra)
library(phytools)
library(igraph)
library(cowplot)

mtbarcodefilename = args[1]
mttreefilename = args[2]
mtgenestatsfilename = args[3]
mtindexfilename = args[4]

ptbarcodefilename = args[5]
pttreefilename = args[6]
ptgenestatsfilename = args[7]
ptindexfilename = args[8]

out.dir = args[9]


#################### MT

#### get barcode data

message("Getting data...")

# read barcode data
mt.df <- read.csv(mtbarcodefilename, header=T)
mt.ngenes = ncol(mt.df)-1

# get species list
mt.species.list = trimws(mt.df$Species)
# create new data frame, storing text versions of barcodes and total gene counts
mt.newdf = data.frame(TextBarcodes = apply(mt.df, 1, function(myrow) paste(myrow[2:ncol(mt.df)], collapse="")), GeneCount = rowSums( mt.df[,2:ncol(mt.df)] ) )

# find and sort unique barcodes in this set
mt.uniques = unique(mt.newdf)
mt.uniques = mt.uniques[order(mt.uniques$GeneCount),]

# add species names to intermediate dataframe
mt.newdf$Species = mt.species.list

mt.org = read.csv(mtgenestatsfilename, header=T, stringsAsFactor=F)
mt.indices = read.csv(mtindexfilename, header=T, stringsAsFactor=F)

mt.org.df = mt.org[mt.org$Protocol=="All",]

#### get phylogenetic information

# read taxonomy tree
mt.treeString = tolower(paste(readLines(mttreefilename), collapse=""))
mt.treeString = gsub(" ", "_", mt.treeString)
mt.tree <- read.tree(text=mt.treeString)
mt.tree$tip.label = gsub("_", " ", mt.tree$tip.label)
mt.to.drop = setdiff(mt.tree$tip.label, mt.species.list)
mt.tree = drop.tip(mt.tree, mt.to.drop)

# list of all node labels
mt.all.labels = c(mt.tree$tip.label, mt.tree$node.label)

# numbers of leaves and internal nodes
mt.n.tips = length(mt.tree$tip.label)
mt.n.nodes = length(mt.tree$node.label)

# identify tree root
mt.luca.node = grep("eukaryota", mt.all.labels)

# identify those clades that are directly descended from root...
mt.clades.refs = mt.tree$edge[which(mt.tree$edge[,1]==mt.luca.node),2]

# ... and which are not leaves
mt.clades.refs = mt.clades.refs[mt.clades.refs > mt.n.tips]

# get and sort labels for these clades
mt.clades = mt.all.labels[mt.clades.refs]
mt.clades = mt.clades[order(mt.clades)]

message("Building phylogenetic info...")

# build list of lists of organisms in each clade
mt.clades.list = list()
for(i in 1:length(mt.clades)) {
  # identify root of clade
  mt.head.node = which(mt.all.labels == mt.clades[i])
  # get all descendants
  mt.head.des = getDescendants(mt.tree, mt.head.node)
  # get labels of descendants that are leaves
  mt.species = gsub("_", " ", mt.all.labels[mt.head.des[mt.head.des <= mt.n.tips]])
  # append these labels to this clade's list
  mt.clades.list = append(mt.clades.list, list(mt.species))
}

# count clade occupancies for each barcode (number of species in the clade with that barcode)
# initialise dataframe
mt.clade.occupancy = data.frame(TextBarcodes = mt.uniques$TextBarcodes)
# run through all clades
for(i in 1:length(mt.clades)) {
  # initialise new column in dataframe for this clade
  mt.clade.occupancy[,mt.clades[i]] = NA
  # run through unique barcodes
  for(j in 1:nrow(mt.uniques)) {
    # find those indices in the full species records that match this barcodes
    mt.refs = which(mt.newdf$TextBarcode == mt.uniques$TextBarcodes[j])
    # get their corresponding species names
    mt.species = trimws(mt.newdf$Species[mt.refs])
    # record the size of the intersection between this barcode's (j) species set and this clade's (i) species set
    mt.clade.occupancy[j,mt.clades[i]] = length(intersect(mt.species, mt.clades.list[[i]]))
  }
}

#### styling and element creation for plotting

label.size = 8
clade.label.size = 5

message("Styling barcode plot...")

# this section identifies groups of gene labels matching patterns, so that we can illustrate genes belonging to the same complex
mt.gene.names <- colnames(mt.df)[2:ncol(mt.df)]
if(length(grep("mt", mtbarcodefilename)) > 0) {
  mt.gois = c("atp", "cox", "cytb", "nad", "rpl", "rps", "sdh")
} else {
  mt.gois = c("acc", "atp", "ndh", "pet", "psa", "psb", "rbc", "rpl", "rps")
}
# initialise list of start and end references
mt.goi.starts = mt.goi.ends = NULL
# loop through patterns
for(i in 1:length(mt.gois)) {
  # record first and last gene label index that match this pattern
  mt.refs = grep(mt.gois[i], mt.gene.names)
  mt.goi.starts = c(mt.goi.starts, min(mt.refs))
  mt.goi.ends = c(mt.goi.ends, max(mt.refs))
}

mt.gene.cols.index = NULL
for(i in 1:length(mt.gene.names)) {
  mt.ref = which(mt.indices$GeneLabel == mt.gene.names[i])
  mt.indexcol = round(100* (mt.indices$Index[mt.ref]))
#  mt.gene.cols.index = c(mt.gene.cols.index, rainbow(250, v=0.75)[mt.indexcol+150])
  mt.gene.cols.index = c(mt.gene.cols.index, rgb(1-0.01*mt.indexcol, 0.5-0.005*mt.indexcol, 0.5-0.005*mt.indexcol))
}

# create lists of x and y coordinates of pixels corresponding to "1"s in this unique set
mt.pix.xs = mt.pix.ys = mt.pix.cols.index = NULL
mt.bar.xs = mt.bar.ys = mt.bar.ws = NULL

# run through unique barcodes
for(y in 1:nrow(mt.uniques))
{
  # run through set of genes of interest, adding pixel coordinates if presence is recorded in this barcode
  for(x in 1:mt.ngenes)
  {
    if(substr(mt.uniques$TextBarcode[y], x, x) == "1")
    {
      mt.pix.xs = c(mt.pix.xs, x)
      mt.pix.ys = c(mt.pix.ys, y)
      mt.pix.cols.index = c(mt.pix.cols.index, mt.gene.cols.index[x])
    }
  }
  # run through list of clades, adding bar coordinates and widths corresponding to occupancy for this barcode
  for(x in 1:length(mt.clades.list))
  {
    if(mt.clade.occupancy[y,mt.clades[x]] >= 1)
    {
      mt.bar.xs = c(mt.bar.xs, x)
      mt.bar.ys = c(mt.bar.ys, y)
      mt.bar.ws = c(mt.bar.ws, log10(mt.clade.occupancy[y,mt.clades[x]]))
    }
  }
}

# various specific styling parameters
mt.clade.x.scale = 1.5
mt.clade.x.offset = -mt.clade.x.scale*length(mt.clades.list)-2
mt.clade.y.scale = 25
mt.gene.y.scale = 35
mt.bar.min.size = 0.5
mt.complex.y.offset = -60
mt.complex.y.scale = 8
mt.complex.bar.width = 4

# produce vertical offsets for complex, clade, and gene labels
# these are staggered for legibility
mt.complex.y.offsets = rep(mt.complex.y.offset+mt.complex.y.scale*c(0,-1), length(mt.gois))[seq(1,length(mt.gois))]
mt.clade.y.offsets = rep(mt.clade.y.scale*c(0,-1), length(mt.clades.list))[seq(1,length(mt.clades.list))]
mt.gene.y.offsets = rep(mt.gene.y.scale*c(0,-0.5,-1), ncol(mt.df)-1)[seq(1,ncol(mt.df)-1)]

# produce dataframes of rectangles, for pixels, clade occupancy, clade background, and complex label bars
mt.pixels = data.frame(x1 = mt.pix.xs, x2 = mt.pix.xs+1, y1 = mt.pix.ys, y2 = mt.pix.ys+1, pcol.index = mt.pix.cols.index, stringsAsFactors=F)
mt.bars = data.frame(x1 = mt.clade.x.scale*(mt.bar.xs+0.5)-mt.bar.ws-mt.bar.min.size+mt.clade.x.offset, x2 = mt.clade.x.scale*(mt.bar.xs+0.5)+mt.bar.ws+mt.bar.min.size+mt.clade.x.offset, y1 = mt.bar.ys, y2 = mt.bar.ys+1)
mt.cladebars = data.frame(x1 = mt.clade.x.scale*seq(1, length(mt.clades.list), by = 1)+mt.clade.x.offset, y1 = rep(0, length(mt.clades.list)), x2 = mt.clade.x.scale*(seq(1, length(mt.clades.list), by = 1)+1)+mt.clade.x.offset, y2 = rep(nrow(mt.uniques), length(mt.clades.list)))
mt.complexbars = data.frame(x1 = mt.goi.starts, y1 = mt.complex.y.offsets+1, x2 = mt.goi.ends, y2 = mt.complex.y.offsets+mt.complex.bar.width)

# data frames of text labels for genes, clades, and complexes
mt.genelabel = data.frame(x = seq(1, ncol(mt.df)-1, by=1), y = mt.gene.y.offsets, label=colnames(mt.df)[2:ncol(mt.df)])
mt.cladelabel = data.frame(x = mt.clade.x.scale*(seq(1, length(mt.clades.list), by = 1)+0.5)+mt.clade.x.offset, y = mt.clade.y.offsets, label=substr(mt.clades, 1, 5))
mt.complexlabel = data.frame(x = mt.goi.starts, y = mt.complex.y.offsets, label = mt.gois)

# colour scheme for clade labels
mt.custom.col = NULL
for(i in 1:length(mt.clades))
{
  if(i %% 2 == 0)
  {
    mt.custom.col = c(mt.custom.col, rgb(0.4, 0.4, 0.4))
  }
  else
  {
    mt.custom.col = c(mt.custom.col, rgb(0.7, 0.7, 0.7))
  }
}

if(length(grep("mt", mtbarcodefilename)) > 0) {
  plot.width = 1000
  plot.height = 750
} else {
  plot.width = 2000
  plot.height = 1500
}
  
mt.ggplot = ggplot() +
  geom_rect(data=mt.cladebars, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=mt.clades), alpha=0.2) +       # background bars for clade section
  geom_rect(data=mt.pixels, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=mt.pixels$pcol.index, color=NA) +           # plot presence/absence pixels
  geom_rect(data=mt.bars, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2)) +                                    # plot clade occupancy bars
  geom_text(data=mt.cladelabel, mapping=aes(x=x,y=y,label=label,color=mt.clades), angle=90, hjust="right", size=clade.label.size) +       # clade labels
  geom_text(data=mt.genelabel, mapping=aes(x=x,y=y,label=label), angle=90, hjust="right") +                     # gene labels
    geom_rect(data=mt.complexbars, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="gray", color="gray") +  # complex bars
  geom_text(data=mt.complexlabel, mapping=aes(x=x,y=y,label=label), hjust="left", size=label.size) +                             # complex labelse
  theme_void() + theme(legend.position = "none") + 
  scale_color_manual(values=mt.custom.col) + scale_fill_manual(values=mt.custom.col)

#################### PT

#### get barcode data

message("Getting data...")

# read barcode data
pt.df <- read.csv(ptbarcodefilename, header=T)
pt.ngenes = ncol(pt.df)-1

# get species list
pt.species.list = trimws(pt.df$Species)
# create new data frame, storing text versions of barcodes and total gene counts
pt.newdf = data.frame(TextBarcodes = apply(pt.df, 1, function(myrow) paste(myrow[2:ncol(pt.df)], collapse="")), GeneCount = rowSums( pt.df[,2:ncol(pt.df)] ) )

# find and sort unique barcodes in this set
pt.uniques = unique(pt.newdf)
pt.uniques = pt.uniques[order(pt.uniques$GeneCount),]

# add species names to intermediate dataframe
pt.newdf$Species = pt.species.list

pt.org = read.csv(ptgenestatsfilename, header=T, stringsAsFactor=F)
pt.indices = read.csv(ptindexfilename, header=T, stringsAsFactor=F)

pt.org.df = pt.org[pt.org$Protocol=="All",]

#### get phylogenetic information

# read taxonomy tree
pt.treeString = tolower(paste(readLines(pttreefilename), collapse=""))
pt.treeString = gsub(" ", "_", pt.treeString)
pt.tree <- read.tree(text=pt.treeString)
pt.tree$tip.label = gsub("_", " ", pt.tree$tip.label)
pt.to.drop = setdiff(pt.tree$tip.label, pt.species.list)
pt.tree = drop.tip(pt.tree, pt.to.drop)

# list of all node labels
pt.all.labels = c(pt.tree$tip.label, pt.tree$node.label)

# numbers of leaves and internal nodes
pt.n.tips = length(pt.tree$tip.label)
pt.n.nodes = length(pt.tree$node.label)

# identify tree root
pt.luca.node = grep("eukaryota", pt.all.labels)

# identify those clades that are directly descended from root...
pt.clades.refs = pt.tree$edge[which(pt.tree$edge[,1]==pt.luca.node),2]

# ... and which are not leaves
pt.clades.refs = pt.clades.refs[pt.clades.refs > pt.n.tips]

# get and sort labels for these clades
pt.clades = pt.all.labels[pt.clades.refs]
pt.clades = pt.clades[order(pt.clades)]

message("Building phylogenetic info...")

# build list of lists of organisms in each clade
pt.clades.list = list()
for(i in 1:length(pt.clades)) {
  # identify root of clade
  pt.head.node = which(pt.all.labels == pt.clades[i])
  # get all descendants
  pt.head.des = getDescendants(pt.tree, pt.head.node)
  # get labels of descendants that are leaves
  pt.species = gsub("_", " ", pt.all.labels[pt.head.des[pt.head.des <= pt.n.tips]])
  # append these labels to this clade's list
  pt.clades.list = append(pt.clades.list, list(pt.species))
}

# count clade occupancies for each barcode (number of species in the clade with that barcode)
# initialise dataframe
pt.clade.occupancy = data.frame(TextBarcodes = pt.uniques$TextBarcodes)
# run through all clades
for(i in 1:length(pt.clades)) {
  # initialise new column in dataframe for this clade
  pt.clade.occupancy[,pt.clades[i]] = NA
  # run through unique barcodes
  for(j in 1:nrow(pt.uniques)) {
    # find those indices in the full species records that match this barcodes
    pt.refs = which(pt.newdf$TextBarcode == pt.uniques$TextBarcodes[j])
    # get their corresponding species names
    pt.species = trimws(pt.newdf$Species[pt.refs])
    # record the size of the intersection between this barcode's (j) species set and this clade's (i) species set
    pt.clade.occupancy[j,pt.clades[i]] = length(intersect(pt.species, pt.clades.list[[i]]))
  }
}

#### styling and element creation for plotting

message("Styling barcode plot...")

# this section identifies groups of gene labels matching patterns, so that we can illustrate genes belonging to the same complex
pt.gene.names <- colnames(pt.df)[2:ncol(pt.df)]
if(length(grep("mt", ptbarcodefilename)) > 0) {
  pt.gois = c("atp", "cox", "cytb", "nad", "rpl", "rps", "sdh")
} else {
  pt.gois = c("acc", "atp", "ndh", "pet", "psa", "psb", "rbc", "rpl", "rps")
}
# initialise list of start and end references
pt.goi.starts = pt.goi.ends = NULL
# loop through patterns
for(i in 1:length(pt.gois)) {
  # record first and last gene label index that match this pattern
  pt.refs = grep(pt.gois[i], pt.gene.names)
  pt.goi.starts = c(pt.goi.starts, min(pt.refs))
  pt.goi.ends = c(pt.goi.ends, max(pt.refs))
}

pt.gene.cols.index = NULL
for(i in 1:length(pt.gene.names)) {
  pt.ref = which(pt.indices$GeneLabel == pt.gene.names[i])
  pt.indexcol = round(100* (pt.indices$Index[pt.ref]))
  pt.gene.cols.index = c(pt.gene.cols.index, rgb(0.5 - 0.005*pt.indexcol, 0.5 - 0.005*pt.indexcol, 1-0.01*pt.indexcol))
}

# create lists of x and y coordinates of pixels corresponding to "1"s in this unique set
pt.pix.xs = pt.pix.ys = pt.pix.cols.index = NULL
pt.bar.xs = pt.bar.ys = pt.bar.ws = NULL

# run through unique barcodes
for(y in 1:nrow(pt.uniques))
{
  # run through set of genes of interest, adding pixel coordinates if presence is recorded in this barcode
  for(x in 1:pt.ngenes)
  {
    if(substr(pt.uniques$TextBarcode[y], x, x) == "1")
    {
      pt.pix.xs = c(pt.pix.xs, x)
      pt.pix.ys = c(pt.pix.ys, y)
      pt.pix.cols.index = c(pt.pix.cols.index, pt.gene.cols.index[x])
    }
  }
  # run through list of clades, adding bar coordinates and widths corresponding to occupancy for this barcode
  for(x in 1:length(pt.clades.list))
  {
    if(pt.clade.occupancy[y,pt.clades[x]] >= 1)
    {
      pt.bar.xs = c(pt.bar.xs, x)
      pt.bar.ys = c(pt.bar.ys, y)
      pt.bar.ws = c(pt.bar.ws, log10(pt.clade.occupancy[y,pt.clades[x]]))
    }
  }
}

# various specific styling parameters
pt.clade.x.scale = 1
pt.clade.x.offset = -length(pt.clades.list)*pt.clade.x.scale-2
pt.clade.y.scale = 60
pt.gene.y.scale = 75
pt.bar.min.size = 0.5
pt.complex.y.offset = -140
pt.complex.y.scale=10
pt.complex.bar.width = 6

# produce vertical offsets for complex, clade, and gene labels
# these are staggered for legibility
pt.complex.y.offsets = rep(pt.complex.y.offset+pt.complex.y.scale*c(0,-1), length(pt.gois))[seq(1,length(pt.gois))]
pt.clade.y.offsets = rep(pt.clade.y.scale*c(0,-1), length(pt.clades.list))[seq(1,length(pt.clades.list))]
pt.gene.y.offsets = rep(pt.gene.y.scale*c(0,-0.5,-1), ncol(pt.df)-1)[seq(1,ncol(pt.df)-1)]

# produce dataframes of rectangles, for pixels, clade occupancy, clade background, and complex label bars
pt.pixels = data.frame(x1 = pt.pix.xs, x2 = pt.pix.xs+1, y1 = pt.pix.ys, y2 = pt.pix.ys+1, pcol.index = pt.pix.cols.index, stringsAsFactors=F)
pt.bars = data.frame(x1 = pt.clade.x.scale*(pt.bar.xs+0.5)-pt.bar.ws-pt.bar.min.size+pt.clade.x.offset, x2 = pt.clade.x.scale*(pt.bar.xs+0.5)+pt.bar.ws+pt.bar.min.size+pt.clade.x.offset, y1 = pt.bar.ys, y2 = pt.bar.ys+1)
pt.cladebars = data.frame(x1 = pt.clade.x.scale*seq(1, length(pt.clades.list), by = 1)+pt.clade.x.offset, y1 = rep(0, length(pt.clades.list)), x2 = pt.clade.x.scale*seq(1, length(pt.clades.list), by = 1)+pt.clade.x.offset+1, y2 = rep(nrow(pt.uniques), length(pt.clades.list)))
pt.complexbars = data.frame(x1 = pt.goi.starts, y1 = pt.complex.y.offsets+1, x2 = pt.goi.ends, y2 = pt.complex.y.offsets+pt.complex.bar.width)

# data frames of text labels for genes, clades, and complexes
pt.genelabel = data.frame(x = seq(1, ncol(pt.df)-1, by=1), y = pt.gene.y.offsets, label=colnames(pt.df)[2:ncol(pt.df)])
pt.cladelabel = data.frame(x = pt.clade.x.scale*(seq(1, length(pt.clades.list), by = 1)+0.5)+pt.clade.x.offset, y = pt.clade.y.offsets, label=substr(pt.clades, 1, 5))
pt.complexlabel = data.frame(x = pt.goi.starts, y = pt.complex.y.offsets, label = pt.gois)

# colour scheme for clade labels
pt.custom.col = NULL
for(i in 1:length(pt.clades))
{
  if(i %% 2 == 0)
  {
    pt.custom.col = c(pt.custom.col, rgb(0.4, 0.4, 0.4))
  }
  else
  {
    pt.custom.col = c(pt.custom.col, rgb(0.7, 0.7, 0.7))
  }
}

if(length(grep("mt", ptbarcodefilename)) > 0) {
  plot.width = 1000
  plot.height = 750
} else {
  plot.width = 2000
  plot.height = 1500
}
  
pt.ggplot = ggplot() +
  geom_rect(data=pt.cladebars, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=pt.clades), alpha=0.2) +       # background bars for clade section
  geom_rect(data=pt.pixels, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill=pt.pixels$pcol.index, color=NA) +           # plot presence/absence pixels
  geom_rect(data=pt.bars, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2)) +                                    # plot clade occupancy bars
  geom_text(data=pt.cladelabel, mapping=aes(x=x,y=y,label=label,color=pt.clades), angle=90, hjust="right", size=clade.label.size) +       # clade labels
  geom_text(data=pt.genelabel, mapping=aes(x=x,y=y,label=label), angle=90, hjust="right") +                     # gene labels
    geom_rect(data=pt.complexbars, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="gray", color="gray") +  # complex bars
  geom_text(data=pt.complexlabel, mapping=aes(x=x,y=y,label=label), hjust="left", size=label.size) +                             # complex labelse
  theme_void() + theme(legend.position = "none") + 
  scale_color_manual(values=pt.custom.col) + scale_fill_manual(values=pt.custom.col)


  outfile = paste(c(out.dir, "both-barcodes.png"), collapse="")
res.factor = 3
png(outfile, width=1600*res.factor,height=800*res.factor, res=72*res.factor)
plot_grid(mt.ggplot, pt.ggplot, rel_widths=c(1/4, 3/4))
dev.off()

###########################

##############

#  taxa.names = tolower(c("Rhodophyta", "Oomycetes", "Apicomplexa", "Gelidium", "Viridiplantae", "Metazoa", "Eukaryota", "Fungi", "Magnoliopsida", "Streptophyta", "Chlorophyta", "Bryopsida", "Chlorophyceae", "Ascomycota", "Basidiomycota", "Sordariomycetes", "Saccharomycetes", "Lecanoromycetes", "Saccharomycetaceae", "Parmeliaceae", "Chordata", "Arthropoda", "Insecta", "Actinopteri", "Mammalia", "Aves", "Mollusca", "Nematoda", "actinopterygii", "sarcopterygii", "spiralia", "hexapoda", "jakobida"))
#  taxa.names = tolower(c("Rhodophyta", "Oomycetes", "Apicomplexa", "Gelidium", "Viridiplantae", "Metazoa", "Eukaryota", "Fungi", "Magnoliopsida", "Streptophyta", "Chlorophyta", "Bryopsida", "Chlorophyceae", "Ascomycota", "Basidiomycota", "Sordariomycetes", "Saccharomycetes", "Lecanoromycetes", "Saccharomycetaceae", "Parmeliaceae", "Chordata", "Arthropoda", "Insecta", "Actinopteri", "Mammalia", "Aves", "Mollusca", "Nematoda", "Poales", "Asparagales", "Asterales", "Solanales", "Rosales", "Caryophyllales", "Polypodiopsida", "Pinopsida", "Phaeophyceae", "Bacillariophyta", "Euglenozoa", "Ericales", "Poaceae", "Lamiales", "Fabales"))

message("Producing tree...")

# produce igraph object from taxonomy tree
pt.tree$tip.label = gsub("_", " ", pt.tree$tip.label)
pt.counts = data.frame(Label = NULL, Count = NULL)
for(label in pt.tree$tip.label) {
  ref = which(pt.df$Species == label)
  if(length(ref) != 0) {
    pt.counts = rbind(pt.counts, data.frame(Label=label, Count = sum(pt.df[ref, 2:ncol(pt.df)])))
  } else {
    print("wtf")
  }
}

# produce igraph object from taxonomy tree
mt.tree$tip.label = gsub("_", " ", mt.tree$tip.label)
mt.counts = data.frame(Label = NULL, Count = NULL)
for(label in mt.tree$tip.label) {
  ref = which(mt.df$Species == label)
  if(length(ref) != 0) {
    mt.counts = rbind(mt.counts, data.frame(Label=label, Count = sum(mt.df[ref, 2:ncol(mt.df)])))
  } else {
    print("wtf")
  }
}

clade.text.size = 10

for(version in 1:2) {

  if(version == 1) {
    layout.type = "circular"
    outfile = paste(c(out.dir, "trees-circ.png"), collapse="")
    horiz.setting = FALSE
    angle.setting = "auto"
  } else {
    layout.type = "slanted"
    outfile = paste(c(out.dir, "trees-slant.png"), collapse="")
    horiz.setting = TRUE
    angle.setting = 0
  }
  
  pt.clades.new = c("Rhodophyta", "Bacillariophyta", "Apicomplexa", "Viridiplantae", "Phaeophyceae")
  pt.clade.refs = data.frame(node=which(pt.all.labels %in% tolower(pt.clades.new)))
  pt.clade.refs$clade.label = pt.all.labels[pt.clade.refs$node]
  
  pt.collapse = which(pt.all.labels == "viridiplantae")
  
  #pt.1.labset = rep(FALSE, length(pt.tree$tip.label))
  #pt.1.labset[sample(length(pt.tree$tip.label), 0.75*length(pt.tree$tip.label))] = TRUE
  pt.1 = ggtree(pt.tree, layout=layout.type, size=0.1, branch.length="none") 
  pt.1.label = collapse(pt.1, node=pt.collapse) +
  #  geom_tiplab(size=3, offset = 1.5, color="#AAAAAA", aes(subset=pt.1.labset)) +
    geom_fruit(data = pt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="blue", breaks=c(0, max(pt.counts$Count))) +
    geom_cladelab(data=pt.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    geom_point2(aes(subset=(node == pt.collapse)), size=2, shape=23, fill="steelblue") +
    theme(legend.position="none")
  
  pt.2.tree = extract.clade(pt.tree, pt.collapse)
  pt.2.all.labels = c(pt.2.tree$tip.label, pt.2.tree$node.label)
  pt.2.collapse = which(pt.2.all.labels == "magnoliopsida")
  
  pt.2.labset = rep(FALSE, length(pt.2.tree$tip.label))
  pt.2.labset[sample(length(pt.2.tree$tip.label), 0.25*length(pt.2.tree$tip.label))] = TRUE
  pt.2.clades.new = c("Chlorophyta", "Streptophyta", "Magnoliopsida")
  pt.2.clade.refs = data.frame(node=which(pt.2.all.labels %in% tolower(pt.2.clades.new)))
  pt.2.clade.refs$clade.label = pt.2.all.labels[pt.2.clade.refs$node]
  pt.2 = ggtree(pt.2.tree, layout=layout.type, size=0.1, branch.length="none")
  pt.2.label = collapse(pt.2, node=pt.2.collapse) +
  #  geom_tiplab(size=3, offset = 1.5, color="#AAAAAA", aes(subset=pt.2.labset)) +
    geom_fruit(data = pt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="blue", breaks=c(0, max(pt.counts$Count))) +
    geom_cladelab(data=pt.2.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    geom_point2(aes(subset=(node == pt.2.collapse)), size=2, shape=23, fill="steelblue") +
    theme(legend.position="none")
  
  pt.3.tree = extract.clade(pt.2.tree, pt.2.collapse)
  pt.3.all.labels = c(pt.3.tree$tip.label, pt.3.tree$node.label)
  
  pt.3.labset = rep(FALSE, length(pt.3.tree$tip.label))
  pt.3.labset[sample(length(pt.3.tree$tip.label), 0.05*length(pt.3.tree$tip.label))] = TRUE
  pt.3.clades.new = c("Asparagales", "Poales", "Asterales", "Lamiales", "Solanales", "Fabales", "Ericales", "Rosales")
  pt.3.clade.refs = data.frame(node=which(pt.3.all.labels %in% tolower(pt.3.clades.new)))
  pt.3.clade.refs$clade.label = pt.3.all.labels[pt.3.clade.refs$node]
  pt.3 = ggtree(pt.3.tree, layout=layout.type, size=0.1, branch.length="none")
  pt.3.label = pt.3 +
  #  geom_tiplab(size=4, offset = 1.5, color="#AAAAAA", aes(subset=pt.3.labset)) +
    geom_fruit(data = pt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="blue", breaks=c(0, max(pt.counts$Count))) +
    geom_cladelab(data=pt.3.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    theme(legend.position="none")
  
  
  
  ###################################
  
  
  mt.clades.new = c("Rhodophyta", "Oomycetes", "Apicomplexa", "Viridiplantae", "Metazoa", "Fungi", "Jakobida")
  mt.clade.refs = data.frame(node=which(mt.all.labels %in% tolower(mt.clades.new)))
  mt.clade.refs$clade.label = mt.all.labels[mt.clade.refs$node]
  
  mt.collapse = which(mt.all.labels == "metazoa")
  mt.1 = ggtree(mt.tree, layout=layout.type, size=0.1, branch.length="none")
  
  #mt.1.labset = rep(FALSE, length(mt.tree$tip.label))
  #mt.1.labset[sample(length(mt.tree$tip.label), 0.2*length(mt.tree$tip.label))] = TRUE
  mt.1.label = collapse(mt.1, node=mt.collapse) +
  #  geom_tiplab(size=3, offset = 1.5, color="#AAAAAA", aes(subset=mt.1.labset)) +
    geom_fruit(data = mt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="red", breaks=c(0, max(mt.counts$Count))) +
    geom_cladelab(data=mt.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    geom_point2(aes(subset=(node == mt.collapse)), size=2, shape=23, fill="steelblue") +
    theme(legend.position="none")
  
  mt.2.tree = extract.clade(mt.tree, mt.collapse)
  mt.2.all.labels = c(mt.2.tree$tip.label, mt.2.tree$node.label)
  mt.2.collapse = which(mt.2.all.labels == "chordata")
  
  mt.2.labset = rep(FALSE, length(mt.2.tree$tip.label))
  mt.2.labset[sample(length(mt.2.tree$tip.label), 0.05*length(mt.2.tree$tip.label))] = TRUE
  mt.2.clades.new = c("Chordata", "Nematoda", "Mollusca", "Arthropoda")
  mt.2.clade.refs = data.frame(node=which(mt.2.all.labels %in% tolower(mt.2.clades.new)))
  mt.2.clade.refs$clade.label = mt.2.all.labels[mt.2.clade.refs$node]
  mt.2 = ggtree(mt.2.tree, layout=layout.type, size=0.1, branch.length="none")
  mt.2.label = collapse(mt.2, node=mt.2.collapse) +
  #  geom_tiplab(size=3, offset = 1.5, color="#AAAAAA", aes(subset=mt.2.labset)) +
    geom_fruit(data = mt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="red", breaks=c(0, max(mt.counts$Count))) +
    geom_cladelab(data=mt.2.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    geom_point2(aes(subset=(node == mt.2.collapse)), size=2, shape=23, fill="steelblue") +
    theme(legend.position="none")
  
  mt.3.tree = extract.clade(mt.2.tree, mt.2.collapse)
  mt.3.all.labels = c(mt.3.tree$tip.label, mt.3.tree$node.label)
  mt.3.collapse = which(mt.3.all.labels == "mammalia")
  
  mt.3.labset = rep(FALSE, length(mt.3.tree$tip.label))
  mt.3.labset[sample(length(mt.3.tree$tip.label), 0.1*length(mt.3.tree$tip.label))] = TRUE
  mt.3.clades.new = c("Mammalia", "Aves", "Actinopteri")
  mt.3.clade.refs = data.frame(node=which(mt.3.all.labels %in% tolower(mt.3.clades.new)))
  mt.3.clade.refs$clade.label = mt.3.all.labels[mt.3.clade.refs$node]
  mt.3 = ggtree(mt.3.tree, layout=layout.type, size=0.1, branch.length="none")
  mt.3.label = collapse(mt.3, node=mt.3.collapse) + 
  #  geom_tiplab(size=3, offset = 1.5, color="#AAAAAA", aes(subset=mt.3.labset)) +
    geom_fruit(data = mt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="red", breaks=c(0, max(mt.counts$Count))) +
    geom_cladelab(data=mt.3.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    geom_point2(aes(subset=(node == mt.3.collapse)), size=2, shape=23, fill="steelblue") +
    theme(legend.position="none")
  
  mt.4.tree = extract.clade(mt.3.tree, mt.3.collapse)
  mt.4.all.labels = c(mt.4.tree$tip.label, mt.4.tree$node.label)
  
  mt.4.labset = rep(FALSE, length(mt.4.tree$tip.label))
  mt.4.labset[sample(length(mt.4.tree$tip.label), 0.05*length(mt.4.tree$tip.label))] = TRUE
  mt.4.clades.new = c("Rodentia", "Chiroptera", "Primates")
  mt.4.clade.refs = data.frame(node=which(mt.4.all.labels %in% tolower(mt.4.clades.new)))
  mt.4.clade.refs$clade.label = mt.4.all.labels[mt.4.clade.refs$node]
  mt.4 = ggtree(mt.4.tree, layout=layout.type, size=0.1, branch.length="none")
  mt.4.label = mt.4 +
  #  geom_tiplab(size=3, offset = 1.5, color="#AAAAAA", aes(subset=mt.4.labset)) +
    geom_fruit(data = mt.counts, geom=geom_bar, mapping = aes(y=Label, x=Count, fill=Count), orientation="y", stat="identity") +
    scale_fill_gradient(low="black", high="red", breaks=c(0, max(mt.counts$Count))) +
    geom_cladelab(data=mt.4.clade.refs, mapping=aes(node=node, label=clade.label), angle=angle.setting, horizontal=horiz.setting, offset = 2, offset.text = 0.5, hjust=0.5, fontsize = clade.text.size) +
    theme(legend.position="none")
  
  
  res.factor=1
  png(outfile, width=3000*res.factor, height=2000*res.factor, res=72*res.factor)
  grid.arrange(mt.1.label, mt.2.label, mt.3.label, mt.4.label, pt.1.label, pt.2.label, pt.3.label, nrow=2)
  dev.off()
}

