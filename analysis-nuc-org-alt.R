#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 8) {
  stop("Need: MT stats file, PT stats file, NU stats file, output label, output path, plot path, test proportion, number of samples")
}

#args = c("Data/mt-stats-manual.csv", "Data/pt-stats-manual.csv", "Data/all-stats.csv", "manual", "Data/", "Plots/", "0.5", "10")

message("Loading libraries...")

library(randomForest)
library(ggplot2)
library(hexbin)
library(ggnewscale)
library(gridExtra)
library(tree)
set.seed(1)

source("lengthNormalise.R")

message("Reading data...")

# read pipeline outputs
mt.df = read.csv(args[1], header=T, stringsAsFactor=T)
pt.df = read.csv(args[2], header=T, stringsAsFactor=T)
nuc.df = read.csv(args[3], header=T, stringsAsFactor=T)

treereducedplotoutput = paste(c(args[6], "nuc-org-", args[4], "-tree-reduced.png"), collapse="")
treecrossplotoutput = paste(c(args[6], "nuc-org-", args[4], "-tree-cross.png"), collapse="")
treeplotoutput = paste(c(args[6], "nuc-org-", args[4], "-tree.png"), collapse="")
rfreducedplotoutput = paste(c(args[6], "nuc-org-", args[4], "-rf-reduced.png"), collapse="")
rfplotoutput = paste(c(args[6], "nuc-org-", args[4], "-rf.png"), collapse="")
rfcrossplotoutput = paste(c(args[6], "nuc-org-", args[4], "-rf-cross.png"), collapse="")
hexesoutput.gc = paste(c(args[6], "nuc-org-", args[4], "-hydro-gc-hexes.png"), collapse="")
hexesoutput.pka1 = paste(c(args[6], "nuc-org-", args[4], "-hydro-pka1-hexes.png"), collapse="")
hexesoutput.pka2 = paste(c(args[6], "nuc-org-", args[4], "-hydro-pka2-hexes.png"), collapse="")
hexesoutput.cw = paste(c(args[6], "nuc-org-", args[4], "-hydro-cw-hexes.png"), collapse="")
hexesoutput.limited = paste(c(args[6], "nuc-org-", args[4], "-hydro-pka1-hexes-limited.png"), collapse="")
resultsoutput = paste(c(args[5], "nuc-org-", args[4], "-results.csv"), collapse="")

testpropn = as.numeric(args[7])
nsamp = as.numeric(args[8])

nuc.df = nuc.df[nuc.df$Compartment=="NU",]

# convert species to lower case
nuc.df$Species = tolower(nuc.df$Species)

# regexs to recognise genes encoding subunits of complexes of interest
mt.regex = c("nad[0-9]", "sdh[0-9]", "cytb", "cox[0-9]", "atp[0-9]", "rp[sl]") 
pt.regex = c("psa[a-x]", "psb[a-z]", "atp[a-z]", "pet[a-z]", "rbc", "rp[sl]") 
nuc.mt.regex = c("mt[-]c1", "mt[-]c2", "mt[-]c3", "mt[-]c4", "mt[-]c5", "mt[-]ri") 
nuc.pt.regex = c("pt[-]pa", "pt[-]pb", "mt[-]c5", "pt[-]cb", "pt[-]rb", "pt[-]ri") 

# create list of dataframes for nuc-org comparison
# one list element for each complex of interest
# only retain entries for which the species appears both in organelle records and nuclear records
mt.compare = list()
for(i in 1:length(mt.regex)) {
  mt.org = mt.df[grepl(mt.regex[i], mt.df$GeneLabel),] 
  mt.nuc = nuc.df[grepl(nuc.mt.regex[i], nuc.df$GeneLabel),]
  mt.compare[[length(mt.compare)+1]] = rbind(mt.org, mt.nuc)
#  intersection = intersect(mt.org$Species, mt.nuc$Species)
#  mt.compare[[length(mt.compare)+1]] = rbind(mt.org[mt.org$Species %in% intersection,], mt.nuc[mt.nuc$Species %in% intersection,])
}

pt.compare = list()
for(i in 1:length(pt.regex)) {
  pt.org = pt.df[grepl(pt.regex[i], pt.df$GeneLabel),]
  pt.nuc = nuc.df[grepl(nuc.pt.regex[i], nuc.df$GeneLabel),]
  pt.compare[[length(pt.compare)+1]] = rbind(pt.org, pt.nuc)
#  intersection = intersect(pt.org$Species, pt.nuc$Species)
#  pt.compare[[length(pt.compare)+1]] = rbind(pt.org[pt.org$Species %in% intersection,], pt.nuc[pt.nuc$Species %in% intersection,])
}

# normalise selected statistics by length
mt.compare = lapply(mt.compare, lengthNormalise)
pt.compare = lapply(pt.compare, lengthNormalise)

#### hex plots

message("Hex plots...")

# labels for plots
plot.labels = c(" MT Complex I", " MT Complex II", " MT Complex III", " MT Complex IV", " MT Complex V", " MT Ribosome", " PT Photosystem I", " PT Photosystem II", " PT ATP Synthase", " PT Cytochrome b6f", " PT Rubisco", " PT Ribosome")

################## create list of plots, one for each complex of interest
####### first: hydro vs gc

my.plots = list()
for(subdf in mt.compare) {
  # two dataframes, one for nuclear-encoded and one for organelle-encoded genes
  df.nu = subdf[subdf$Compartment == "NU",]
  df.mt = subdf[subdf$Compartment == "MT",]
  ref = length(my.plots)+1
  
  my.plots[[ref]] = ggplot(mapping = aes(x=GC, y=Hydro)) +                                      # plot GC vs hydrophobicity
    geom_hex(data = df.nu, alpha=0.5, col=NA) +                                                 # use geom_hex to reflect point density for nuc genes
    scale_fill_gradientn(colours=c("#DDDDDD","#333333"),name='Frequency',na.value=NA) +		# blue scale for "frequency" (point density)
    new_scale("fill") +
    geom_hex(data = df.mt, alpha = 0.5, col=NA) +						# same for org genes
    scale_fill_gradientn(colours=c("#FFCCCC","#FF0000")) +                                      # red scale
    ggtitle(plot.labels[ref]) +                                                                 # styling elements follow beneath
    theme_classic() +
    theme(axis.line = element_line(color = "#AAAAAA"),
    axis.ticks = element_line(color = "#AAAAAA"),
    axis.text = element_text(color = "#AAAAAA"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(margin = margin(b = -10)))
}

for(subdf in pt.compare) {
  # two dataframes, one for nuclear-encoded and one for organelle-encoded genes
  df.nu = subdf[subdf$Compartment == "NU",]
  df.pt = subdf[subdf$Compartment == "PT",]
  ref = length(my.plots)+1
  
  my.plots[[ref]] = ggplot(mapping = aes(x=GC, y=Hydro)) +                                      # plot GC vs hydrophobicity
    geom_hex(data = df.nu, alpha=0.5, col=NA) +                                                 # use geom_hex to reflect point density for nuc genes
    scale_fill_gradientn(colours=c("#DDDDDD","#333333"),name='Frequency',na.value=NA) +		# blue scale for "frequency" (point density)
    new_scale("fill") +
    geom_hex(data = df.pt, alpha = 0.5, col=NA) +						# same for org genes
    scale_fill_gradientn(colours=c("#CCCCFF","#0000FF")) +                                      # red scale
    ggtitle(plot.labels[ref]) +                                                                 # styling elements follow beneath
    theme_classic() +
    theme(axis.line = element_line(color = "#AAAAAA"),
    axis.ticks = element_line(color = "#AAAAAA"),
    axis.text = element_text(color = "#AAAAAA"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(margin = margin(b = -10)))
}

res.factor = 3
png(hexesoutput.gc, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)

grid.arrange(my.plots[[1]], my.plots[[2]], my.plots[[3]], my.plots[[4]], my.plots[[5]], my.plots[[6]], my.plots[[7]], my.plots[[8]], my.plots[[9]], my.plots[[10]], my.plots[[11]], my.plots[[12]], ncol=3)

dev.off()

################ now hydro-pka1

my.plots = list()
for(subdf in mt.compare) {
  # two dataframes, one for nuclear-encoded and one for organelle-encoded genes
  df.nu = subdf[subdf$Compartment == "NU",]
  df.mt = subdf[subdf$Compartment == "MT",]
  ref = length(my.plots)+1
  
  my.plots[[ref]] = ggplot(mapping = aes(x=pKa1, y=Hydro)) +                                      # plot GC vs hydrophobicity
    geom_hex(data = df.nu, col=NA) +                                                 # use geom_hex to reflect point density for nuc genes
    scale_fill_gradientn(colours=c("#DDDDDD22","#333333FF"),name='Frequency',na.value=NA) +		# blue scale for "frequency" (point density)
    new_scale("fill") +
    geom_hex(data = df.mt, col=NA) +						# same for org genes
    scale_fill_gradientn(colours=c("#FFCCCC22","#FF0000FF")) +                                      # red scale
    ggtitle(plot.labels[ref]) +                                                                 # styling elements follow beneath
    theme_classic() +
    theme(axis.line = element_line(color = "#AAAAAA"),
    axis.ticks = element_line(color = "#AAAAAA"),
    axis.text = element_text(color = "#AAAAAA"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(margin = margin(b = -10)))
}

for(subdf in pt.compare) {
  # two dataframes, one for nuclear-encoded and one for organelle-encoded genes
  df.nu = subdf[subdf$Compartment == "NU",]
  df.pt = subdf[subdf$Compartment == "PT",]
  ref = length(my.plots)+1
  
  my.plots[[ref]] = ggplot(mapping = aes(x=pKa1, y=Hydro)) +                                      # plot GC vs hydrophobicity
    geom_hex(data = df.nu, col=NA) +                                                 # use geom_hex to reflect point density for nuc genes
    scale_fill_gradientn(colours=c("#DDDDDD22","#333333FF"),name='Frequency',na.value=NA) +		# blue scale for "frequency" (point density)
    new_scale("fill") +
    geom_hex(data = df.pt, col=NA) +						# same for org genes
    scale_fill_gradientn(colours=c("#CCCCFF22","#0000FFFF")) +                                      # red scale
    ggtitle(plot.labels[ref]) +                                                                 # styling elements follow beneath
    theme_classic() +
    theme(axis.line = element_line(color = "#AAAAAA"),
    axis.ticks = element_line(color = "#AAAAAA"),
    axis.text = element_text(color = "#AAAAAA"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(margin = margin(b = -10)))
}

res.factor = 3
png(hexesoutput.pka1, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)

grid.arrange(my.plots[[1]], my.plots[[2]], my.plots[[3]], my.plots[[4]], my.plots[[5]], my.plots[[6]], my.plots[[7]], my.plots[[8]], my.plots[[9]], my.plots[[10]], my.plots[[11]], my.plots[[12]], ncol=3)

dev.off()

res.factor=3
png(hexesoutput.limited, width=300*res.factor, height=400*res.factor, res=72*res.factor)
grid.arrange(my.plots[[1]], my.plots[[7]], ncol=1)
dev.off()

################ now hydro-pka2

my.plots = list()
for(subdf in mt.compare) {
  # two dataframes, one for nuclear-encoded and one for organelle-encoded genes
  df.nu = subdf[subdf$Compartment == "NU",]
  df.mt = subdf[subdf$Compartment == "MT",]
  ref = length(my.plots)+1
  
  my.plots[[ref]] = ggplot(mapping = aes(x=pKa2, y=Hydro)) +                                      # plot GC vs hydrophobicity
    geom_hex(data = df.nu, alpha=0.5, col=NA) +                                                 # use geom_hex to reflect point density for nuc genes
    scale_fill_gradientn(colours=c("#DDDDDD","#333333"),name='Frequency',na.value=NA) +		# blue scale for "frequency" (point density)
    new_scale("fill") +
    geom_hex(data = df.mt, alpha = 0.5, col=NA) +						# same for org genes
    scale_fill_gradientn(colours=c("#FFCCCC","#FF0000")) +                                      # red scale
    ggtitle(plot.labels[ref]) +                                                                 # styling elements follow beneath
    theme_classic() +
    theme(axis.line = element_line(color = "#AAAAAA"),
    axis.ticks = element_line(color = "#AAAAAA"),
    axis.text = element_text(color = "#AAAAAA"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(margin = margin(b = -10)))
}

for(subdf in pt.compare) {
  # two dataframes, one for nuclear-encoded and one for organelle-encoded genes
  df.nu = subdf[subdf$Compartment == "NU",]
  df.pt = subdf[subdf$Compartment == "PT",]
  ref = length(my.plots)+1
  
  my.plots[[ref]] = ggplot(mapping = aes(x=pKa2, y=Hydro)) +                                      # plot GC vs hydrophobicity
    geom_hex(data = df.nu, alpha=0.5, col=NA) +                                                 # use geom_hex to reflect point density for nuc genes
    scale_fill_gradientn(colours=c("#DDDDDD","#333333"),name='Frequency',na.value=NA) +		# blue scale for "frequency" (point density)
    new_scale("fill") +
    geom_hex(data = df.pt, alpha = 0.5, col=NA) +						# same for org genes
    scale_fill_gradientn(colours=c("#CCCCFF","#0000FF")) +                                      # red scale
    ggtitle(plot.labels[ref]) +                                                                 # styling elements follow beneath
    theme_classic() +
    theme(axis.line = element_line(color = "#AAAAAA"),
    axis.ticks = element_line(color = "#AAAAAA"),
    axis.text = element_text(color = "#AAAAAA"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "white"),
    legend.position = "none",
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    plot.title = element_text(margin = margin(b = -10)))
}

res.factor=3
png(hexesoutput.pka2, width=1000*res.factor, height=1000*res.factor, res=72*res.factor)

grid.arrange(my.plots[[1]], my.plots[[2]], my.plots[[3]], my.plots[[4]], my.plots[[5]], my.plots[[6]], my.plots[[7]], my.plots[[8]], my.plots[[9]], my.plots[[10]], my.plots[[11]], my.plots[[12]], ncol=3)

dev.off()

##################

message("Normal trees...")

results = data.frame(complex=NULL,model=NULL,training.acc=NULL,test.acc=NULL,balance=NULL)

#### usual predictors

png(treeplotoutput, width=800, height=800)
par(mfrow=c(4,3))

for(j in 1:length(mt.regex)) {
  working.df = mt.compare[[j]]
  working.df = working.df[working.df$Compartment != "PT",]
  mt.test.acc = mt.training.acc = NULL
  for(i in 1:nsamp) {
    mt.sample.n = nrow(working.df)
    mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
    mt.training.set = working.df[mt.training.refs,]
    mt.test.set = working.df[-mt.training.refs,]
    mt.trained.tree = tree(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, mt.training.set)
    mt.training.predictions = predict(mt.trained.tree, mt.training.set)
    mt.test.predictions = predict(mt.trained.tree, mt.test.set)
    mt.test.predicted.class = ifelse(mt.test.predictions[,which(colnames(mt.test.predictions) == "MT")] < 0.5, "NU", "MT")
    mt.training.predictions = predict(mt.trained.tree, mt.training.set)
    mt.training.predicted.class = ifelse(mt.training.predictions[,which(colnames(mt.training.predictions) == "MT")] < 0.5, "NU", "MT")
    mt.training.accuracy = sum(mt.training.predicted.class == mt.training.set$Compartment)/nrow(mt.training.set)
    mt.training.acc = c(mt.training.acc, mt.training.accuracy)
    mt.test.accuracy = sum(mt.test.predicted.class == mt.test.set$Compartment)/nrow(mt.test.set)
    mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=mt.regex[j], model="tree",training.acc=mean(mt.training.acc, na.rm=T), test.acc=mean(mt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="MT")/nrow(working.df)))
  plot(mt.trained.tree)
  text(mt.trained.tree, pretty=1)
  title(mt.regex[j])
}

for(j in 1:length(pt.regex)) {
  working.df = pt.compare[[j]]
  working.df = working.df[working.df$Compartment != "MT",]
  pt.test.acc = pt.training.acc = NULL
  for(i in 1:nsamp) {
    pt.sample.n = nrow(working.df)
    pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
    pt.training.set = working.df[pt.training.refs,]
    pt.test.set = working.df[-pt.training.refs,]
    pt.trained.tree = tree(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, pt.training.set)
    pt.training.predictions = predict(pt.trained.tree, pt.training.set)
    pt.test.predictions = predict(pt.trained.tree, pt.test.set)
    pt.test.predicted.class = ifelse(pt.test.predictions[,which(colnames(pt.test.predictions) == "PT")] < 0.5, "NU", "PT")
    pt.training.predictions = predict(pt.trained.tree, pt.training.set)
    pt.training.predicted.class = ifelse(pt.training.predictions[,which(colnames(pt.training.predictions) == "PT")] < 0.5, "NU", "PT")
    pt.training.accuracy = sum(pt.training.predicted.class == pt.training.set$Compartment)/nrow(pt.training.set)
    pt.training.acc = c(pt.training.acc, pt.training.accuracy)
    pt.test.accuracy = sum(pt.test.predicted.class == pt.test.set$Compartment)/nrow(pt.test.set)
    pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=pt.regex[j], model="tree",training.acc=mean(pt.training.acc, na.rm=T), test.acc=mean(pt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="PT")/nrow(working.df)))
  plot(pt.trained.tree)
  text(pt.trained.tree, pretty=1)
  title(pt.regex[j])
}
dev.off()

message("Reduced trees...")

#### reduced

png(treereducedplotoutput, width=3000, height=600)
par(mfrow=c(2,max(length(mt.regex), length(pt.regex))))


for(j in 1:length(mt.regex)) {
  working.df = mt.compare[[j]]
  working.df = working.df[working.df$Compartment != "PT",]
  mt.test.acc = mt.training.acc = NULL
  for(i in 1:nsamp) {
    mt.sample.n = nrow(working.df)
    mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
    mt.training.set = working.df[mt.training.refs,]
    mt.test.set = working.df[-mt.training.refs,]
    mt.trained.tree = tree(Compartment ~ Hydro + pKa1, mt.training.set)
    mt.training.predictions = predict(mt.trained.tree, mt.training.set)
    mt.test.predictions = predict(mt.trained.tree, mt.test.set)
    mt.test.predicted.class = ifelse(mt.test.predictions[,which(colnames(mt.test.predictions) == "MT")] < 0.5, "NU", "MT")
    mt.training.predictions = predict(mt.trained.tree, mt.training.set)
    mt.training.predicted.class = ifelse(mt.training.predictions[,which(colnames(mt.training.predictions) == "MT")] < 0.5, "NU", "MT")
    mt.training.accuracy = sum(mt.training.predicted.class == mt.training.set$Compartment)/nrow(mt.training.set)
    mt.training.acc = c(mt.training.acc, mt.training.accuracy)
    mt.test.accuracy = sum(mt.test.predicted.class == mt.test.set$Compartment)/nrow(mt.test.set)
    mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=mt.regex[j], model="tree-reduced",training.acc=mean(mt.training.acc, na.rm=T), test.acc=mean(mt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="MT")/nrow(working.df)))
  plot(mt.trained.tree)
  text(mt.trained.tree, pretty=1)
  title(mt.regex[j])
}

for(j in 1:length(pt.regex)) {
  working.df = pt.compare[[j]]
  working.df = working.df[working.df$Compartment != "MT",]
  pt.test.acc = pt.training.acc = NULL
  for(i in 1:nsamp) {
    pt.sample.n = nrow(working.df)
    pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
    pt.training.set = working.df[pt.training.refs,]
    pt.test.set = working.df[-pt.training.refs,]
    pt.trained.tree = tree(Compartment ~ Hydro + pKa1, pt.training.set)
    pt.training.predictions = predict(pt.trained.tree, pt.training.set)
    pt.test.predictions = predict(pt.trained.tree, pt.test.set)
    pt.test.predicted.class = ifelse(pt.test.predictions[,which(colnames(pt.test.predictions) == "PT")] < 0.5, "NU", "PT")
    pt.training.predictions = predict(pt.trained.tree, pt.training.set)
    pt.training.predicted.class = ifelse(pt.training.predictions[,which(colnames(pt.training.predictions) == "PT")] < 0.5, "NU", "PT")
    pt.training.accuracy = sum(pt.training.predicted.class == pt.training.set$Compartment)/nrow(pt.training.set)
    pt.training.acc = c(pt.training.acc, pt.training.accuracy)
    pt.test.accuracy = sum(pt.test.predicted.class == pt.test.set$Compartment)/nrow(pt.test.set)
    pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=pt.regex[j], model="tree-reduced",training.acc=mean(pt.training.acc, na.rm=T), test.acc=mean(pt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="PT")/nrow(working.df)))
  plot(pt.trained.tree)
  text(pt.trained.tree, pretty=1)
  title(pt.regex[j])
}
dev.off()

message("Total trees cross organelles...")

png(treecrossplotoutput, width=1000, height=600)
par(mfrow=c(2, 2))

all.mt = mt.compare[[1]]
for(j in 2:length(mt.regex)) {
  all.mt = rbind(all.mt, mt.compare[[j]])
}
all.pt = pt.compare[[1]]
for(j in 2:length(pt.regex)) {
  all.pt = rbind(all.pt, pt.compare[[j]])
}

mt.cross.tree = tree(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, all.mt)
mt.cross.training.predictions = predict(mt.cross.tree, all.mt)
mt.cross.test.predictions = predict(mt.cross.tree, all.pt)
mt.cross.training.predicted.class = ifelse(mt.cross.training.predictions[,which(colnames(mt.cross.training.predictions) == "MT")] < 0.5, "NU", "MT")
mt.cross.test.predicted.class = ifelse(mt.cross.test.predictions[,which(colnames(mt.cross.test.predictions) == "MT")] < 0.5, "NU", "PT")
mt.cross.training.accuracy = sum(mt.cross.training.predicted.class == all.mt$Compartment)/nrow(all.mt)
mt.cross.test.accuracy = sum(mt.cross.test.predicted.class == all.pt$Compartment)/nrow(all.pt)
results = rbind(results, data.frame(complex="allPT", model="tree-cross",training.acc=mt.cross.training.accuracy, test.acc=mt.cross.test.accuracy, balance=sum(all.pt$Compartment!="PT")/nrow(all.pt)))

plot(mt.cross.tree)
text(mt.cross.tree, pretty=1)
title("MT predict PT")
  
pt.cross.tree = tree(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, all.pt)
pt.cross.training.predictions = predict(pt.cross.tree, all.pt)
pt.cross.test.predictions = predict(pt.cross.tree, all.mt)
pt.cross.training.predicted.class = ifelse(pt.cross.training.predictions[,which(colnames(pt.cross.training.predictions) == "PT")] < 0.5, "NU", "PT")
pt.cross.test.predicted.class = ifelse(pt.cross.test.predictions[,which(colnames(pt.cross.test.predictions) == "PT")] < 0.5, "NU", "MT")
pt.cross.training.accuracy = sum(pt.cross.training.predicted.class == all.pt$Compartment)/nrow(all.pt)
pt.cross.test.accuracy = sum(pt.cross.test.predicted.class == all.mt$Compartment)/nrow(all.mt)
results = rbind(results, data.frame(complex="allMT", model="tree-cross",training.acc=pt.cross.training.accuracy, test.acc=pt.cross.test.accuracy, balance=sum(all.mt$Compartment!="MT")/nrow(all.mt)))

plot(pt.cross.tree)
text(pt.cross.tree, pretty=1)
title("PT predict MT")

mt.cross.tree = tree(Compartment ~ Hydro + pKa1, all.mt)
mt.cross.training.predictions = predict(mt.cross.tree, all.mt)
mt.cross.test.predictions = predict(mt.cross.tree, all.pt)
mt.cross.training.predicted.class = ifelse(mt.cross.training.predictions[,which(colnames(mt.cross.training.predictions) == "MT")] < 0.5, "NU", "MT")
mt.cross.test.predicted.class = ifelse(mt.cross.test.predictions[,which(colnames(mt.cross.test.predictions) == "MT")] < 0.5, "NU", "PT")
mt.cross.training.accuracy = sum(mt.cross.training.predicted.class == all.mt$Compartment)/nrow(all.mt)
mt.cross.test.accuracy = sum(mt.cross.test.predicted.class == all.pt$Compartment)/nrow(all.pt)
results = rbind(results, data.frame(complex="allPT", model="tree-cross-reduced",training.acc=mt.cross.training.accuracy, test.acc=mt.cross.test.accuracy, balance=sum(all.pt$Compartment!="PT")/nrow(all.pt)))

plot(mt.cross.tree)
text(mt.cross.tree, pretty=1)
title("MT reduced predict PT")

pt.cross.tree = tree(Compartment ~ Hydro + pKa1, all.pt)
pt.cross.training.predictions = predict(pt.cross.tree, all.pt)
pt.cross.test.predictions = predict(pt.cross.tree, all.mt)
pt.cross.training.predicted.class = ifelse(pt.cross.training.predictions[,which(colnames(pt.cross.training.predictions) == "PT")] < 0.5, "NU", "PT")
pt.cross.test.predicted.class = ifelse(pt.cross.test.predictions[,which(colnames(pt.cross.test.predictions) == "PT")] < 0.5, "NU", "MT")
pt.cross.training.accuracy = sum(pt.cross.training.predicted.class == all.pt$Compartment)/nrow(all.pt)
pt.cross.test.accuracy = sum(pt.cross.test.predicted.class == all.mt$Compartment)/nrow(all.mt)
results = rbind(results, data.frame(complex="allMT", model="tree-cross-reduced",training.acc=pt.cross.training.accuracy, test.acc=pt.cross.test.accuracy, balance=sum(all.mt$Compartment!="MT")/nrow(all.mt)))

plot(pt.cross.tree)
text(pt.cross.tree, pretty=1)
title("PT reduced predict MT")

dev.off()

############## random forests for classification

message("The random forest calculations may take several minutes.")
message("Normal RFs...")

#### usual predictors

png(rfplotoutput, width=3000, height=600)
par(mfrow=c(2, max(length(mt.regex), length(pt.regex))))

for(j in 1:length(mt.regex)) {
  working.df = mt.compare[[j]]
  working.df = droplevels(working.df[working.df$Compartment != "PT",])
  mt.test.acc = mt.training.acc = NULL
  for(i in 1:nsamp) {
    mt.sample.n = nrow(working.df)
    mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
    mt.training.set = working.df[mt.training.refs,]
    mt.test.set = working.df[-mt.training.refs,]
    mt.trained.rf = randomForest(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, mt.training.set)
    mt.training.predictions = predict(mt.trained.rf, mt.training.set)
    mt.test.predictions = predict(mt.trained.rf, mt.test.set)
    mt.test.predicted.class = mt.test.predictions
    mt.training.predictions = predict(mt.trained.rf, mt.training.set)
    mt.training.predicted.class = mt.training.predictions
    mt.training.accuracy = sum(mt.training.predicted.class == mt.training.set$Compartment)/nrow(mt.training.set)
    mt.training.acc = c(mt.training.acc, mt.training.accuracy)
    mt.test.accuracy = sum(mt.test.predicted.class == mt.test.set$Compartment)/nrow(mt.test.set)
    mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=mt.regex[j], model="RF",training.acc=mean(mt.training.acc, na.rm=T), test.acc=mean(mt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="MT")/nrow(working.df)))
  varImpPlot(mt.trained.rf, main=mt.regex[j])
}

for(j in 1:length(pt.regex)) {
  working.df = pt.compare[[j]]
  working.df = droplevels(working.df[working.df$Compartment != "MT",])
  pt.test.acc = pt.training.acc = NULL
  for(i in 1:nsamp) {
    pt.sample.n = nrow(working.df)
    pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
    pt.training.set = working.df[pt.training.refs,]
    pt.test.set = working.df[-pt.training.refs,]
    pt.trained.rf = randomForest(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, pt.training.set)
    pt.training.predictions = predict(pt.trained.rf, pt.training.set)
    pt.test.predictions = predict(pt.trained.rf, pt.test.set)
    pt.test.predicted.class = pt.test.predictions
    pt.training.predictions = predict(pt.trained.rf, pt.training.set)
    pt.training.predicted.class = pt.training.predictions
    pt.training.accuracy = sum(pt.training.predicted.class == pt.training.set$Compartment)/nrow(pt.training.set)
    pt.training.acc = c(pt.training.acc, pt.training.accuracy)
    pt.test.accuracy = sum(pt.test.predicted.class == pt.test.set$Compartment)/nrow(pt.test.set)
    pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=pt.regex[j], model="RF",training.acc=mean(pt.training.acc, na.rm=T), test.acc=mean(pt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="PT")/nrow(working.df)))
  varImpPlot(pt.trained.rf, main=pt.regex[j])
}
dev.off()

message("Reduced RFs...")

#### reduced

png(rfreducedplotoutput, width=3000, height=600)
par(mfrow=c(2,max(length(mt.regex), length(pt.regex))))


for(j in 1:length(mt.regex)) {
  working.df = mt.compare[[j]]
  working.df = droplevels(working.df[working.df$Compartment != "PT",])
  mt.test.acc = mt.training.acc = NULL
  for(i in 1:nsamp) {
    mt.sample.n = nrow(working.df)
    mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
    mt.training.set = working.df[mt.training.refs,]
    mt.test.set = working.df[-mt.training.refs,]
    mt.trained.rf = randomForest(Compartment ~ Hydro + pKa1, mt.training.set)
    mt.training.predictions = predict(mt.trained.rf, mt.training.set)
    mt.test.predictions = predict(mt.trained.rf, mt.test.set)
    mt.test.predicted.class = mt.test.predictions
    mt.training.predictions = predict(mt.trained.rf, mt.training.set)
    mt.training.predicted.class = mt.training.predictions
    mt.training.accuracy = sum(mt.training.predicted.class == mt.training.set$Compartment)/nrow(mt.training.set)
    mt.training.acc = c(mt.training.acc, mt.training.accuracy)
    mt.test.accuracy = sum(mt.test.predicted.class == mt.test.set$Compartment)/nrow(mt.test.set)
    mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=mt.regex[j], model="RF-reduced",training.acc=mean(mt.training.acc, na.rm=T), test.acc=mean(mt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="MT")/nrow(working.df)))
  varImpPlot(mt.trained.rf, main=mt.regex[j])
}

for(j in 1:length(pt.regex)) {
  working.df = pt.compare[[j]]
  working.df = droplevels(working.df[working.df$Compartment != "MT",])
  pt.test.acc = pt.training.acc = NULL
  for(i in 1:nsamp) {
    pt.sample.n = nrow(working.df)
    pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
    pt.training.set = working.df[pt.training.refs,]
    pt.test.set = working.df[-pt.training.refs,]
    pt.trained.rf = randomForest(Compartment ~ Hydro + pKa1, pt.training.set)
    pt.training.predictions = predict(pt.trained.rf, pt.training.set)
    pt.test.predictions = predict(pt.trained.rf, pt.test.set)
    pt.test.predicted.class = pt.test.predictions
    pt.training.predictions = predict(pt.trained.rf, pt.training.set)
    pt.training.predicted.class = pt.training.predictions
    pt.training.accuracy = sum(pt.training.predicted.class == pt.training.set$Compartment)/nrow(pt.training.set)
    pt.training.acc = c(pt.training.acc, pt.training.accuracy)
    pt.test.accuracy = sum(pt.test.predicted.class == pt.test.set$Compartment)/nrow(pt.test.set)
    pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  }
  results = rbind(results, data.frame(complex=pt.regex[j], model="RF-reduced",training.acc=mean(pt.training.acc, na.rm=T), test.acc=mean(pt.test.acc, na.rm=T), balance=sum(working.df$Compartment!="PT")/nrow(working.df)))
  varImpPlot(pt.trained.rf, main=pt.regex[j])
}

dev.off()

message("Total RFs cross organelles...")

all.mt = mt.compare[[1]]
for(j in 2:length(mt.regex)) {
  all.mt = rbind(all.mt, mt.compare[[j]])
}
all.pt = pt.compare[[1]]
for(j in 2:length(pt.regex)) {
  all.pt = rbind(all.pt, pt.compare[[j]])
}

# downsample to save memory
all.mt = all.mt[sample(nrow(all.mt), 0.1*nrow(all.mt)),]
all.pt = all.pt[sample(nrow(all.pt), 0.1*nrow(all.pt)),]
all.mt$Compartment = droplevels(all.mt$Compartment)
all.mt$Compartment = as.factor(ifelse(all.mt$Compartment == "MT", "O", "NU"))
all.pt$Compartment = droplevels(all.pt$Compartment)
all.pt$Compartment = as.factor(ifelse(all.pt$Compartment == "PT" | all.pt$Compartment == "CP", "O", "NU"))

mt.cross.rf = randomForest(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, droplevels(all.mt))
mt.cross.training.predictions = predict(mt.cross.rf, all.mt)
mt.cross.test.predictions = predict(mt.cross.rf, all.pt)
mt.cross.training.predicted.class = mt.cross.training.predictions
mt.cross.test.predicted.class = mt.cross.test.predictions
mt.cross.training.accuracy = sum(mt.cross.training.predicted.class == all.mt$Compartment)/nrow(all.mt)
mt.cross.test.accuracy = sum(mt.cross.test.predicted.class == all.pt$Compartment)/nrow(all.pt)
results = rbind(results, data.frame(complex="allPT", model="RF-cross",training.acc=mt.cross.training.accuracy, test.acc=mt.cross.test.accuracy, balance=sum(all.pt$Compartment!="PT")/nrow(all.pt)))

pt.cross.rf = randomForest(Compartment ~ . - Compartment - GeneLabel - Species - GC - GC12 - GC3 - Robust - Uni1 - Uni2, droplevels(all.pt))
pt.cross.training.predictions = predict(pt.cross.rf, all.pt)
pt.cross.test.predictions = predict(pt.cross.rf, all.mt)
pt.cross.training.predicted.class = pt.cross.training.predictions
pt.cross.test.predicted.class = pt.cross.test.predictions
pt.cross.training.accuracy = sum(pt.cross.training.predicted.class == all.pt$Compartment)/nrow(all.pt)
pt.cross.test.accuracy = sum(pt.cross.test.predicted.class == all.mt$Compartment)/nrow(all.mt)
results = rbind(results, data.frame(complex="allMT", model="RF-cross",training.acc=pt.cross.training.accuracy, test.acc=pt.cross.test.accuracy, balance=sum(all.mt$Compartment!="MT")/nrow(all.mt)))

png(rfcrossplotoutput, width=800, height=400)
par(mfrow=c(1,2))
varImpPlot(mt.cross.rf)
varImpPlot(pt.cross.rf)
dev.off()  
  
mt.cross.rf = randomForest(Compartment ~ Hydro + pKa1, droplevels(all.mt))
mt.cross.training.predictions = predict(mt.cross.rf, all.mt)
mt.cross.test.predictions = predict(mt.cross.rf, all.pt)
mt.cross.training.predicted.class = mt.cross.training.predictions
mt.cross.test.predicted.class = mt.cross.test.predictions
mt.cross.training.accuracy = sum(mt.cross.training.predicted.class == all.mt$Compartment)/nrow(all.mt)
mt.cross.test.accuracy = sum(mt.cross.test.predicted.class == all.pt$Compartment)/nrow(all.pt)
results = rbind(results, data.frame(complex="allPT", model="RF-cross-reduced",training.acc=mt.cross.training.accuracy, test.acc=mt.cross.test.accuracy, balance=sum(all.pt$Compartment!="PT")/nrow(all.pt)))

pt.cross.rf = randomForest(Compartment ~ Hydro + pKa1, droplevels(all.pt))
pt.cross.training.predictions = predict(pt.cross.rf, all.pt)
pt.cross.test.predictions = predict(pt.cross.rf, all.mt)
pt.cross.training.predicted.class = pt.cross.training.predictions
pt.cross.test.predicted.class = pt.cross.test.predictions
pt.cross.training.accuracy = sum(pt.cross.training.predicted.class == all.pt$Compartment)/nrow(all.pt)
pt.cross.test.accuracy = sum(pt.cross.test.predicted.class == all.mt$Compartment)/nrow(all.mt)
results = rbind(results, data.frame(complex="allMT", model="RF-cross-reduced",training.acc=pt.cross.training.accuracy, test.acc=pt.cross.test.accuracy, balance=sum(all.mt$Compartment!="MT")/nrow(all.mt)))

write.csv(results, resultsoutput, row.names=F, quote=F)



