#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 5) {
  stop("Need: phylogeny, NU stats file, MT stats file, PT stats file, NU stats file, output label")
}

message("Loading libraries...")
library(ggplot2)
library(ape)
library(ggtree)
library(phytools)

# function to normalise selected statistics by gene length
# some, like molecular weight and production energy, make more sense unnormalised
lengthNormalise = function(df) {
  df$Hydro = df$Hydro/df$Length
  df$Hydro_i = df$Hydro_i/df$Length
  df$pKa1 = df$pKa1/df$Length
  df$pKa2 = df$pKa2/df$Length
  df$pI = df$pI/df$Length
  df$Uni1 = df$Uni1/df$Length
  df$Uni2 = df$Uni2/df$Length
  df$Robust = df$Robust/df$Length
  df$GC = df$GC/df$Length
  df$GC12 = df$GC12/df$Length
  df$GC3 = df$GC3/df$Length
  return(df)
}

#args = c("Prelims/whole-genome-species.phy", "Data/all-stats.csv", "Data/mt-stats-manual.csv", "Data/pt-stats-manual.csv", "Data/phy-plot")

phy.file = args[1]
nuc.file = args[2]
mt.file = args[3]
pt.file = args[4]
out.label = args[5]

message("Loading phylogeny...")
treeString = tolower(paste(readLines(phy.file), collapse=""))
treeString = gsub(" ", "_", treeString)
tree <- read.tree(text=treeString)
tree$tip.label = gsub("_", " ", tree$tip.label)
tree$node.label = gsub("_", " ", tree$node.label)

metazoa.node = which(tree$node.label=="metazoa")+length(tree$tip.label)
plants.node = which(tree$node.label=="viridiplantae")+length(tree$tip.label)
fungi.node = which(tree$node.label=="fungi")+length(tree$tip.label)
hs.node = which(tree$tip.label=="homo sapiens")
at.node = which(tree$tip.label=="arabidopsis thaliana")
pf.node = which(tree$tip.label=="plasmodium falciparum")
dm.node = which(tree$tip.label=="drosophila melanogaster")
dr.node = which(tree$tip.label=="danio rerio")
ce.node = which(tree$tip.label=="caenorhabditis elegans")
gm.node = which(tree$tip.label=="glycine max")
pp.node = which(tree$tip.label=="physcomitrium patens")
dd.node = which(tree$tip.label=="dictyostelium discoideum")
sp.node = which(tree$tip.label=="schizosaccharomyces pombe")

plot.just.tree = ggtree(tree) +
    theme_tree2() +
    geom_cladelabel(node=metazoa.node, label="Metazoa", color="#000000", offset=10, align=T, angle=90, offset.text = 3) +
    geom_cladelabel(node=plants.node, label="Viridiplantae", color="#000000", offset=10, align=T, angle=90, offset.text = 3) +
    geom_cladelabel(node=fungi.node, label="Fungi", color="#000000", offset=10, align=T, angle=90, offset.text = 3) +
    geom_cladelabel(node=hs.node, label="HS", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=at.node, label="AT", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=pf.node, label="PF", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=dm.node, label="DM", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=dr.node, label="DR", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=ce.node, label="CE", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=gm.node, label="GM", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=pp.node, label="PP", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=dd.node, label="DD", color="#444444", align=F, fontsize=3) +
    geom_cladelabel(node=sp.node, label="SP", color="#444444", align=F, fontsize=3) +

    xlim_tree(60)

plot.label.tree = 
    ggtree(tree) +
    theme_tree2() +
    geom_tiplab(align=F, size=3, color="#888888") +
    geom_cladelabel(node=metazoa.node, label="Metazoa", color="#000088", offset=10, align=T) +
    geom_cladelabel(node=plants.node, label="Viridiplantae", color="#008800", offset=10, align=T) +
    geom_cladelabel(node=fungi.node, label="Fungi", color="#444444", offset=10, align=T) +
    xlim_tree(60)

message("Loading gene data...")
df = read.csv(nuc.file, header=T)
mt.df = read.csv(mt.file, header=T)
pt.df = read.csv(pt.file, header=T)

df = rbind(df, mt.df, pt.df)
df = lengthNormalise(df)

nuc.species = tolower(unique(df$Species[which(df$Compartment == "NU")]))
df$Species = tolower(df$Species)
plotnames = colnames(df)[4:ncol(df)]

df2 = df[which(df$Species %in% tree$tip.label),]
nuc.species = unique(df2$Species)

message("Building comparative plot...")
plot.sd = list()
plot.sem = list()
plot.sd[[1]] = plot.just.tree
plot.sem[[1]] = plot.just.tree
for(j in 4:12) {
#for(j in 4:ncol(df2)) {
  blank = rep(0, length(nuc.species))
  stats = data.frame(species = nuc.species, mean.n = blank, sd.n = blank, n.n = blank, mean.m = blank, sd.m = blank, n.m = blank, mean.p = blank, sd.p = blank, n.p = blank)
  for(i in 1:length(nuc.species))
  {
    vals.n = df2[df2$Species == nuc.species[i] & df2$Compartment == "NU",j]
    vals.m = df2[df2$Species == nuc.species[i] & df2$Compartment == "MT",j]
    vals.p = df2[df2$Species == nuc.species[i] & (df2$Compartment == "CP" | df2$Compartment == "PT"),j]

    stats$mean.n[i] = mean(vals.n)
    stats$sd.n[i] = sd(vals.n)
    stats$n.n[i] = length(vals.n)
    stats$mean.m[i] = mean(vals.m)
    stats$sd.m[i] = sd(vals.m)
    stats$n.m[i] = length(vals.m)
    stats$mean.p[i] = mean(vals.p)
    stats$sd.p[i] = sd(vals.p)
    stats$n.p[i] = length(vals.p)
  }
  
  rownames(stats) = nuc.species
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.n, xmin=mean.n-sd.n, xmax=mean.n+sd.n), color='red3', alpha=0.5)
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.m, xmin=mean.m-sd.m, xmax=mean.m+sd.m), color='blue3', alpha=0.5)
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.p, xmin=mean.p-sd.p, xmax=mean.p+sd.p), color='green3', alpha=0.5)

  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.n, xmin=mean.n-sd.n/sqrt(n.n), xmax=mean.n+sd.n/sqrt(n.n)), color='red3', alpha=0.5)
  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.m, xmin=mean.m-sd.m/sqrt(n.m), xmax=mean.m+sd.m/sqrt(n.m)), color='blue3', alpha=0.5)
  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.p, xmin=mean.p-sd.p/sqrt(n.p), xmax=mean.p+sd.p/sqrt(n.p)), color='green3', alpha=0.5)
}

message("Plotting...")
res.factor = 3
png(paste(c(out.label, "-sd.png"), collapse=""), width=1000*res.factor, height=600*res.factor, res=72*res.factor)
plot.sd[[length(plot.sd)]]
dev.off()

png(paste(c(out.label, "-sem.png"), collapse=""), width=1000*res.factor, height=600*res.factor, res=72*res.factor)
plot.sem[[length(plot.sem)]]
dev.off()

### now with sequence features

message("Building comparative plot with sequence features...")
plot.sd = list()
plot.sem = list()
plot.sd[[1]] = plot.just.tree
plot.sem[[1]] = plot.just.tree
for(j in 4:ncol(df2)) {
  blank = rep(0, length(nuc.species))
  stats = data.frame(species = nuc.species, mean.n = blank, sd.n = blank, n.n = blank, mean.m = blank, sd.m = blank, n.m = blank, mean.p = blank, sd.p = blank, n.p = blank)
  for(i in 1:length(nuc.species))
  {
    vals.n = df2[df2$Species == nuc.species[i] & df2$Compartment == "NU",j]
    vals.m = df2[df2$Species == nuc.species[i] & df2$Compartment == "MT",j]
    vals.p = df2[df2$Species == nuc.species[i] & (df2$Compartment == "CP" | df2$Compartment == "PT"),j]

    stats$mean.n[i] = mean(vals.n)
    stats$sd.n[i] = sd(vals.n)
    stats$n.n[i] = length(vals.n)
    stats$mean.m[i] = mean(vals.m)
    stats$sd.m[i] = sd(vals.m)
    stats$n.m[i] = length(vals.m)
    stats$mean.p[i] = mean(vals.p)
    stats$sd.p[i] = sd(vals.p)
    stats$n.p[i] = length(vals.p)
  }
  
  rownames(stats) = nuc.species
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.n, xmin=mean.n-sd.n, xmax=mean.n+sd.n), color='red3', alpha=0.5)
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.m, xmin=mean.m-sd.m, xmax=mean.m+sd.m), color='blue3', alpha=0.5)
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.p, xmin=mean.p-sd.p, xmax=mean.p+sd.p), color='green3', alpha=0.5)

  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.n, xmin=mean.n-sd.n/sqrt(n.n), xmax=mean.n+sd.n/sqrt(n.n)), color='red3', alpha=0.5)
  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.m, xmin=mean.m-sd.m/sqrt(n.m), xmax=mean.m+sd.m/sqrt(n.m)), color='blue3', alpha=0.5)
  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.p, xmin=mean.p-sd.p/sqrt(n.p), xmax=mean.p+sd.p/sqrt(n.p)), color='green3', alpha=0.5)
}

message("Plotting with sequence features...")
res.factor = 3
png(paste(c(out.label, "-sd-full.png"), collapse=""), width=1000*res.factor, height=600*res.factor, res=72*res.factor)
plot.sd[[length(plot.sd)]]
dev.off()

png(paste(c(out.label, "-sem-full.png"), collapse=""), width=1000*res.factor, height=600*res.factor, res=72*res.factor)
plot.sem[[length(plot.sem)]]
dev.off()


### now with selected features

message("Building comparative plot with sequence features...")
plot.sd = list()
plot.sem = list()
plot.sd[[1]] = plot.just.tree
plot.sem[[1]] = plot.just.tree
for(j in c(5, 8, 12)) {
  blank = rep(0, length(nuc.species))
  stats = data.frame(species = nuc.species, mean.n = blank, sd.n = blank, n.n = blank, mean.m = blank, sd.m = blank, n.m = blank, mean.p = blank, sd.p = blank, n.p = blank)
  for(i in 1:length(nuc.species))
  {
    vals.n = df2[df2$Species == nuc.species[i] & df2$Compartment == "NU",j]
    vals.m = df2[df2$Species == nuc.species[i] & df2$Compartment == "MT",j]
    vals.p = df2[df2$Species == nuc.species[i] & (df2$Compartment == "CP" | df2$Compartment == "PT"),j]

    stats$mean.n[i] = mean(vals.n)
    stats$sd.n[i] = sd(vals.n)
    stats$n.n[i] = length(vals.n)
    stats$mean.m[i] = mean(vals.m)
    stats$sd.m[i] = sd(vals.m)
    stats$n.m[i] = length(vals.m)
    stats$mean.p[i] = mean(vals.p)
    stats$sd.p[i] = sd(vals.p)
    stats$n.p[i] = length(vals.p)
  }
  
  rownames(stats) = nuc.species
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.n, xmin=mean.n-sd.n, xmax=mean.n+sd.n), color='red3', alpha=0.5)
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.m, xmin=mean.m-sd.m, xmax=mean.m+sd.m), color='blue3', alpha=0.5)
  plot.sd[[length(plot.sd)+1]] <- facet_plot(plot.sd[[length(plot.sd)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.p, xmin=mean.p-sd.p, xmax=mean.p+sd.p), color='green3', alpha=0.5)

  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.n, xmin=mean.n-sd.n/sqrt(n.n), xmax=mean.n+sd.n/sqrt(n.n)), color="#888888", alpha=1)
  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.m, xmin=mean.m-sd.m/sqrt(n.m), xmax=mean.m+sd.m/sqrt(n.m)), color="#FF8888", alpha=1)
  plot.sem[[length(plot.sem)+1]] <- facet_plot(plot.sem[[length(plot.sem)]], panel=plotnames[j-3], data=stats, geom=geom_errorbar, aes(x = mean.p, xmin=mean.p-sd.p/sqrt(n.p), xmax=mean.p+sd.p/sqrt(n.p)), color="#8888FF", alpha=1)
}

message("Plotting with selected features...")
res.factor = 3
png(paste(c(out.label, "-sd-select.png"), collapse=""), width=600*res.factor, height=300*res.factor, res=72*res.factor)
plot.sd[[length(plot.sd)]]
dev.off()

png(paste(c(out.label, "-sem-select.png"), collapse=""), width=600*res.factor, height=300*res.factor, res=72*res.factor)
plot.sem[[length(plot.sem)]]
dev.off()
