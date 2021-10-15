library(ggplot2)
library(GGally)

source("lengthNormalise.R")

mt = read.csv("Data/mt-stats-manual.csv", header=T, stringsAsFactors=F)
pt = read.csv("Data/pt-stats-manual.csv", header=T, stringsAsFactors=F)

mt.df = mt[1,]
mt.df = mt.df[-1,]
mt.df$Occurrence = NULL
mt.df$Species = mt.df$Compartment = NULL
for(genelabel in unique(mt$GeneLabel)) {
  subset = mt[mt$GeneLabel == genelabel,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  mt.df = rbind(mt.df, as.data.frame(genestats))
}
mt.df = lengthNormalise(mt.df)

pt.df = pt[1,]
pt.df = pt.df[-1,]
pt.df$Occurrence = NULL
pt.df$Species = pt.df$Compartment = NULL
for(genelabel in unique(pt$GeneLabel)) {
  subset = pt[pt$GeneLabel == genelabel,]
  subset$Species = subset$GeneLabel = subset$Compartment = NULL
  genestats = colMeans(subset)
  genestats$GeneLabel = genelabel
  pt.df = rbind(pt.df, as.data.frame(genestats))
}
pt.df = lengthNormalise(pt.df)

df = read.csv("Data/complex-data.csv", header=T, stringsAsFactor=F)

# read identities of genes and complexes
aliases = read.csv("Prelims/label-aliases.csv", header=T, stringsAsFactor=F)
complexes = read.csv("Prelims/complex-labels.csv", header=T, stringsAsFactor=F)


# construct dataframe for analysis and plotting
# introducing occurrence counts and complex labels to the energetic dataframe
for(i in 1:nrow(df)) {
  # find and replace gene labels
  ref = which(aliases$PDBLabel == df$GeneLabel[i])
  if(length(ref) == 1) {
    df$GeneLabel[i] = aliases$AltLabel[ref]
  }
}

df$GeneLabel = tolower(df$GeneLabel)
mt.df$GeneLabel = tolower(mt.df$GeneLabel)
pt.df$GeneLabel = tolower(pt.df$GeneLabel)

mt.set = data.frame(GeneLabel=intersect(df$GeneLabel, mt.df$GeneLabel), Hydro=0, MeanEnergywLigand=0, EnergywLigand=0, InterfaceswLigand=0)
for(i in 1:nrow(mt.set)) {
  ref.1 = which(mt.df$GeneLabel == mt.set$GeneLabel[i])
  ref.2 = which(df$GeneLabel == mt.set$GeneLabel[i])
  mt.set$Hydro[i] = mt.df$Hydro[i]
  mt.set$MeanEnergywLigand[i] = df$MeanEnergywLigand[i]
  mt.set$EnergywLigand[i] = df$EnergywLigand[i]
  mt.set$InterfaceswLigand[i] = df$InterfaceswLigand[i]
}
pt.set = data.frame(GeneLabel=intersect(df$GeneLabel, pt.df$GeneLabel), Hydro=0, MeanEnergywLigand=0, EnergywLigand=0, InterfaceswLigand=0)
for(i in 1:nrow(pt.set)) {
  ref.1 = which(pt.df$GeneLabel == pt.set$GeneLabel[i])
  ref.2 = which(df$GeneLabel == pt.set$GeneLabel[i])
  pt.set$Hydro[i] = pt.df$Hydro[i]
  pt.set$MeanEnergywLigand[i] = df$MeanEnergywLigand[i]
  pt.set$EnergywLigand[i] = df$EnergywLigand[i]
  pt.set$InterfaceswLigand[i] = df$InterfaceswLigand[i]
}
mt.set$Compartment = "MT"
pt.set$Compartment = "PT"
both.set = rbind(mt.set, pt.set)

#both.test = both.set[both.set$Energy > -1000,]
#ggpairs(subset(both.test, select=-c(GeneLabel, Compartment)))

png("Plots/energy-hydro.png", width=800, height=800)
ggpairs(subset(both.set, select=-c(GeneLabel, Compartment)))
dev.off()