library(ggplot2)
library(GGally)

df = read.csv("feature-corr.csv", header=T)
ggpairs(subset(df[df$Residue!="*",], select=-c(Codon, Residue)))
