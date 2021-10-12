#library(party)
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2) {
  stop("Need: symbiont stats file, output path")
}

#args = c("Data/symbiont-all-stats.csv", "Plots/")

message("Loading libraries...")
library(ggplot2)
library(tree)
library(gridExtra)
library(ggnewscale)
library(ggpval)
library(ggpubr)

stats.file = args[1]
out.path = args[2]

# function to normalise selected statistics by gene length
# some, like molecular weight and production energy, make more sense unnormalised
lengthNormalise = function(df) {
  df$Hydro = df$Hydro/df$Length
  df$Hydro_i = df$Hydro_i/df$Length
  df$pKa1 = df$pKa1/df$Length
  df$pKa2 = df$pKa2/df$Length
  df$Uni1 = df$Uni1/df$Length
  df$Uni2 = df$Uni2/df$Length
  df$Robust = df$Robust/df$Length
  df$GC = df$GC/df$Length
  df$GC12 = df$GC12/df$Length
  df$GC3 = df$GC3/df$Length
  return(df)
}

df = read.csv(stats.file, header=T)
df = lengthNormalise(df)
df$SystemLabel = as.factor(gsub("Downloads/", "", df$SystemLabel))

out.file = paste(c(out.path, "oo-hydro.png"), collapse="")
res.factor = 3
png(out.file, width=600*res.factor,height=600*res.factor, res=72*res.factor)
hydro.plot = ggplot(df, aes(x=factor(Partner), y=Hydro)) +
  geom_jitter(alpha = 0.1, width=0.25) +
  geom_boxplot(alpha = 0.5) +
  theme_light() +
  xlab("") +
  facet_wrap(~ SystemLabel)
add_pval(hydro.plot)
dev.off()

out.file = paste(c(out.path, "oo-pka1.png"), collapse="")
res.factor = 3
png(out.file, width=600*res.factor,height=600*res.factor, res=72*res.factor)
pKa1.plot = ggplot(df, aes(x=factor(Partner), y=pKa1)) +
  geom_jitter(alpha = 0.1, width=0.25) +
  geom_boxplot(alpha = 0.5) +
  theme_light() +
  xlab("") +
  facet_wrap(~ SystemLabel)
add_pval(pKa1.plot)
dev.off()

out.file = paste(c(out.path, "oo-gc.png"), collapse="")
res.factor = 3
png(out.file, width=600*res.factor,height=600*res.factor, res=72*res.factor)
gc.plot = ggplot(df, aes(x=factor(Partner), y=GC)) +
  geom_jitter(alpha = 0.1, width=0.25) +
  geom_boxplot(alpha = 0.5) +
  theme_light() +
  xlab("") +
  facet_wrap(~ SystemLabel)
add_pval(gc.plot)
dev.off()

reduced.set = c("Azoamicus", "Azolla", "Chromatophore")

df.red = df[df$SystemLabel %in% reduced.set,]
df.red$SystemLabel = droplevels(df.red$SystemLabel)

hydro.reduced.plot = ggplot(df.red, aes(x=factor(Partner), y=Hydro)) +
  geom_jitter(alpha = 0.1, width=0.25) +
  geom_boxplot(alpha = 0.5) +
  theme_light() +
  xlab("") +
  facet_wrap(~ SystemLabel)
hydro.red = add_pval(hydro.reduced.plot)
pKa1.reduced.plot = ggplot(df.red, aes(x=factor(Partner), y=pKa1)) +
  geom_jitter(alpha = 0.1, width=0.25) +
  geom_boxplot(alpha = 0.5) +
  theme_light() +
  xlab("") +
  facet_wrap(~ SystemLabel)
pKa1.red = add_pval(pKa1.reduced.plot)

#res.factor = 3
#out.file = paste(c(out.path, "oo-reduced.png"), collapse="")
#png(out.file, width=600*res.factor,height=400*res.factor, res=72*res.factor)
#grid.arrange(hydro.red, pKa1.red, nrow=2)
#dev.off()

my.round = function(p) {
  if(p < 0.05) {
    return(formatC(p, format="e", digits=1))
  } else {
    return(round(p, digits =2))
  }
}


hydro.new = ggplot(df, aes(x=SystemLabel,y=Hydro, fill=factor(Partner))) +
  geom_point(position=position_jitterdodge(), size=0.2, alpha = 0.1) +
  geom_boxplot(outlier.shape=NA, alpha = 0.75) +
  stat_compare_means(aes(label = my.round(..p..)), method = "wilcox.test", size=3) +
  theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45), axis.title.x=element_blank()) +
  scale_fill_manual(values=c("#8888FF", "#FF8888")) 
pKa1.new = ggplot(df, aes(x=SystemLabel,y=pKa1, fill=factor(Partner))) +
  geom_point(position=position_jitterdodge(), size=0.2, alpha = 0.1) +
  geom_boxplot(outlier.shape=NA, alpha = 0.75) +
  stat_compare_means(aes(label = my.round(..p..)), method = "wilcox.test", size=3) +
  theme_light() + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45), axis.title.x=element_blank()) +
  scale_fill_manual(values=c("#8888FF", "#FF8888")) 

res.factor = 3
out.file = paste(c(out.path, "oo-reduced.png"), collapse="")
png(out.file, width=400*res.factor,height=400*res.factor, res=72*res.factor)

grid.arrange(hydro.new, pKa1.new, nrow=2)

dev.off()

 