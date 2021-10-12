#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 8) {
  stop("Need: complex data file, MT gene occurrence, PT gene occurrence, label aliases, complex labels, output plot file, output model summaries, threshold")
}

#args = c("Data/complex-data.csv", "Data/mt-gene-occurrence.csv", "Data/pt-gene-occurrence.csv", "Prelims/label-aliases.csv", "Prelims/complex-labels.csv", "Plots/", "Data/", "1000")

message("Loading libraries...")

library(arm)
library(blme)
library(lme4)
library(ggplot2)
library(ggrepel)

# threshold of species count for high/low occurrence binarisation
threshold = as.numeric(args[8])

message("Getting data...")

# read pipeline outputs
df = read.csv(args[1], header=T, stringsAsFactor=F)
mt.occur = read.csv(args[2], header=T, stringsAsFactor=F)
pt.occur = read.csv(args[3], header=T, stringsAsFactor=F)

# read identities of genes and complexes
aliases = read.csv(args[4], header=T, stringsAsFactor=F)
complexes = read.csv(args[5], header=T, stringsAsFactor=F)

outputplot.all = paste(c(args[6], "complexes-", args[8], "-all.png"), collapse="")
outputplot.single = paste(c(args[6], "complexes-", args[8], "-single.png"), collapse="")
modeloutput = paste(c(args[7], "complexes-", args[8], "-model-out.txt"), collapse="")

# construct dataframe for analysis and plotting
# introducing occurrence counts and complex labels to the energetic dataframe
df$Occurrence = 0
df$Complex = NA
for(i in 1:nrow(df)) {
  # find and replace gene labels
  ref = which(aliases$PDBLabel == df$GeneLabel[i])
  if(length(ref) == 1) {
    df$GeneLabel[i] = aliases$AltLabel[ref]
  }
  df$GeneLabel[i] = tolower(df$GeneLabel[i])
  # look up occurrence counts
  ref = which(mt.occur$GeneLabel == df$GeneLabel[i])
  if(length(ref) == 1) {
    df$Occurrence[i] = mt.occur$Occurrence[ref]
  }
  # look up complex label
  ref = which(pt.occur$GeneLabel == df$GeneLabel[i])
  if(length(ref) == 1) {
    df$Occurrence[i] = pt.occur$Occurrence[ref]
  }
  ref = which(complexes$PDBLabel == df$PDBLabel[i])
  df$Complex[i] = complexes$Complex[ref]
}

# discard rna from test data set
test.df = df[grep("rna", df$GeneLabel, invert=T),]

# binarise data with high threshold
test.df$Binarised = ifelse(test.df$Occurrence>threshold, 1, 0)

message("Fitting models...")

# bayesian logit model 
bin.model = bayesglm(Binarised ~ MeanEnergywLigand*Complex, data = test.df, family=binomial)
test.df$fit.bin.model = fitted(bin.model)

# bayesian logit mixed model
bin.mixed.model = bglmer(Binarised ~ MeanEnergywLigand + (MeanEnergywLigand | Complex), data=test.df, binomial, cov.prior=NULL)
test.df$fit.bin.mixed.model = fitted(bin.mixed.model)

# compare Occurrence models
ggplot(test.df) +
  geom_point(aes(x = MeanEnergywLigand, y = Binarised)) +
  geom_line(aes(x = MeanEnergywLigand, y = fit.bin.model), color = "black") +
  geom_line(aes(x = MeanEnergywLigand, y = fit.bin.mixed.model), color = "red") +
  facet_wrap(~ Complex)

message("Computing uncertainties (this may take several minutes)...")

fit.frame = data.frame(MeanEnergywLigand = NULL, Complex = NULL)
for(complex in unique(df$Complex)) {
  tmp.frame = data.frame(MeanEnergywLigand = seq(from=-300,to=0, by=10), Complex = complex)
  fit.frame = rbind(fit.frame, tmp.frame)
}
predictions = predict(bin.model, fit.frame, se.fit=TRUE)
linkinv = family(bin.model)$linkinv
raw.predictions = predictions$fit
fit.frame$bin.model.mean = linkinv(predictions$fit)
fit.frame$bin.model.hi = linkinv(raw.predictions+1.96*predictions$se.fit)
fit.frame$bin.model.lo = linkinv(raw.predictions-1.96*predictions$se.fit)

fitfn = function(model) {
   return(predict(model, fit.frame, type="response"))
}
cifn = function(x) {
  tx = x[order(x)]
  lci = round(length(x)*0.05)
  uci = round(length(x)*0.95)
  return(c(tx[lci],tx[uci]))
}
boot<-bootMer(bin.mixed.model, FUN=fitfn, nsim=1000) 
boot.se<-apply(boot$t, 2, cifn)
fit.frame$bin.mixed.model.mean = predict(bin.mixed.model, fit.frame, type="response")
fit.frame$bin.mixed.model.lo = boot.se[1,]
fit.frame$bin.mixed.model.hi = boot.se[2,]

res.factor = 3
png(outputplot.all, width=800*res.factor, height=400*res.factor, res=72*res.factor)

ggplot() + 
# geom_line(aes(x = MeanEnergywLigand, y = fit.bin.mixed.model), color = "#000000") +
  geom_ribbon(data = fit.frame, aes(x = MeanEnergywLigand, ymin = bin.mixed.model.lo, ymax = bin.mixed.model.hi), fill = "#0000FF", alpha = 0.1) +
  geom_ribbon(data = fit.frame, aes(x = MeanEnergywLigand, ymin = bin.model.lo, ymax = bin.model.hi), fill = "#FF0000", alpha = 0.1) + 
  geom_line(data = fit.frame, size = 1, aes(x = MeanEnergywLigand, y = bin.mixed.model.mean), color = "#8888FF") +
  geom_line(data = fit.frame, size = 1, aes(x = MeanEnergywLigand, y = bin.model.mean), color = "#FF8888") +
#  geom_line(data = fit.frame, aes(x = MeanEnergywLigand, y = bin.mixed.model.lo), color = "#8888FF") +
#  geom_line(data = fit.frame, aes(x = MeanEnergywLigand, y = bin.mixed.model.hi), color = "#8888FF") +
  geom_text_repel(data = test.df, aes(x = MeanEnergywLigand, y = Binarised, label=GeneLabel, angle=90), size = 3, segment.color = "#AAAAAA", color="#888888", max.overlaps=20) +
  geom_point(data = test.df, size = 2, alpha = 0.7, aes(x = MeanEnergywLigand, y = Binarised)) +
  scale_y_continuous(name = "Retention", breaks = c(0, 1)) +
  scale_x_continuous(name = "Binding energy / kJmol-1") +
  theme_light() +
  facet_wrap(~ Complex, nrow=2)

dev.off()


relabel = function(str) {
  if(str == "mt-c1") return("MT CI")
  if(str == "mt-c2") return("MT CII")
  if(str == "mt-c3") return("MT CIII")
  if(str == "mt-c4") return("MT CIV")
  if(str == "mt-c5") return("MT CV")
  if(str == "pt-as") return("PT ATP Synthase")
  if(str == "pt-cb") return("PT Cytochrome b6f")
  if(str == "pt-pa") return("PT PSI")
  if(str == "pt-pb") return("PT PSII")
  if(str == "pt-rb") return("PT Rubisco")
  if(str == "pt-ri") return("PT Ribosome")
}

fit.frame$Organelle = "PT"
fit.frame$Organelle[grep("mt", fit.frame$Complex)] = "MT"
fit.frame$Label = unlist(lapply(fit.frame$Complex, relabel))
test.df$Label = unlist(lapply(test.df$Complex, relabel))

res.factor = 3
png(outputplot.single, width=800*res.factor, height=300*res.factor, res=72*res.factor)

ggplot() +
  geom_line(data = fit.frame, size = 1, aes(x = MeanEnergywLigand, y = bin.mixed.model.mean, color = factor(Organelle))) +
  geom_ribbon(data = fit.frame, aes(x = MeanEnergywLigand, ymin = bin.mixed.model.lo, ymax = bin.mixed.model.hi, fill = factor(Organelle)), alpha = 0.1)+  facet_wrap(~ Complex, nrow=2) +
  geom_point(data = test.df, size = 2, alpha = 0.7, aes(x = MeanEnergywLigand, y = Binarised)) +
  geom_text_repel(data = test.df, aes(x = MeanEnergywLigand, y = Binarised, label=GeneLabel, angle=90), size = 2, segment.color = "#AAAAAA", color="#888888", max.overlaps=20) +
  scale_y_continuous(name = "Retention", breaks = c(0, 1)) +
  scale_x_continuous(name = "Binding energy / kJmol-1") +
  theme_light() +
  scale_color_manual(values = c("#FF0000", "#0000FF")) +  scale_fill_manual(values = c("#FF0000", "#0000FF")) +
  facet_wrap(~ Label, nrow=2) +
  theme(legend.position = "none")

dev.off()

sink(modeloutput)
print(summary(bin.model))
print(summary(bin.mixed.model))
sink()


