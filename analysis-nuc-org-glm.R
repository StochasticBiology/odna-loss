#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 8) {
  stop("Need: MT stats file, PT stats file, NU stats file, output label, output path, plot path, test proportion, number of samples")
}

#args = c("Data/mt-stats-manual.csv", "Data/pt-stats-manual.csv", "Data/all-stats.csv", "manual", "Data/", "Plots/", "0.5", "10")

message("Loading libraries...")

library(ggplot2)
library(ggnewscale)
library(gridExtra)
library(tree)
library(BMA)
library(cowplot)

outputplot = paste(c(args[6], "nuc-org-glm-", args[4], ".png"), collapse="")
outputtext = paste(c(args[5], "nuc-org-glm-", args[4], ".csv"), collapse="")

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

feature.relabel = function(s) {
  s = gsub("Length.x", "Len", s)
  s = gsub("Hydro.x", "Hyd", s)
  s = gsub("Hydro_i.x", "HydI", s)
  s = gsub("MolWeight.x", "MW", s)
  s = gsub("pKa1.x", "pK1", s)
  s = gsub("pKa2.x", "pK2", s)
  s = gsub("A_Glu.x", "AG", s)
  s = gsub("CW.x", "CW", s)
  return(s)
}
  
message("Reading data...")

# read pipeline outputs
mt.df = read.csv(args[1], header=T, stringsAsFactor=T)
pt.df = read.csv(args[2], header=T, stringsAsFactor=T)
nuc.df = read.csv(args[3], header=T, stringsAsFactor=T)


################ MT BEHAVIOUR

all.df = rbind(mt.df, nuc.df)
all.df = lengthNormalise(all.df)

all.df$Encoding = ifelse(all.df$Compartment == "NU", 0, 1)


results = NULL
modelset = data.frame(Model = NULL, Prob = NULL)
allmodelset = modelset
species.set = unique(nuc.df$Species)
for(species in species.set) {
  test.df = all.df[all.df$Species==species,]
  glm.x = subset(test.df, select=-c(Species, Compartment, GeneLabel, GC, Uni1, Uni2, Robust, GC12, GC3, Encoding))
  glm.y = test.df$Encoding
  if(length(unique(glm.y)) < 2 | length(glm.y) < 10) {
    print(paste("Only one compartment, or insufficient records, found in ", species))
  } else {
    bglm = bic.glm(glm.x, glm.y, glm.family="binomial")
    posts = as.data.frame(t(bglm$probne0))
    posts$Species = species
    if(length(results) == 0) {
      results = posts
    } else {
      results = rbind(results, posts)
    }
    for(i in 1:length(bglm$label)) {
      allmodelset = rbind(allmodelset, data.frame(Model = bglm$label, Prob = bglm$postprob))
      ref = which(modelset$Model == feature.relabel(bglm$label[i]))
      if(length(ref) != 0) {
        modelset$Prob[ref] = modelset$Prob[ref] + bglm$postprob[i]
      } else {
        modelset = rbind(modelset, data.frame(Model = feature.relabel(bglm$label[i]), Prob = bglm$postprob[i]))
      }
    }
  }
}


maxresults = 6
mt.top10 = modelset[order(-modelset$Prob)[1:maxresults],]
alltoplot = allmodelset[which(allmodelset$Model %in% mt.top10$Model),]

alltoplot$Model = gsub(",", "\n", alltoplot$Model)
mt.top10$Model = gsub(",", "\n", mt.top10$Model)
mt.top10$Prob = mt.top10$Prob/length(unique(results$Species))

mt.tmp.valid.df = data.frame(Species=NULL, tp=NULL, tn=NULL, fp=NULL, fn=NULL)
species.set = results$Species
for(species in species.set) {
  current.df = all.df[all.df$Species==species,]
  if(length(unique(current.df$Encoding)) < 2) {
    print(paste("Only one compartment found in ", species))
  } else {
    org.refs = which(current.df$Encoding == 1)
    nuc.refs = which(current.df$Encoding == 0)
    test.refs = c(sample(org.refs, ceiling(length(org.refs)/2)), sample(nuc.refs, ceiling(length(nuc.refs)/2)))
    training.set = current.df[test.refs,]
    test.set = current.df[-test.refs,]
    
    fglm = glm(Encoding ~ Hydro + pKa1, data = training.set, family="binomial")
    training.predictions = predict(fglm, type="response")
    test.predictions = predict(fglm, newdata = test.set, type="response")
    plot(training.set$Encoding, training.predictions)
    points(test.set$Encoding, test.predictions, col="blue")
    
    tp = sum(training.set$Encoding == 1 & training.predictions > 0.5)/ceiling(length(org.refs)/2)
    tn = sum(training.set$Encoding == 0 & training.predictions < 0.5)/ceiling(length(nuc.refs)/2)
    fp = sum(training.set$Encoding == 0 & training.predictions > 0.5)/ceiling(length(nuc.refs)/2)
    fn = sum(training.set$Encoding == 1 & training.predictions < 0.5)/ceiling(length(org.refs)/2)

   mt.tmp.valid.df = rbind(mt.tmp.valid.df, data.frame(Species = species, tp=tp, tn=tn, fp=fp, fn=fn))
   }
   }

mt.valid.df = rbind(data.frame(Species=mt.tmp.valid.df$Species, Stat="TP", Value=mt.tmp.valid.df$tp), data.frame(Species=mt.tmp.valid.df$Species, Stat="TN", Value=mt.tmp.valid.df$tn), data.frame(Species=mt.tmp.valid.df$Species, Stat="FP", Value=mt.tmp.valid.df$fp), data.frame(Species=mt.tmp.valid.df$Species, Stat="FN", Value=mt.tmp.valid.df$fn))

################## PT BEHAVIOUR

all.df = rbind(pt.df, nuc.df)
all.df = lengthNormalise(all.df)

all.df = all.df[all.df$Compartment != "MT",]

all.df$Encoding = ifelse(all.df$Compartment == "NU", 0, 1)


results = NULL
modelset = data.frame(Model = NULL, Prob = NULL)
allmodelset = modelset
species.set = unique(nuc.df$Species)
for(species in species.set) {
  test.df = all.df[all.df$Species==species,]
  glm.x = subset(test.df, select=-c(Species, Compartment, GeneLabel, GC, Uni1, Uni2, Robust, GC12, GC3, Encoding))
  glm.y = test.df$Encoding
  if(length(unique(glm.y)) < 2 | length(glm.y) < 10) {
    print(paste("Only one compartment, or insufficient records, found in ", species))
  } else {
    bglm = bic.glm(glm.x, glm.y, glm.family="binomial")
    posts = as.data.frame(t(bglm$probne0))
    posts$Species = species
    if(length(results) == 0) {
      results = posts
    } else {
      results = rbind(results, posts)
    }
    for(i in 1:length(bglm$label)) {
      allmodelset = rbind(allmodelset, data.frame(Model = bglm$label, Prob = bglm$postprob))
      ref = which(modelset$Model == feature.relabel(bglm$label[i]))
      if(length(ref) != 0) {
        modelset$Prob[ref] = modelset$Prob[ref] + bglm$postprob[i]
      } else {
        modelset = rbind(modelset, data.frame(Model = feature.relabel(bglm$label[i]), Prob = bglm$postprob[i]))
      }
    }
  }
}


pt.top10 = modelset[order(-modelset$Prob)[1:maxresults],]
alltoplot = allmodelset[which(allmodelset$Model %in% pt.top10$Model),]

alltoplot$Model = gsub(",", "\n", alltoplot$Model)
pt.top10$Model = gsub(",", "\n", pt.top10$Model)
pt.top10$Prob = pt.top10$Prob/length(unique(results$Species))



pt.tmp.valid.df = data.frame(Species=NULL, tp=NULL, tn=NULL, fp=NULL, fn=NULL)
species.set = results$Species
for(species in species.set) {
  current.df = all.df[all.df$Species==species,]
  if(length(unique(current.df$Encoding)) < 2) {
    print(paste("Only one compartment found in ", species))
  } else {
    org.refs = which(current.df$Encoding == 1)
    nuc.refs = which(current.df$Encoding == 0)
    test.refs = c(sample(org.refs, ceiling(length(org.refs)/2)), sample(nuc.refs, ceiling(length(nuc.refs)/2)))
    training.set = current.df[test.refs,]
    test.set = current.df[-test.refs,]
    
    fglm = glm(Encoding ~ Hydro + pKa1, data = training.set, family="binomial")
    training.predictions = predict(fglm, type="response")
    test.predictions = predict(fglm, newdata = test.set, type="response")
    plot(training.set$Encoding, training.predictions)
    points(test.set$Encoding, test.predictions, col="blue")
    
    tp = sum(training.set$Encoding == 1 & training.predictions > 0.5)/ceiling(length(org.refs)/2)
    tn = sum(training.set$Encoding == 0 & training.predictions < 0.5)/ceiling(length(nuc.refs)/2)
    fp = sum(training.set$Encoding == 0 & training.predictions > 0.5)/ceiling(length(nuc.refs)/2)
    fn = sum(training.set$Encoding == 1 & training.predictions < 0.5)/ceiling(length(org.refs)/2)

   pt.tmp.valid.df = rbind(pt.tmp.valid.df, data.frame(Species = species, tp=tp, tn=tn, fp=fp, fn=fn))
   }
}

pt.valid.df = rbind(data.frame(Species=pt.tmp.valid.df$Species, Stat="TP", Value=pt.tmp.valid.df$tp), data.frame(Species=pt.tmp.valid.df$Species, Stat="TN", Value=pt.tmp.valid.df$tn), data.frame(Species=pt.tmp.valid.df$Species, Stat="FP", Value=pt.tmp.valid.df$fp), data.frame(Species=pt.tmp.valid.df$Species, Stat="FN", Value=pt.tmp.valid.df$fn))

mt.top10.plot = ggplot(mt.top10, aes(x = factor(Model, levels=Model), y = Prob)) +
  geom_col(fill="#FF8888", colour="#000000") +
  theme_light() + theme(axis.text.x = element_text(size = 10, margin = unit(c(t = -5, r = 0, b = 0, l = 0), "cm"))) +
  xlab("Model structure") + ylab("Posterior probability")
pt.top10.plot = ggplot(pt.top10, aes(x = factor(Model, levels=Model), y = Prob)) +
  geom_col(fill="#8888FF", colour="#000000") +
  theme_light() + theme(axis.text.x = element_text(size = 10, margin = unit(c(t = -5, r = 0, b = 0, l = 0), "cm"))) +
  xlab("Model structure") + ylab("Posterior probability")

mt.valid.plot = ggplot(mt.valid.df, aes(x=Stat, y=Value)) +
  geom_jitter(width=0.25) +
  theme_light() +
  xlab("") + ylab("Rate")
pt.valid.plot = ggplot(pt.valid.df, aes(x=Stat, y=Value)) +
  geom_jitter(width=0.25) +
  theme_light() +
  xlab("") + ylab("Rate")

res.factor=3
png(outputplot, width=1000*res.factor,height=200*res.factor, res=72*res.factor)
plot_grid(mt.top10.plot, pt.top10.plot, mt.valid.plot, pt.valid.plot, nrow=1, align="h")
dev.off()

all.output = rbind(mt.tmp.valid.df, data.frame(Species="Mean MT", tp = mean(mt.tmp.valid.df$tp), tn = mean(mt.tmp.valid.df$tn), fp = mean(mt.tmp.valid.df$fp), fn = mean(mt.tmp.valid.df$fn)), data.frame(Species="SD MT", tp = sd(mt.tmp.valid.df$tp), tn = sd(mt.tmp.valid.df$tn), fp = sd(mt.tmp.valid.df$fp), fn = sd(mt.tmp.valid.df$fn)), pt.tmp.valid.df, data.frame(Species="Mean PT", tp = mean(pt.tmp.valid.df$tp), tn = mean(pt.tmp.valid.df$tn), fp = mean(pt.tmp.valid.df$fp), fn = mean(pt.tmp.valid.df$fn)), data.frame(Species="SD PT", tp = sd(pt.tmp.valid.df$tp), tn = sd(pt.tmp.valid.df$tn), fp = sd(pt.tmp.valid.df$fp), fn = sd(pt.tmp.valid.df$fn)))
write.csv(all.output, outputtext, row.names=F, quote=F)




