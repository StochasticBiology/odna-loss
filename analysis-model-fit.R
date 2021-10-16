#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 10) {
  stop("Need: MT stats, occurence, indices; PT stats, occurrence, indices; average type; plot and data directories, and label")
}


#args = c("Data/mt-stats-means-manual.csv", "Data/mt-gene-occurrence-manual.csv", "Data/mt-simple-manual-indices.csv", "Data/pt-stats-means-manual.csv", "Data/pt-gene-occurrence-manual.csv", "Data/pt-simple-manual-indices.csv", "AltSampled", "Plots/", "Data/", "simple-alt")
#args = c("Data/mt-stats-means-manual.csv", "Data/mt-gene-occurrence-manual.csv", "Data/mt-simple-manual-indices.csv", "Data/pt-stats-means-manual.csv", "Data/pt-gene-occurrence-manual.csv", "Data/pt-simple-manual-indices.csv", "Sampled", "Plots/", "Data/", "simple")

message("Loading libraries...")

library(mombf)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(cowplot)

use.ranks = F              # whether to rank retention indices or use straight value
clean.genes = T            # whether to remove msh1, matk, etc
mean.protocol = args[7]    

message("Getting data...")

mt = read.csv(args[1], header=T, stringsAsFactor=F)
occur.mt = read.csv(args[2], header=T, stringsAsFactor=F)
index.mt = read.csv(args[3], header=T, stringsAsFactor=F)
pt = read.csv(args[4], header=T, stringsAsFactor=F)
occur.pt = read.csv(args[5], header=T, stringsAsFactor=F)
index.pt = read.csv(args[6], header=T, stringsAsFactor=F)
output.plot = paste(c(args[8], "model-sel-", args[10], "-full-plot.png"), collapse="")
output.reduced.plot = paste(c(args[8], "model-sel-", args[10], "-reduced-plot.png"), collapse="")
output.others.plot = paste(c(args[8], "model-sel-", args[10], "-others-plot.png"), collapse="")
output.text = paste(c(args[9], "model-sel-", args[10], "-bayeslm-stats.csv"), collapse="")

# get stats from the averaging protocol we need, and add index values

mt.df = subset(mt[mt$Protocol==mean.protocol & !is.na(mt$Length),], select=-Protocol)
mt.df$Occurrence = NULL
for(i in 1:nrow(mt.df)) {
  mt.df$Occurrence = occur.mt$Occurrence[occur.mt$GeneLabel == mt.df$GeneLabel[i]]
}
mt.df$Index = index.mt$Index[match(mt.df$GeneLabel, index.mt$GeneLabel)]

pt.df = subset(pt[pt$Protocol==mean.protocol & !is.na(pt$Length),], select=-Protocol)
pt.df$Occurrence = NULL
for(i in 1:nrow(pt.df)) {
  pt.df$Occurrence = occur.pt$Occurrence[occur.pt$GeneLabel == pt.df$GeneLabel[i]]
}
pt.df$Index = index.pt$Index[match(pt.df$GeneLabel, index.pt$GeneLabel)]

if(use.ranks == T) {
  mt.df$Index = rank(mt.df$Index)
  pt.df$Index = rank(pt.df$Index)
}

if(clean.genes == T) {
  clean.refs = grep("matr|mttb|msh1|muts", mt.df$GeneLabel)
  if(length(clean.refs) > 0) {
    mt.df = mt.df[-clean.refs,]
  }
}

################ first do model selection

# settings, and data structure, for bayesian model selection

nf = ncol(mt.df)
priorCoef <- momprior(tau=0.348)
priorDelta <- modelbbprior(1,1)

mt.x = mt.df[,1:nf]
mt.y = mt.df$Index
pt.x = pt.df[,1:nf]
pt.y = pt.df$Index

# perform model selection and grab posterior probabilities and signs of coefficients

mt.fit <- modelSelection(mt.y ~ mt.x[,1]+mt.x[,2]+mt.x[,3]+mt.x[,4]+mt.x[,5]+mt.x[,6]+mt.x[,7]+mt.x[,8]+mt.x[,9]+mt.x[,10], priorCoef=priorCoef, priorDelta=priorDelta)
pt.fit <- modelSelection(pt.y ~ pt.x[,1]+pt.x[,2]+pt.x[,3]+pt.x[,4]+pt.x[,5]+pt.x[,6]+pt.x[,7]+pt.x[,8]+pt.x[,9]+pt.x[,10], priorCoef=priorCoef, priorDelta=priorDelta)

feature.labels = model.fit.labels()

mt.pp = postProb(mt.fit)
pt.pp = postProb(pt.fit)
mt.pp$modelid = as.character(mt.pp$modelid)
pt.pp$modelid = as.character(pt.pp$modelid)

mt.signs = ifelse(coef(mt.fit)[,1] < 0, "-", "+")
pt.signs = ifelse(coef(pt.fit)[,1] < 0, "-", "+")

# convert model labels to human-readable format using names of features from dataframe, and get top 10

mt.results = data.frame(Model=NULL,pp=NULL)
maxresults = 6 #length(pp$modelid)
for(i in 1:maxresults) {
  this.label = NULL
  this.result = strsplit(mt.pp$modelid[i], ",")[[1]]
  if(length(this.result) > 0) {
    for(j in 1:length(this.result)) {
      this.trait = as.numeric(this.result[j])
      this.label = c(this.label, paste(c(feature.labels[this.trait], mt.signs[this.trait]), collapse=""))
    }
  }
  mt.results = rbind(mt.results, data.frame(Model=paste(this.label,collapse="\n"), pp=mt.pp$pp[i]))
}

pt.results = data.frame(Model=NULL,pp=NULL)
for(i in 1:maxresults) {
  this.label = NULL
  this.result = strsplit(pt.pp$modelid[i], ",")[[1]]
  if(length(this.result) > 0) {
    for(j in 1:length(this.result)) {
      this.trait = as.numeric(this.result[j])
      this.label = c(this.label, paste(c(feature.labels[this.trait], pt.signs[this.trait]), collapse=""))
    }
  }
  pt.results = rbind(pt.results, data.frame(Model=paste(this.label,collapse="\n"), pp=pt.pp$pp[i]))
}

# produce subplots for final figure

mt.model.sel.plot = ggplot(mt.results, aes(x=factor(Model, levels=Model),y=mt.pp$pp[1:maxresults])) +
  geom_col(fill="#FF8888", colour="#000000") +
  theme_light() + theme(axis.text.x = element_text(size = 10, margin = unit(c(t = -4, r = 0, b = 0, l = 0), "cm"))) +
  xlab("Model structure") + ylab("Posterior probability")
pt.model.sel.plot = ggplot(pt.results, aes(x=factor(Model, levels=Model),y=pt.pp$pp[1:maxresults])) +
  geom_col(fill="#8888FF", colour="#000000") +
  theme_light() + theme(axis.text.x = element_text(size = 10, margin = unit(c(t = -4, r = 0, b = 0, l = 0), "cm"))) +
  xlab("Model structure") + ylab("Posterior probability")

################ now do some model validation

# same structure as elsewhere; do lots of training-test splits and compute correlations between model prediction and test values for each, reporting summary stats
# here we focus on Hydro+GC model structure

mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
nsamp = 100
testpropn = 0.5
cor.method = "spearman"
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  mt.test.set = mt.df[-mt.training.refs,]
  mt.trained.lm = lm(Index ~ Hydro+GC, mt.training.set)
  mt.training.predictions = predict(mt.trained.lm, mt.training.set)
  mt.test.predictions = predict(mt.trained.lm, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  pt.test.set = pt.df[-pt.training.refs,]
  pt.trained.lm = lm(Index ~ Hydro+GC, pt.training.set)
  pt.training.predictions = predict(pt.trained.lm, pt.training.set)
  pt.test.predictions = predict(pt.trained.lm, pt.test.set)

  mt.cross.predictions = predict(pt.trained.lm, mt.df)
  pt.cross.predictions = predict(mt.trained.lm, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index)
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method=cor.method)
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method=cor.method)
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method=cor.method)
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method=cor.method)

  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)
}
mean(mt.training.acc)
mean(mt.test.acc)
mean(pt.training.acc)
mean(pt.test.acc)
mean(mt.cross.acc)
mean(pt.cross.acc)
results = data.frame(method="LM-reduced", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc))

# dataframes for training-test data plots

mt.plot = rbind(data.frame(Predicted = mt.test.predictions, True = mt.test.set$Index, Class = "Test", GeneLabel = mt.test.set$GeneLabel), data.frame(Predicted = mt.training.predictions, True = mt.training.set$Index, Class = "Training", GeneLabel = mt.training.set$GeneLabel))
pt.plot = rbind(data.frame(Predicted = pt.test.predictions, True = pt.test.set$Index, Class = "Test", GeneLabel = pt.test.set$GeneLabel), data.frame(Predicted = pt.training.predictions, True = pt.training.set$Index, Class = "Training", GeneLabel = pt.training.set$GeneLabel))
mt.plot$Class = factor(mt.plot$Class, levels = c("Training", "Test"))
pt.plot$Class = factor(pt.plot$Class, levels = c("Training", "Test"))

mt.cross.plot = data.frame(Predicted = mt.cross.predictions, True = mt.df$Index, GeneLabel = mt.df$GeneLabel) 
pt.cross.plot = data.frame(Predicted = pt.cross.predictions, True = pt.df$Index, GeneLabel = pt.df$GeneLabel) 

# produce plots

mt.model.test.plot = ggplot(mt.plot, aes(x=Predicted, y=True, col=Class, fill=Class)) +
  geom_smooth(fullrange=TRUE, method="lm") +
  geom_point() +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA") +
  theme_light() + theme(legend.position = c(0.15,0.8), legend.background = element_rect(fill=alpha("#FFFFFF", 0.8))) +
  scale_fill_manual(values=c("#AAAAAA", "#FF8888")) +
  scale_color_manual(values=c("#AAAAAA", "#FF8888")) +
  xlab("Predicted retention index") + ylab("True retention index")
pt.model.test.plot = ggplot(pt.plot, aes(x=Predicted, y=True, col=Class, fill=Class)) +
  geom_smooth(fullrange=TRUE, method="lm") +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA") +
  geom_point() +
  theme_light() + theme(legend.position = c(0.15,0.8), legend.background = element_rect(fill=alpha("#FFFFFF", 0.8))) +
  scale_fill_manual(values=c("#AAAAAA", "#8888FF")) +
  scale_color_manual(values=c("#AAAAAA", "#8888FF")) +
  xlab("Predicted retention index") + ylab("True retention index")

mt.model.test.cross.plot = ggplot(mt.cross.plot, aes(x=Predicted, y=True)) +
  geom_smooth(method="lm", fullrange=TRUE, color="#FF8888", fill="#FFCCCC") +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA", color="#888888") +
  geom_point()  +
  theme_light() + xlab("Predicted retention index (from pt fit)") + ylab("True retention index")
pt.model.test.cross.plot = ggplot(pt.cross.plot, aes(x=Predicted, y=True)) +
  geom_smooth(method="lm", fullrange=TRUE, color="#8888FF", fill="#CCCCFF") +
  geom_text_repel(aes(label=GeneLabel), size = 3, segment.color = "#AAAAAA", color="#888888") +
  geom_point()  +
  theme_light() + xlab("Predicted retention index (from mt fit)") + ylab("True retention index")

# lay out overall figure

#png(output.plot, width=800, height=800)
#plot_grid(mt.model.sel.plot, pt.model.sel.plot, mt.model.test.plot, pt.model.test.plot, mt.model.test.cross.plot, pt.model.test.cross.plot, nrow=3, align="v")
#dev.off()

#png(output.plot, width=800, height=400)
res.factor = 3
png(output.plot, width=800*res.factor, height=400*res.factor, res=72*res.factor)
plot_grid(mt.model.sel.plot, mt.model.test.plot, mt.model.test.cross.plot, pt.model.sel.plot, pt.model.test.plot, pt.model.test.cross.plot, nrow=2, align="h")
dev.off()

#grid.arrange(mt.model.sel.plot, pt.model.sel.plot, mt.model.test.plot, pt.model.test.plot, nrow=2)





################# now work with the reduced dataset containing only genes encoding "key" complexes

# grep out subunits of interest
mt.red.df = mt.df[grep("nad|sdh|atp|cox|cytb|rp", mt.df$GeneLabel),]
pt.red.df = pt.df[grep("psa|psb|rp|rbc|ndh|atp|pet", pt.df$GeneLabel),]

# exactly the same model validation process as before

mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
nsamp = 100
testpropn = 0.5
cor.method = "spearman"
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.red.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.red.df[mt.training.refs,]
  mt.test.set = mt.red.df[-mt.training.refs,]
  mt.trained.lm = lm(Index ~ Hydro+GC, mt.training.set)
  mt.training.predictions = predict(mt.trained.lm, mt.training.set)
  mt.test.predictions = predict(mt.trained.lm, mt.test.set)
  
  pt.sample.n = nrow(pt.red.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.red.df[pt.training.refs,]
  pt.test.set = pt.red.df[-pt.training.refs,]
  pt.trained.lm = lm(Index ~ Hydro+GC, pt.training.set)
  pt.training.predictions = predict(pt.trained.lm, pt.training.set)
  pt.test.predictions = predict(pt.trained.lm, pt.test.set)

  mt.cross.predictions = predict(pt.trained.lm, mt.red.df)
  pt.cross.predictions = predict(mt.trained.lm, pt.red.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.red.df$Index)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.red.df$Index)
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method=cor.method)
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method=cor.method)
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method=cor.method)
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method=cor.method)

  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)
}
mean(mt.training.acc)
mean(mt.test.acc)
mean(pt.training.acc)
mean(pt.test.acc)
mean(mt.cross.acc)
mean(pt.cross.acc)
results = rbind(results, data.frame(method="LM-reduced-pruned", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))

mt.plot = rbind(data.frame(Predicted = mt.test.predictions, True = mt.test.set$Index, Class = "Test", GeneLabel = mt.test.set$GeneLabel), data.frame(Predicted = mt.training.predictions, True = mt.training.set$Index, Class = "Training", GeneLabel = mt.training.set$GeneLabel))
pt.plot = rbind(data.frame(Predicted = pt.test.predictions, True = pt.test.set$Index, Class = "Test", GeneLabel = pt.test.set$GeneLabel), data.frame(Predicted = pt.training.predictions, True = pt.training.set$Index, Class = "Training", GeneLabel = pt.training.set$GeneLabel))

mt.model.red.test.plot = ggplot(mt.plot, aes(x=Predicted, y=True, col=Class)) + geom_smooth(method="lm") + geom_point()
pt.model.red.test.plot = ggplot(pt.plot, aes(x=Predicted, y=True, col=Class)) + geom_smooth(method="lm") + geom_point()

res.factor = 3
png(output.reduced.plot, width=800*res.factor, height=800*res.factor, res=72*res.factor)
grid.arrange(mt.model.sel.plot, pt.model.sel.plot, mt.model.red.test.plot, pt.model.red.test.plot, nrow=2)
dev.off()

write.csv(results, output.text, row.names=F, quote=F)

################################## different priors

# settings, and data structure, for bayesian model selection

nf = 14
priorCoef.1 <- imomprior(tau=0.133)
priorDelta.1 <- modelbbprior(1,1)
priorCoef.2 <- momprior(tau=0.348)
priorDelta.2 <- modelunifprior()
priorCoef.3 <- imomprior(tau=0.133)
priorDelta.3 <- modelunifprior()

mt.x = mt.df[,1:nf]
mt.y = mt.df$Index
pt.x = pt.df[,1:nf]
pt.y = pt.df$Index

# perform model selection and grab posterior probabilities and signs of coefficients

mt.fit.1 <- modelSelection(mt.y ~ mt.x[,1]+mt.x[,2]+mt.x[,3]+mt.x[,4]+mt.x[,5]+mt.x[,6]+mt.x[,7]+mt.x[,8]+mt.x[,9], priorCoef=priorCoef.1, priorDelta=priorDelta.1)
pt.fit.1 <- modelSelection(pt.y ~ pt.x[,1]+pt.x[,2]+pt.x[,3]+pt.x[,4]+pt.x[,5]+pt.x[,6]+pt.x[,7]+pt.x[,8]+pt.x[,9], priorCoef=priorCoef.1, priorDelta=priorDelta.1)
mt.fit.2 <- modelSelection(mt.y ~ mt.x[,1]+mt.x[,2]+mt.x[,3]+mt.x[,4]+mt.x[,5]+mt.x[,6]+mt.x[,7]+mt.x[,8]+mt.x[,9], priorCoef=priorCoef.2, priorDelta=priorDelta.2)
pt.fit.2 <- modelSelection(pt.y ~ pt.x[,1]+pt.x[,2]+pt.x[,3]+pt.x[,4]+pt.x[,5]+pt.x[,6]+pt.x[,7]+pt.x[,8]+pt.x[,9], priorCoef=priorCoef.2, priorDelta=priorDelta.2)
mt.fit.3 <- modelSelection(mt.y ~ mt.x[,1]+mt.x[,2]+mt.x[,3]+mt.x[,4]+mt.x[,5]+mt.x[,6]+mt.x[,7]+mt.x[,8]+mt.x[,9], priorCoef=priorCoef.3, priorDelta=priorDelta.3)
pt.fit.3 <- modelSelection(pt.y ~ pt.x[,1]+pt.x[,2]+pt.x[,3]+pt.x[,4]+pt.x[,5]+pt.x[,6]+pt.x[,7]+pt.x[,8]+pt.x[,9], priorCoef=priorCoef.3, priorDelta=priorDelta.3)

mt.pp.1 = postProb(mt.fit.1)
pt.pp.1 = postProb(pt.fit.1)
mt.pp.1$modelid = as.character(mt.pp.1$modelid)
pt.pp.1$modelid = as.character(pt.pp.1$modelid)
mt.pp.2 = postProb(mt.fit.2)
pt.pp.2 = postProb(pt.fit.2)
mt.pp.2$modelid = as.character(mt.pp.2$modelid)
pt.pp.2$modelid = as.character(pt.pp.2$modelid)
mt.pp.3 = postProb(mt.fit.3)
pt.pp.3 = postProb(pt.fit.3)
mt.pp.3$modelid = as.character(mt.pp.3$modelid)
pt.pp.3$modelid = as.character(pt.pp.3$modelid)

mt.results.1 = data.frame(Model = mt.pp.1$modelid[1:maxresults], Prob = mt.pp.1$pp[1:maxresults])
pt.results.1 = data.frame(Model = pt.pp.1$modelid[1:maxresults], Prob = pt.pp.1$pp[1:maxresults])
mt.results.2 = data.frame(Model = mt.pp.2$modelid[1:maxresults], Prob = mt.pp.2$pp[1:maxresults])
pt.results.2 = data.frame(Model = pt.pp.2$modelid[1:maxresults], Prob = pt.pp.2$pp[1:maxresults])
mt.results.3 = data.frame(Model = mt.pp.3$modelid[1:maxresults], Prob = mt.pp.3$pp[1:maxresults])
pt.results.3 = data.frame(Model = pt.pp.3$modelid[1:maxresults], Prob = pt.pp.3$pp[1:maxresults])

# produce subplots for final figure

mt.model.sel.plot.1 = ggplot(mt.results.1, aes(x=factor(Model, levels=Model),y=Prob)) + geom_col() + xlab("Model structure") + ylab("Posterior probability")
pt.model.sel.plot.1 = ggplot(pt.results.1, aes(x=factor(Model, levels=Model),y=Prob)) + geom_col() + xlab("Model structure") + ylab("Posterior probability")
mt.model.sel.plot.2 = ggplot(mt.results.2, aes(x=factor(Model, levels=Model),y=Prob)) + geom_col() + xlab("Model structure") + ylab("Posterior probability")
pt.model.sel.plot.2 = ggplot(pt.results.2, aes(x=factor(Model, levels=Model),y=Prob)) + geom_col() + xlab("Model structure") + ylab("Posterior probability")
mt.model.sel.plot.3 = ggplot(mt.results.3, aes(x=factor(Model, levels=Model),y=Prob)) + geom_col() + xlab("Model structure") + ylab("Posterior probability")
pt.model.sel.plot.3 = ggplot(pt.results.3, aes(x=factor(Model, levels=Model),y=Prob)) + geom_col() + xlab("Model structure") + ylab("Posterior probability")

res.factor=3
png(output.others.plot, width=800*res.factor, height=800*res.factor, res=72*res.factor)
grid.arrange(mt.model.sel.plot.1, pt.model.sel.plot.1, mt.model.sel.plot.2, pt.model.sel.plot.2, mt.model.sel.plot.3, pt.model.sel.plot.3, nrow = 3)
dev.off()