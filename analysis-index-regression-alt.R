#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 11) {
  stop("Need: MT stats file, MT occurrence file, MT indices file, PT stats file, PT occurrence file, PT indices file, mean protocol, label for outputs, results path, figures path, test proportion, samples")
}

#args = c("Data/mt-stats-means-manual.csv", "Data/mt-gene-occurrence-manual.csv", "Data/mt-simple-manual-indices.csv", "Data/pt-stats-means-manual.csv", "Data/pt-gene-occurrence-manual.csv", "Data/pt-simple-manual-indices.csv", "Sampled", "simple-manual-sampled", "Data/", "Plots/", "0.5", "100")

message("Loading libraries...")

library(logistf)
library(arm)
library(tree)
library(GGally)
library(randomForest)
library(glmnet)
library(phytools)
library(e1071)
set.seed(121)

# lots in here so far!
######### Section 0. get data.
######### Section 1. predictions of different responses: tree for occurrence index, LM for occurrence index
# 1.1 all features
# 1.2 reduced features
######### Section 2. use glmnet for ridge regression and LASSO
######### Section 3. SVR
######### Section 4. use random forests with different feature subsets

mean.protocol = args[7]
testpropn = as.numeric(args[11])
nsamp = as.numeric(args[12])
tree.examples.plot = paste(c(args[10], "index-regression-tree-examples-", args[8], ".png"), collapse="")
rf.examples.plot = paste(c(args[10], "index-regression-rf-examples-", args[8], ".png"), collapse="")

#testpropn = 0.5
#nsamp = 1000

# 90s for nsamp = 100, 690 for nsamp = 1000

startpoint = proc.time()

################################################
#################
######### Section 0. get data.

message("Getting data...")

mt = read.csv(args[1], header=T, stringsAsFactor=F)
occur.mt = read.csv(args[2], header=T, stringsAsFactor=F)
index.mt = read.csv(args[3], header=T, stringsAsFactor=F)
pt = read.csv(args[4], header=T, stringsAsFactor=F)
occur.pt = read.csv(args[5], header=T, stringsAsFactor=F)
index.pt = read.csv(args[6], header=T, stringsAsFactor=F)

#mt = read.csv("Data/mt-stats-manual.csv", header=T, stringsAsFactor=F)
#occur.mt = read.csv("Data/mt-gene-occurrence.csv", header=T, stringsAsFactor=F)
#index.mt = read.csv("Data/mt-indices.csv", header=T, stringsAsFactor=F)
#index.mt = read.csv("Data/mt-indices-network-manual-prune-indices.csv", header=T, stringsAsFactor=F)

#pt = read.csv("Data/pt-stats-manual.csv", header=T, stringsAsFactor=F)
#occur.pt = read.csv("Data/pt-gene-occurrence.csv", header=T, stringsAsFactor=F)
#index.pt = read.csv("Data/pt-indices.csv", header=T, stringsAsFactor=F)
#index.pt = read.csv("Data/pt-indices-network-manual-prune-indices.csv", header=T, stringsAsFactor=F)

outputfilename = paste(c(args[9], "index-regression-", args[8], ".csv"), collapse="")
outputmtpairs = paste(c(args[10], "index-regression-mt-stats-", args[8], ".png"), collapse="")
outputptpairs = paste(c(args[10], "index-regression-pt-stats-", args[8], ".png"), collapse="")

mt.df = subset(mt[mt$Protocol==mean.protocol & !is.na(mt$Length),], select=-Protocol)
mt.df$Occurrence = NULL
for(i in 1:nrow(mt.df)) {
  mt.df$Occurrence = occur.mt$Occurrence[occur.mt$GeneLabel == mt.df$GeneLabel[i]]
}
mt.df$Index = index.mt$Index[match(mt.df$GeneLabel, index.mt$GeneLabel)]

res.factor = 3
png(outputmtpairs, width=1000*res.factor, height=1000*res.factor, res = 72*res.factor)
ggpairs(subset(mt.df, select=-c(Occurrence, GeneLabel)), upper = list(continuous = wrap("cor", size = 5))) + theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
dev.off()

pt.df = subset(pt[pt$Protocol==mean.protocol & !is.na(pt$Length),], select=-Protocol)
pt.df$Occurrence = NULL
for(i in 1:nrow(pt.df)) {
  pt.df$Occurrence = occur.pt$Occurrence[occur.pt$GeneLabel == pt.df$GeneLabel[i]]
}
pt.df$Index = index.pt$Index[match(pt.df$GeneLabel, index.pt$GeneLabel)]

res.factor = 3
png(outputptpairs, width=1000*res.factor, height=1000*res.factor, res = 72*res.factor)
ggpairs(subset(pt.df, select=-c(Occurrence, GeneLabel)),upper = list(continuous = wrap("cor", size = 5))) + theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))
dev.off()

#par(mfrow=c(2,1))
#plot(mt.df$Index, sqrt(mt.df$Occurrence))
#plot(pt.df$Index, sqrt(pt.df$Occurrence))

#pt.df = pt.df[grep("ndh", pt.df$GeneLabel, invert=T),]

results = data.frame(method=NULL, mt.training=NULL, mt.test=NULL, pt.training=NULL, pt.test=NULL, mt.cross=NULL, pt.cross=NULL)

################################################
#################
######### Section 1. predictions of different responses: tree for occurrence index, LM for occurrence index
## Section 1.1. all features

message("Simple models...")

res.factor = 3
png(tree.examples.plot, width=800*res.factor, height=400*res.factor, res = 72*res.factor)
n.plot.tree = 3
par(mfrow=c(2,n.plot.tree))

# full trees for occurrence index
mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  mt.test.set = mt.df[-mt.training.refs,]
  mt.trained.tree = tree(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, mt.training.set)
  mt.training.predictions = predict(mt.trained.tree, mt.training.set)
  mt.test.predictions = predict(mt.trained.tree, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  pt.test.set = pt.df[-pt.training.refs,]
  pt.trained.tree = tree(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, pt.training.set)
  pt.training.predictions = predict(pt.trained.tree, pt.training.set)
  pt.test.predictions = predict(pt.trained.tree, pt.test.set)

  mt.cross.predictions = predict(pt.trained.tree, mt.df)
  pt.cross.predictions = predict(mt.trained.tree, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index)
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)
  if(i <= n.plot.tree) {
    plot(mt.trained.tree)
    text(mt.trained.tree, pretty=1)
    title("MT")
    plot(pt.trained.tree)
    text(pt.trained.tree, pretty=1)
    title("PT")
  }
}
mean(mt.training.acc)
mean(mt.test.acc)
mean(pt.training.acc)
mean(pt.test.acc)
mean(mt.cross.acc)
mean(pt.cross.acc)
results = rbind(results, data.frame(method="Tree", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))
dev.off()



mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = subset(mt.df[mt.training.refs,], select=-GeneLabel)
  mt.test.set = subset(mt.df[-mt.training.refs,], select=-GeneLabel)
  mt.trained.lm = lm(Index ~ . - Index - Occurrence - Uni1 - Uni2 - GC12 - GC3 - Robustness, mt.training.set)
  mt.training.predictions = predict(mt.trained.lm, mt.training.set)
  mt.test.predictions = predict(mt.trained.lm, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = subset(pt.df[pt.training.refs,], select=-GeneLabel)
  pt.test.set = subset(pt.df[-pt.training.refs,], select=-GeneLabel)
  pt.trained.lm = lm(Index ~ . - Index - Occurrence - Uni1 - Uni2 - GC12 - GC3 - Robustness, pt.training.set)
  pt.training.predictions = predict(pt.trained.lm, pt.training.set)
  pt.test.predictions = predict(pt.trained.lm, pt.test.set)

  mt.cross.predictions = predict(pt.trained.lm, mt.df)
  pt.cross.predictions = predict(mt.trained.lm, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index)
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

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
results = rbind(results, data.frame(method="LM", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))


### Section 1.2 reduced simple models

# full trees for occurrence index
mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  mt.test.set = mt.df[-mt.training.refs,]
  mt.trained.tree = tree(Index ~ Hydro + GC, mt.training.set)
  mt.training.predictions = predict(mt.trained.tree, mt.training.set)
  mt.test.predictions = predict(mt.trained.tree, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  pt.test.set = pt.df[-pt.training.refs,]
  pt.trained.tree = tree(Index ~ Hydro + GC, pt.training.set)
  pt.training.predictions = predict(pt.trained.tree, pt.training.set)
  pt.test.predictions = predict(pt.trained.tree, pt.test.set)

  mt.cross.predictions = predict(pt.trained.tree, mt.df)
  pt.cross.predictions = predict(mt.trained.tree, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index)
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)

  if(i <= 4) {
    plot(mt.trained.tree)
    text(mt.trained.tree, pretty=1)
    title("MT")
    plot(pt.trained.tree)
    text(pt.trained.tree, pretty=1)
    title("PT")
  }
}
mean(mt.training.acc)
mean(mt.test.acc)
mean(pt.training.acc)
mean(pt.test.acc)
mean(mt.cross.acc)
mean(pt.cross.acc)
results = rbind(results, data.frame(method="Tree-reduced", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))


mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = subset(mt.df[mt.training.refs,], select=-GeneLabel)
  mt.test.set = subset(mt.df[-mt.training.refs,], select=-GeneLabel)
  mt.trained.lm = lm(Index ~ Hydro + GC, mt.training.set)
  mt.training.predictions = predict(mt.trained.lm, mt.training.set)
  mt.test.predictions = predict(mt.trained.lm, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = subset(pt.df[pt.training.refs,], select=-GeneLabel)
  pt.test.set = subset(pt.df[-pt.training.refs,], select=-GeneLabel)
  pt.trained.lm = lm(Index ~ Hydro + GC, pt.training.set)
  pt.training.predictions = predict(pt.trained.lm, pt.training.set)
  pt.test.predictions = predict(pt.trained.lm, pt.test.set)

  mt.cross.predictions = predict(pt.trained.lm, mt.df)
  pt.cross.predictions = predict(mt.trained.lm, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index)
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index)
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

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
results = rbind(results, data.frame(method="LM-Reduced", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))

################################################
#################
######### Section 2. use glmnet for ridge regression and LASSO
## Section 2.1. MT+PT

message("Ridge and LASSO...")

mt.reduced.predictors = as.matrix(subset(mt.df, select=-c(Index, Occurrence, GeneLabel, Uni1, Uni2, GC12, GC3, Robustness)))
mt.reduced.response = as.matrix(subset(mt.df, select=Index))
pt.reduced.predictors = as.matrix(subset(pt.df, select=-c(Index, Occurrence, GeneLabel, Uni1, Uni2, GC12, GC3, Robustness)))
pt.reduced.response = as.matrix(subset(pt.df, select=Index))

# ridge regression for occurrence index

mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)

  mt.training.predictors = mt.reduced.predictors[mt.training.refs,]
  mt.test.predictors = mt.reduced.predictors[-mt.training.refs,]
  mt.training.response = mt.reduced.response[mt.training.refs]
  mt.test.response = mt.reduced.response[-mt.training.refs]
  mt.lambdas <- 10^seq(2, -3, by = -.1)
  mt.cv_ridge <- cv.glmnet(mt.training.predictors, mt.training.response, alpha = 0, lambda = mt.lambdas)
  mt.optimal_lambda <- mt.cv_ridge$lambda.min
  mt.ridge.fit = glmnet(mt.training.predictors, mt.training.response, lambda=mt.optimal_lambda)

  pt.training.predictors = pt.reduced.predictors[pt.training.refs,]
  pt.test.predictors = pt.reduced.predictors[-pt.training.refs,]
  pt.training.response = pt.reduced.response[pt.training.refs]
  pt.test.response = pt.reduced.response[-pt.training.refs]
  pt.lambdas <- 10^seq(2, -3, by = -.1)
  pt.cv_ridge <- cv.glmnet(pt.training.predictors, pt.training.response, alpha = 0, lambda = pt.lambdas)
  pt.optimal_lambda <- pt.cv_ridge$lambda.min
  pt.ridge.fit = glmnet(pt.training.predictors, pt.training.response, lambda=pt.optimal_lambda)

  mt.training.predictions = predict(mt.ridge.fit, mt.training.predictors, s = mt.optimal_lambda)
  mt.test.predictions = predict(mt.ridge.fit, mt.test.predictors, s = mt.optimal_lambda)
  pt.training.predictions = predict(pt.ridge.fit, pt.training.predictors, s = pt.optimal_lambda)
  pt.test.predictions = predict(pt.ridge.fit, pt.test.predictors, s = pt.optimal_lambda)
  mt.cross.predictions = predict(pt.ridge.fit, mt.reduced.predictors, s=pt.optimal_lambda)
  pt.cross.predictions = predict(mt.ridge.fit, pt.reduced.predictors, s=mt.optimal_lambda)

  mt.training.accuracy = cor(mt.training.predictions, mt.training.response, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.response, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.response, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.response, method="spearman")
  mt.cross.accuracy = cor(mt.cross.predictions, mt.reduced.response, method="spearman")
  pt.cross.accuracy = cor(pt.cross.predictions, pt.reduced.response, method="spearman")
  
  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)
}
mean(mt.training.acc, na.rm=T)
mean(mt.test.acc, na.rm=T)
mean(pt.training.acc, na.rm=T)
mean(pt.test.acc, na.rm=T)
mean(mt.cross.acc, na.rm=T)
mean(pt.cross.acc, na.rm=T)
results = rbind(results, data.frame(method="Ridge", mt.training=mean(mt.training.acc, na.rm=T), mt.test=mean(mt.test.acc, na.rm=T), pt.training=mean(pt.training.acc, na.rm=T), pt.test=mean(pt.test.acc, na.rm=T), mt.cross=mean(mt.cross.acc, na.rm=T), pt.cross=mean(pt.cross.acc, na.rm=T)))

# LASSO for occurrence index

mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)

  mt.training.predictors = mt.reduced.predictors[mt.training.refs,]
  mt.test.predictors = mt.reduced.predictors[-mt.training.refs,]
  mt.training.response = mt.reduced.response[mt.training.refs]
  mt.test.response = mt.reduced.response[-mt.training.refs]
  mt.lambdas <- 10^seq(2, -3, by = -.1)
  mt.cv_ridge <- cv.glmnet(mt.training.predictors, mt.training.response, alpha = 1, lambda = mt.lambdas, standardize = TRUE, nfolds = 5)
  mt.optimal_lambda <- mt.cv_ridge$lambda.min
  mt.ridge.fit = glmnet(mt.training.predictors, mt.training.response, alpha = 1, lambda=mt.optimal_lambda, standardize = TRUE)

  pt.training.predictors = pt.reduced.predictors[pt.training.refs,]
  pt.test.predictors = pt.reduced.predictors[-pt.training.refs,]
  pt.training.response = pt.reduced.response[pt.training.refs]
  pt.test.response = pt.reduced.response[-pt.training.refs]
  pt.lambdas <- 10^seq(2, -3, by = -.1)
  pt.cv_ridge <- cv.glmnet(pt.training.predictors, pt.training.response, alpha = 1, lambda = pt.lambdas, standardize = TRUE, nfolds = 5)
  pt.optimal_lambda <- pt.cv_ridge$lambda.min
  pt.ridge.fit = glmnet(pt.training.predictors, pt.training.response, alpha = 1, lambda=pt.optimal_lambda, standardize = TRUE)

  mt.training.predictions = predict(mt.ridge.fit, mt.training.predictors, s = mt.optimal_lambda)
  mt.test.predictions = predict(mt.ridge.fit, mt.test.predictors, s = mt.optimal_lambda)
  pt.training.predictions = predict(pt.ridge.fit, pt.training.predictors, s = pt.optimal_lambda)
  pt.test.predictions = predict(pt.ridge.fit, pt.test.predictors, s = pt.optimal_lambda)
  mt.cross.predictions = predict(pt.ridge.fit, mt.reduced.predictors, s=pt.optimal_lambda)
  pt.cross.predictions = predict(mt.ridge.fit, pt.reduced.predictors, s=mt.optimal_lambda)

  mt.training.accuracy = cor(mt.training.predictions, mt.training.response, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.response, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.response, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.response, method="spearman")
  mt.cross.accuracy = cor(mt.cross.predictions, mt.reduced.response, method="spearman")
  pt.cross.accuracy = cor(pt.cross.predictions, pt.reduced.response, method="spearman")
  
  # store statistics
  mt.test.acc = c(mt.test.acc, mt.test.accuracy)
  mt.training.acc = c(mt.training.acc, mt.training.accuracy)
  pt.test.acc = c(pt.test.acc, pt.test.accuracy)
  pt.training.acc = c(pt.training.acc, pt.training.accuracy)
  mt.cross.acc = c(mt.cross.acc, mt.cross.accuracy)
  pt.cross.acc = c(pt.cross.acc, pt.cross.accuracy)
}
mean(mt.training.acc, na.rm=T)
mean(mt.test.acc, na.rm=T)
mean(pt.training.acc, na.rm=T)
mean(pt.test.acc, na.rm=T)
mean(mt.cross.acc, na.rm=T)
mean(pt.cross.acc, na.rm=T)
results = rbind(results, data.frame(method="LASSO", mt.training=mean(mt.training.acc, na.rm=T), mt.test=mean(mt.test.acc, na.rm=T), pt.training=mean(pt.training.acc, na.rm=T), pt.test=mean(pt.test.acc, na.rm=T), mt.cross=mean(mt.cross.acc, na.rm=T), pt.cross=mean(pt.cross.acc, na.rm=T)))

proc.time()-startpoint

################################################
#################
######### Section 3. SVR
## Section 3.1. MT+PT

message("SVR...")

tune = F

mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
mt.importance = pt.importance = 0
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  mt.test.set = mt.df[-mt.training.refs,]
  if(tune == T) {
    mt.trained.set = tune(svm, Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, data=mt.training.set, ranges=list(epsilon=seq(0,1,0.1), cost=1:10))
    mt.trained.svm = mt.trained.set$best.model
  } else {
    mt.trained.svm = svm(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, mt.training.set)
  }
  mt.training.predictions = predict(mt.trained.svm, mt.training.set)
  mt.test.predictions = predict(mt.trained.svm, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  pt.test.set = pt.df[-pt.training.refs,]
  if(tune == T) {
    pt.trained.set = tune(svm, Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, data=pt.training.set, ranges=list(epsilon=seq(0,1,0.1), cost=1:10))
    pt.trained.svm = pt.trained.set$best.model
  } else {
    pt.trained.svm = svm(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, pt.training.set)
  }
  pt.training.predictions = predict(pt.trained.svm, pt.training.set)
  pt.test.predictions = predict(pt.trained.svm, pt.test.set)

  mt.cross.predictions = predict(pt.trained.svm, mt.df)
  pt.cross.predictions = predict(mt.trained.svm, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index, method="spearman")
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index, method="spearman")
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

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

if(tune == T) {
  results = rbind(results, data.frame(method="SVR-tune", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))
} else {
  results = rbind(results, data.frame(method="SVR", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))
}

proc.time()-startpoint

################################################
#################
######### Section 4. random forests for occurrence index
## Section 4.1. MT+PT

message("Random forests...")
res.factor = 3
png(rf.examples.plot, width=1000*res.factor, height=600*res.factor, res = 72*res.factor)
par(mfrow=c(2,2))

mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
mt.importance = pt.importance = 0
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  mt.test.set = mt.df[-mt.training.refs,]
  mt.trained.rf = randomForest(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, mt.training.set)
  mt.training.predictions = predict(mt.trained.rf, mt.training.set)
  mt.test.predictions = predict(mt.trained.rf, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  pt.test.set = pt.df[-pt.training.refs,]
  pt.trained.rf = randomForest(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, pt.training.set)
  pt.training.predictions = predict(pt.trained.rf, pt.training.set)
  pt.test.predictions = predict(pt.trained.rf, pt.test.set)

  mt.cross.predictions = predict(pt.trained.rf, mt.df)
  pt.cross.predictions = predict(mt.trained.rf, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index, method="spearman")
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index, method="spearman")
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

  # store statistics
  mt.importance = mt.importance + importance(mt.trained.rf)
  pt.importance = pt.importance + importance(pt.trained.rf)
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

results = rbind(results, data.frame(method="RF", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))

varImpPlot(mt.trained.rf, main="MT")
varImpPlot(pt.trained.rf, main="PT")

#plot(test.predictions, mt.test.set$Index)

######### 
## Section 4.2. reduced


mt.test.acc = mt.training.acc = mt.cross.acc = pt.test.acc = pt.training.acc = pt.cross.acc = NULL
for(i in 1:nsamp) {
  mt.sample.n = nrow(mt.df)
  mt.training.refs = sample(seq(from=1,to=mt.sample.n), size=mt.sample.n*testpropn)
  mt.training.set = mt.df[mt.training.refs,]
  mt.test.set = mt.df[-mt.training.refs,]
  mt.trained.rf = randomForest(Index ~ Hydro + GC, mt.training.set)
  mt.training.predictions = predict(mt.trained.rf, mt.training.set)
  mt.test.predictions = predict(mt.trained.rf, mt.test.set)
  
  pt.sample.n = nrow(pt.df)
  pt.training.refs = sample(seq(from=1,to=pt.sample.n), size=pt.sample.n*testpropn)
  pt.training.set = pt.df[pt.training.refs,]
  pt.test.set = pt.df[-pt.training.refs,]
  pt.trained.rf = randomForest(Index ~ Hydro + GC, pt.training.set)
  pt.training.predictions = predict(pt.trained.rf, pt.training.set)
  pt.test.predictions = predict(pt.trained.rf, pt.test.set)

  mt.cross.predictions = predict(pt.trained.rf, mt.df)
  pt.cross.predictions = predict(mt.trained.rf, pt.df)
  mt.cross.accuracy = cor(mt.cross.predictions, mt.df$Index, method="spearman")
  pt.cross.accuracy = cor(pt.cross.predictions, pt.df$Index, method="spearman")
  
  mt.training.accuracy = cor(mt.training.predictions, mt.training.set$Index, method="spearman")
  mt.test.accuracy = cor(mt.test.predictions, mt.test.set$Index, method="spearman")
  pt.training.accuracy = cor(pt.training.predictions, pt.training.set$Index, method="spearman")
  pt.test.accuracy = cor(pt.test.predictions, pt.test.set$Index, method="spearman")

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


results = rbind(results, data.frame(method="RF-Reduced", mt.training=mean(mt.training.acc), mt.test=mean(mt.test.acc), pt.training=mean(pt.training.acc), pt.test=mean(pt.test.acc), mt.cross=mean(mt.cross.acc), pt.cross=mean(pt.cross.acc)))


#plot(test.predictions, mt.test.set$Index)
#

#################################
#################################

mt.trained.rf = randomForest(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, mt.df)
pt.trained.rf = randomForest(Index ~ . - Index - Occurrence - GeneLabel - Uni1 - Uni2 - GC12 - GC3 - Robustness, pt.df)

mt.training.predictions = predict(mt.trained.rf, mt.df)
pt.predictions = predict(mt.trained.rf, pt.df)
pt.training.predictions = predict(pt.trained.rf, pt.df)
mt.predictions = predict(pt.trained.rf, mt.df)
mt.accuracy = cor(mt.training.predictions, mt.df$Index, method="spearman")
pt.accuracy = cor(pt.training.predictions, pt.df$Index, method="spearman")
mt.cross.accuracy = cor(mt.predictions, mt.df$Index, method="spearman")
pt.cross.accuracy = cor(pt.predictions, pt.df$Index, method="spearman")

results = rbind(results, data.frame(method="RF-Cross", mt.training=mt.accuracy, mt.test=0, pt.training=pt.accuracy, pt.test=0, mt.cross=mt.cross.accuracy, pt.cross=pt.cross.accuracy))

varImpPlot(mt.trained.rf, main="MT (all genes)")
varImpPlot(pt.trained.rf, main="PT (all genes)")
dev.off()

mt.trained.red.rf = randomForest(Index ~ Hydro + GC, mt.df)
pt.trained.red.rf = randomForest(Index ~ Hydro + GC, pt.df)

mt.training.predictions = predict(mt.trained.red.rf, mt.df)
pt.predictions = predict(mt.trained.red.rf, pt.df)
pt.training.predictions = predict(pt.trained.red.rf, pt.df)
mt.predictions = predict(pt.trained.red.rf, mt.df)
mt.accuracy = cor(mt.training.predictions, mt.df$Index, method="spearman")
pt.accuracy = cor(pt.training.predictions, pt.df$Index, method="spearman")
mt.cross.accuracy = cor(mt.predictions, mt.df$Index, method="spearman")
pt.cross.accuracy = cor(pt.predictions, pt.df$Index, method="spearman")

results = rbind(results, data.frame(method="RF-Cross-Reduced", mt.training=mt.accuracy, mt.test=0, pt.training=pt.accuracy, pt.test=0, mt.cross=mt.cross.accuracy, pt.cross=pt.cross.accuracy))

varImpPlot(mt.trained.rf)
varImpPlot(pt.trained.rf)
varImpPlot(mt.trained.rf)
varImpPlot(pt.trained.rf)
dev.off()

message("Output...")

write.csv(results, outputfilename, row.names=F, quote=F)

proc.time()-startpoint

