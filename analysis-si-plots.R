library(ggplot2)
library(ggrepel)
library(gridExtra)

mt.blast = read.csv("Data/mt-barcodes-blast.csv", header=T)
pt.blast = read.csv("Data/pt-barcodes-blast.csv", header=T)

mt.manual = read.csv("Data/mt-barcodes-manual.csv", header=T)
pt.manual = read.csv("Data/pt-barcodes-manual.csv", header=T)


#### MT

mt.blast.labels = colnames(mt.blast)[2:ncol(mt.blast)]
mt.manual.labels = colnames(mt.manual)[2:ncol(mt.manual)]

mt.both.labels = intersect(mt.blast.labels, mt.manual.labels)

mt.diffs = union(setdiff(mt.blast.labels, mt.manual.labels), setdiff(mt.manual.labels, mt.blast.labels))

mt.compare.df = data.frame(GeneLabel = NULL, ManualCount = NULL, BLASTCount = NULL)
for(i in 1:length(mt.both.labels)) {
  ref.blast = which(colnames(mt.blast) == mt.both.labels[i])
  ref.manual = which(colnames(mt.manual) == mt.both.labels[i])
  mt.compare.df = rbind(mt.compare.df, data.frame(GeneLabel = mt.both.labels[i], ManualCount = sum(mt.manual[,ref.manual]), BLASTCount = sum(mt.blast[,ref.blast])))
}

#### PT

pt.blast.labels = colnames(pt.blast)[2:ncol(pt.blast)]
pt.manual.labels = colnames(pt.manual)[2:ncol(pt.manual)]

pt.both.labels = intersect(pt.blast.labels, pt.manual.labels)

pt.diffs = union(setdiff(pt.blast.labels, pt.manual.labels), setdiff(pt.manual.labels, pt.blast.labels))

pt.compare.df = data.frame(GeneLabel = NULL, ManualCount = NULL, BLASTCount = NULL)
for(i in 1:length(pt.both.labels)) {
  ref.blast = which(colnames(pt.blast) == pt.both.labels[i])
  ref.manual = which(colnames(pt.manual) == pt.both.labels[i])
  pt.compare.df = rbind(pt.compare.df, data.frame(GeneLabel = pt.both.labels[i], ManualCount = sum(pt.manual[,ref.manual]), BLASTCount = sum(pt.blast[,ref.blast])))
}

res.factor = 3
png("Plots/plot-blast-manual.png", width=800*res.factor,height=300*res.factor, res=72*res.factor)
p.mt = ggplot(mt.compare.df, aes(x = log(ManualCount), y = log(BLASTCount))) +
  geom_point() +
  geom_text_repel(aes(label=GeneLabel)) +
  geom_abline(a = 1, b = 0) +
  xlab("log Manual count") + ylab("log BLAST count") +
  theme_light()
p.pt = ggplot(pt.compare.df, aes(x = log(ManualCount), y = log(BLASTCount))) +
  geom_point() +
  geom_text_repel(aes(label=GeneLabel)) +
  geom_abline(a = 1, b = 0) +
  xlab("log Manual count") + ylab("log BLAST count") +
  theme_light()
  
  
grid.arrange(p.mt, p.pt, nrow=1)
dev.off()

sink("Data/blast-manual-diffs.txt")
cat(cor(mt.compare.df$ManualCount, mt.compare.df$BLASTCount))
cat("\n")
cat(cor(pt.compare.df$ManualCount, pt.compare.df$BLASTCount))
cat("\n")
cat(mt.diffs)
cat("\n")
cat(pt.diffs)
sink()




