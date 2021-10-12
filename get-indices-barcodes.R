#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#args = c("Data/mt-barcodes-manual.csv", "Data/mt-barcode-manual-indices.csv")

if(length(args) < 2) {
  stop("Need barcode file and output file")
}

df = read.csv(args[1], header=T, stringsAsFactor=F)
df = subset(df, select=-Species)

ngenes = ncol(df)
counts = colSums(df)

output.df = data.frame(GeneLabel=colnames(df), Index=log(counts)/max(log(counts)))

write.csv(output.df, args[2], row.names = FALSE, quote= FALSE)

  
