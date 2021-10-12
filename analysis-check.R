df = read.csv("Data/mt-barcodes-manual.csv", header=T, stringsAsFactor=F)
orgs = c("homo sapiens", "plasmodium vivax", "reclinomonas americana", "arabidopsis thaliana", "zea mays", "fucus vesiculosus", "saccharomyces cerevisiae", "malawimonas californiana", "malawimonas jakobiformis")
tests = df[which(sapply(df, function(y) y %in% orgs)),]
for(i in 1:nrow(tests)) {
  print(paste(c(tests$Species[i], "--", sum(tests[i,2:ncol(tests)]), ":", colnames(tests)[tests[i,]==1]), collapse = " "))
}
df$gsum = rowSums(df[,2:ncol(df)])
df[which(df$gsum == max(df$gsum)),]
df[which(df$gsum == min(df$gsum)),]

dbfilename = "Prelims/MTFull14.txt"
# read Kostas' dataset
kostas.df = read.table(dbfilename, sep="\t", header=T)
mismatch = 0
for(i in 1:nrow(df)) {
  ref = which(kostas.df$Scientific.Name == df$Species[i])
  if(length(ref) > 0) {
    if(kostas.df$NCBIcount[ref] != df$gsum[i] & kostas.df$NCBIcount[ref] != 0) {
      print(paste(c(df$Species[i], kostas.df$NCBIcount[ref], df$gsum[i]), collapse=" "))
      mismatch = mismatch+1
    }
  }
}
print(paste(c(mismatch, "of", nrow(df), "=", mismatch/nrow(df)), collapse = " "))