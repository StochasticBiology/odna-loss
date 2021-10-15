sink("Outputs/latex-summaries.tex")

df = read.csv("Data/model-sel-simple-bayeslm-stats.csv", header=T)

cat(paste(c("Average Spearman correlations were $\\rho = ", sprintf("%.2f", df$mt.training[1]), "$ and $\\rho = ", sprintf("%.2f", df$pt.training[1]), "$ for training mt and pt sets respectively, and $\\rho = ", sprintf("%.2f", df$mt.test[1]), "$ and $\\rho = ", sprintf("%.2f", df$pt.test[1]), "$ for test mt and pt sets respectively"), collapse=""))
cat("\n\n")

cat(paste(c("$\\rho = ", sprintf("%.2f", df$mt.cross[1]), "$ for pt predicting mt; $\\rho = ", sprintf("%.2f", df$pt.cross[1]), "$ for mt predicting pt"), collapse=""))
cat("\n\n")

cat(paste(c("Considering only genes encoding subunits of central bioenergetic complexes, correlations were $\\rho = ", sprintf("%.2f", df$mt.training[2]), "$ and $\\rho = ", sprintf("%.2f", df$pt.training[2]), "$ for training mt and pt sets respectively, $\\rho = ", sprintf("%.2f", df$mt.test[2]), "$ and $\\rho = ", sprintf("%.2f", df$pt.test[2]), "$ for test mt and pt sets respectively, and $\\rho = ", sprintf("%.2f", df$mt.cross[2]), "$ for pt predicting mt; $\\rho = ", sprintf("%.2f", df$pt.cross[2]), "$ for mt predicting pt"), collapse=""))
cat("\n\n")

df = read.csv("Data/nuc-org-glm-manual.csv", header=T)
mean.mt = which(df$Species == "Mean MT")
sd.mt = which(df$Species == "SD MT")
mean.pt = which(df$Species == "Mean PT")
sd.pt = which(df$Species == "SD PT")

cat(paste(c("(True Positive/Negative rates: mt TP $", sprintf("%.2f", df$tp[mean.mt]), " \\pm ", sprintf("%.2f", df$tp[sd.mt]), "$, TN $", sprintf("%.2f", df$tn[mean.mt]), " \\pm ", sprintf("%.2f", df$tn[sd.mt]), "$, pt TP $", sprintf("%.2f", df$tp[mean.pt]), " \\pm ", sprintf("%.2f", df$tp[sd.pt]), "$, TN $", sprintf("%.2f", df$tn[mean.pt]), " \\pm ", sprintf("%.2f", df$tn[sd.pt]), "$, mean and s.d. across species)"), collapse=""))

sink()