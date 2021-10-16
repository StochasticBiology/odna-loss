# function to normalise selected statistics by gene length
# some, like molecular weight and production energy, make more sense unnormalised
# depends on the spelling of features in the original data
lengthNormalise = function(df) {
  df$Hydro = df$Hydro/df$Length
  df$Hydro_i = df$Hydro_i/df$Length
  df$pKa1 = df$pKa1/df$Length
  df$pKa2 = df$pKa2/df$Length
  df$pI = df$pI/df$Length
  df$Uni1 = df$Uni1/df$Length
  df$Uni2 = df$Uni2/df$Length
  df$Robust = df$Robust/df$Length
  df$GC = df$GC/df$Length
  df$GC12 = df$GC12/df$Length
  df$GC3 = df$GC3/df$Length
  return(df)
}

# these two functions provide shortened labels for the different features in the model selection process
# one based on relabelling existing labels, one constructed de novo (both depend on the spelling and ordering of features in the original data)
feature.relabel = function(s) {
  s = gsub("Length.x", "Len", s)
  s = gsub("Hydro.x", "Hyd", s)
  s = gsub("Hydro_i.x", "HydI", s)
  s = gsub("MolWeight.x", "MW", s)
  s = gsub("pKa1.x", "pK1", s)
  s = gsub("pKa2.x", "pK2", s)
  s = gsub("pI.x", "pI", s)
  s = gsub("A_Glu.x", "AG", s)
  s = gsub("CW.x", "CW", s)
  return(s)
}

model.fit.labels = function() {
  return(c("C", "Len", "Hyd", "HydI", "MW", "pK1", "pK2", "AG", "CW", "pI", "GC"))
}