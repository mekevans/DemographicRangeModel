#### download FIA census data from FIA DataMart
#### https://apps.fs.usda.gov/fia/datamart/CSV/datamart_csv.html

### code developed by Michiel Pillet
### modified by M. Evans, L. Huelsmann


plot.path <- "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/plots/"
tree.path <- "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/trees/"
cond.path <- "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/conds/"

# Combine plots
plot_list <- list.files(path=paste(plot.path,sep=''), pattern="*.csv")
plots_combined <- do.call("rbind",
                         lapply(plot_list,
                                function(x)
                                read.csv(paste(plot.path, x, sep=''),
                                           header = T, stringsAsFactors = F)))
write.csv(plots_combined, "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/PLOT_COMBINED.csv")


# Combine trees
tree_list <- list.files(path=paste(tree.path,sep=''), pattern="*.csv")
trees_combined <- do.call("rbind",
                          lapply(tree_list,
                                 function(x)
                                   read.csv(paste(tree.path, x, sep=''),
                                            header = T, stringsAsFactors = F)))
trees_combined$PLT_CN <- as.character(trees_combined$PLT_CN)
write.csv(trees_combined, "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/TREE_COMBINED.csv")

# Combine conditions
cond_list <- list.files(path=paste(cond.path,sep=''), pattern="*.csv")
conds_combined <- do.call("rbind",
                          lapply(cond_list,
                                 function(x)
                                   read.csv(paste(cond.path, x, sep=''),
                                            header = T, stringsAsFactors = F)))
conds_combined$PLT_CN <- as.character(conds_combined$PLT_CN)
write.csv(conds_combined, "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/COND_COMBINED.csv")
