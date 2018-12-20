# Subset AZ, CO, UT, NM plots
plots <- read.csv("D:/EvansLab/Final/Data/FIA/PLOT.csv", header = T, stringsAsFactors = F)
plots_subset <- subset(plots, STATECD %in% c(4, 8, 49, 35))
write.csv(plots_subset, "D:/EvansLab/Final/Data/FIA/PLOT_SUBSET.csv")

# Combine trees
az_trees <- read.csv("D:/EvansLab/Final/Data/FIA/AZ_TREE.csv", header = T, stringsAsFactors = F)
co_trees <- read.csv("D:/EvansLab/Final/Data/FIA/CO_TREE.csv", header = T, stringsAsFactors = F)
nm_trees <- read.csv("D:/EvansLab/Final/Data/FIA/NM_TREE.csv", header = T, stringsAsFactors = F)
ut_trees <- read.csv("D:/EvansLab/Final/Data/FIA/UT_TREE.csv", header = T, stringsAsFactors = F)
trees <- rbind(az_trees, co_trees, nm_trees, ut_trees)
trees$PLT_CN <- as.character(trees$PLT_CN)
write.csv(trees, "D:/EvansLab/Final/Data/FIA/TREE_COMBINED.csv")

# Combine conditions
az_conds <- read.csv("D:/EvansLab/Final/Data/FIA/AZ_COND.csv", header = T, stringsAsFactors = F)
co_conds <- read.csv("D:/EvansLab/Final/Data/FIA/CO_COND.csv", header = T, stringsAsFactors = F)
nm_conds <- read.csv("D:/EvansLab/Final/Data/FIA/NM_COND.csv", header = T, stringsAsFactors = F)
ut_conds <- read.csv("D:/EvansLab/Final/Data/FIA/UT_COND.csv", header = T, stringsAsFactors = F)
conds <- rbind(az_conds, co_conds, nm_conds, ut_conds)
conds$PLT_CN <- as.character(conds$PLT_CN)
write.csv(conds, "D:/EvansLab/Final/Data/FIA/COND_COMBINED.csv")