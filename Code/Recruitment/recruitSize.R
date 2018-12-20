library(raster)
library(car)
library(rgdal)

plots_subset <- read.csv("D:/EvansLab/Final/Data/FIA/PLOT_SUBSET.csv")

# Create output data frame
out <- data.frame(plot = as.character(plots_subset$CN), 
                  lat = plots_subset$LAT, 
                  lon = plots_subset$LON, 
                  PIED = integer(length(plots_subset$CN)), 
                  measYear = plots_subset$MEASYEAR, 
                  plotPrev = plots_subset$PREV_PLT_CN)
out$plot <- as.character(out$plot)
out$plotPrev <- as.character(out$plotPrev)

# Get rid of plots that were not remeasured
out <- subset(out, !is.na(out$plotPrev))

trees <- read.csv("D:/EvansLab/Final/Data/FIA/TREE_COMBINED.csv")

# Get rid of all plots with no TREE records
out <- out[!is.na(out$PIED),]

# Get recruit stats
out$recruits <- 0
recruits <- subset(trees, RECONCILECD %in% c(1,2) & SPCD == 106)
sizemean <- mean(log(recruits$DRYBIO_AG))
sizesd <- sd(log(recruits$DRYBIO_AG))
save(sizemean, sizesd, file = "D:/EvansLab/Final/Models/B/recrstats.rda")
save(sizemean, sizesd, file = "D:/EvansLab/Final/Models/C/recrstats.rda")
save(sizemean, sizesd, file = "D:/EvansLab/Final/Models/BC/recrstats.rda")
