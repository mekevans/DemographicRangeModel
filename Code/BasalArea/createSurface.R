library(sp)
library(raster)
library(rgdal)
library(automap)

# Read in FIA condition table
cond <- read.csv("D:/EvansLab/Final/Data/FIA/COND_COMBINED.csv", header = T)

# Remove non-forest conditions (missing BALIVE, coded as 0 so causes underestimation)
cond <- subset(cond, COND_STATUS_CD == 1)

# Read in FIA plot table
plots <- read.csv("D:/EvansLab/Final/Data/FIA/PLOT_SUBSET.csv", header = T)
  
# Loop through all plots and calculate plot basal area as mean condition basal area
for (i in 1:nrow(plots)) {
  tmpPlot <- plots[i, "CN"]
  condsMatch <- subset(cond, PLT_CN == tmpPlot)
  balive <- condsMatch$BALIVE
  props <- condsMatch$CONDPROP_UNADJ
  wmean <- weighted.mean(balive, props, na.rm = T)
  print(i)
  # if (length(balive) > 1) {
  #   print(i)
  #   print(wmean)
  # }
  print(wmean)
  if (length(balive) == 0) plots[i, "baLive"] <- NA
  else plots[i, "baLive"] <- wmean
}

# Create output data frame
output <- plots[, c("CN", "LON", "LAT", "baLive")]
output <- output[complete.cases(output),]

# Make output spatial
proj <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"
spOutput <- SpatialPointsDataFrame(output[, c("LON", "LAT")], as.data.frame(output[, "baLive"]), proj4string = CRS(proj))
names(spOutput) <- "BA"

# Write as shapefile
writeOGR(spOutput, "D:/EvansLab/Final/Data/BA", "BasalAreaPoints", "ESRI Shapefile")