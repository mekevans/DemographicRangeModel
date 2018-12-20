# Read in data
setwd("D:/EvansLab/Final/Data/FIA/")
UT <- read.csv("UT_TREE.csv", header = T, stringsAsFactors = F)
AZ <- read.csv("AZ_TREE.csv", header = T, stringsAsFactors = F)
NM <- read.csv("NM_TREE.csv", header = T, stringsAsFactors = F)
CO <- read.csv("CO_TREE.csv", header = T, stringsAsFactors = F)
all <- rbind(AZ, CO, NM, UT)

# Keep only plot and species code
all <- all[, c("PLT_CN", "SPCD")]

# Read in plots
UT <- read.csv("UT_PLOT.csv", header = T, stringsAsFactors = F)
AZ <- read.csv("AZ_PLOT.csv", header = T, stringsAsFactors = F)
NM <- read.csv("NM_PLOT.csv", header = T, stringsAsFactors = F)
CO <- read.csv("CO_PLOT.csv", header = T, stringsAsFactors = F)
plots <- rbind(AZ, CO, NM, UT)
plots <- plots[, c("CN", "LAT", "LON")]

# Processing
all$LAT <- NA
all$LON <- NA
progress <- 1
no <- length(plots$CN)
for (i in plots$CN) {
  if (nrow(all[all$PLT_CN == i,]) > 0) {
    all[all$PLT_CN == i,]$LAT <- subset(plots, CN == i)$LAT
    all[all$PLT_CN == i,]$LON <- subset(plots, CN == i)$LON
  }
  progress <- progress + 1
  print(progress / no * 100)
}

# Export
setwd("D:/EvansLab/Final/Data/Processed/Validation")
write.csv(all, "validationRecords.csv", row.names = F)

# Read in a raster for alignment
pres <- raster("D:/EvansLab/Final/Output/B/lambda.tif")
pres <- setValues(pres, NA)

# Make records spatial
proj <- "+init=epsg:4269 +proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs +towgs84=0,0,0"
allSp <- SpatialPointsDataFrame(all[, c("LON", "LAT")], as.data.frame(all[, "SPCD"]), proj4string = CRS("+proj=longlat +datum=NAD83"))

# For each tree, find in which raster cell it falls
allSp$cellNumber <- NA
progress <- 1
no <- nrow(allSp)
for (i in 1:nrow(allSp)) {
 tmp <- as.integer(extract(pres, allSp[i,], cellnumbers = T)[[1]])
 allSp[i, "cellNumber"] <- tmp
 progress <- progress + 1
 print(progress / no * 100)
} 

# For each raster cell, find whether it contains PIED (SPCD = 106)
names(allSp) <- c("SPCD", "cellNumber")
PIED <- 106
for (i in 1:ncell(pres)) {
  tmp <- subset(allSp, cellNumber == i)
  if (nrow(tmp) > 0) {
    if (PIED %in% unique(tmp$SPCD)) pres[i] <- 1
    else pres[i] <- 0
  }
}

# Export presence/absence raster
writeRaster(pres, "presenceAbsenceRaster.tif", overwrite = T)