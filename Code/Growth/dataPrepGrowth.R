library(sp)
library(raster)

# Read data and subset PIED
grData <- read.csv("D:/EvansLab/Final/Data/FIA/TREE_COMBINED.csv", header = T, stringsAsFactors = F)
grData <- subset(grData, SPCD == 106)

# Only keep trees that didn't die
grData <- subset(grData, STATUSCD == 1)

# Only keep remeasured trees
grData_remeas <- subset(grData, !is.na(PREVDIA))

# Look up previous AGB
for (i in 1:nrow(grData_remeas)) {
  print(i)
  tmpCN <- grData_remeas[i, "PREV_TRE_CN"]
  tmpSubset <- subset(grData, CN == tmpCN)
  if (nrow(tmpSubset) == 1) {
    grData_remeas[i, "PREV_DRYBIO_AG"] <- tmpSubset$DRYBIO_AG
    grData_remeas[i, "PREV_PLT_CN"] <- tmpSubset$PLT_CN
  }
  else {
    grData_remeas[i, "PREV_DRYBIO_AG"] <- NA
    grData_remeas[i, "PREV_PLT_CN"] <- NA
  }
}

# Subset records with previous AGB found
grData_remeas <- subset(grData_remeas, !is.na(PREV_DRYBIO_AG))
grData_remeas$DRYBIO_AG_DIFF <- grData_remeas$DRYBIO_AG - grData_remeas$PREV_DRYBIO_AG
grData_remeas <- subset(grData_remeas, DRYBIO_AG_DIFF >= 0)

# Read in plot data and get coordinates and previous measurement year
plots <- read.csv("D:/EvansLab/Final/Data/FIA/PLOT_SUBSET.csv", header = T, stringsAsFactors = F)
for (i in 1:nrow(grData_remeas)) {
  print(i)
  # Get latitude and longitude
  tmpPlot <- grData_remeas[i, "PLT_CN"]
  tmpSubset <- subset(plots, CN == tmpPlot)
  grData_remeas[i, "LAT"] <- tmpSubset$LAT
  grData_remeas[i, "LON"] <- tmpSubset$LON
  grData_remeas[i, "MEASYEAR"] <- tmpSubset$MEASYEAR
  # Get previous measurement year
  tmpPlot <- grData_remeas[i, "PREV_PLT_CN"]
  tmpSubset <- subset(plots, CN == tmpPlot)
  grData_remeas[i, "PREV_MEASYEAR"] <- tmpSubset$MEASYEAR
}

# Calculate census interval
grData_remeas$CENSUS_INTERVAL <- grData_remeas$MEASYEAR - grData_remeas$PREV_MEASYEAR

# Get plot BA
conds <- read.csv("D:/EvansLab/Final/Data/FIA/COND_COMBINED.csv", header = T, stringsAsFactors = F)
conds <- subset(conds, COND_STATUS_CD == 1)
for (i in 1:nrow(grData_remeas)) {
  print(i)
  tmpPlot <- grData_remeas[i, "PLT_CN"]
  condsMatch <- subset(conds, PLT_CN == tmpPlot)
  balive <- condsMatch$BALIVE
  props <- condsMatch$CONDPROP_UNADJ
  wmean <- weighted.mean(balive, props, na.rm = T)
  if (length(balive) == 0) grData_remeas[i, "baLive"] <- NA
  else grData_remeas[i, "baLive"] <- wmean
}
grData_remeas[is.nan(grData_remeas$baLive), "baLive"] <- NA
grData_remeas <- subset(grData_remeas, !is.na(baLive))

# Make data spatial
grSpat <- SpatialPointsDataFrame(coords = cbind(grData_remeas$LON, grData_remeas$LAT), 
                                 data = grData_remeas, 
                                 proj4string = CRS("+proj=longlat +datum=NAD83"))

# Read and extract elevation
SRTM <- raster("D:/EvansLab/Final/Data/SRTM/SRTM_all.tif")
elev.extr <- raster::extract(SRTM, grSpat)
grData_remeas$elev <- elev.extr

# Read and extract climate

# Read in PRISM climate stacks
ppt <- stack("D:/EvansLab/Final/Data/PRISM/Historic/pptStack.tif")
tmp <- stack("D:/EvansLab/Final/Data/PRISM/Historic/tmpStack.tif")
vpd <- stack("D:/EvansLab/Final/Data/PRISM/Historic/vpdStack.tif")

# raster::extract PRISM data
ppt.extr <- raster::extract(ppt, grSpat)
tmp.extr <- raster::extract(tmp, grSpat)
vpd.extr <- raster::extract(vpd, grSpat)

# Add sensible column names for raster::extracted climate data
ppt.extr <- as.data.frame(ppt.extr)
tmp.extr <- as.data.frame(tmp.extr)
vpd.extr <- as.data.frame(vpd.extr)
climateDir <- "D:/EvansLab/Final/Data/PRISM/Historic/"
pptFiles <- list.files(path = climateDir, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
tmpFiles <- list.files(path = climateDir, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
vpdFiles <- list.files(path = climateDir, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
colNames <- lapply(strsplit(pptFiles, "4kmM._"), function (x) x[2])
colNames <- unlist(colNames)
colNames <- lapply(strsplit(colNames, "_"), function (x) x[1])
colNames <- unlist(colNames)
colnames(ppt.extr) <- paste0("ppt_", colNames)
colnames(tmp.extr) <- paste0("tmp_", colNames)
colnames(vpd.extr) <- paste0("vpd_", colNames[1:ncol(vpd.extr)])

# Remove partial 2016 data
ppt.extr <- ppt.extr[, 1:1452]
tmp.extr <- tmp.extr[, 1:1452]
vpd.extr <- vpd.extr[, 1:1452]

# Export climate data
write.csv(ppt.extr, "D:/EvansLab/Final/Data/Processed/Growth/ppt_extr.csv", row.names = F)
write.csv(tmp.extr, "D:/EvansLab/Final/Data/Processed/Growth/tmp_extr.csv", row.names = F)
write.csv(vpd.extr, "D:/EvansLab/Final/Data/Processed/Growth/vpd_extr.csv", row.names = F)

# Calculate seasonal climate for each year
for (i in 1896:2015) {
  print(i)
  ppt.extr[, paste0("PPT_c_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i-1, "10"),
                                                          paste0("ppt_", i-1, "11"), 
                                                          paste0("ppt_", i-1, "12"), 
                                                          paste0("ppt_", i, "01"),
                                                          paste0("ppt_", i, "02"),
                                                          paste0("ppt_", i, "03"))])
  tmp.extr[, paste0("T_c_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i-1, "10"),
                                                         paste0("tmp_", i-1, "11"), 
                                                         paste0("tmp_", i-1, "12"), 
                                                         paste0("tmp_", i, "01"),
                                                         paste0("tmp_", i, "02"),
                                                         paste0("tmp_", i, "03"))])
  vpd.extr[, paste0("VPD_c_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i-1, "10"),
                                                           paste0("vpd_", i-1, "11"), 
                                                           paste0("vpd_", i-1, "12"), 
                                                           paste0("vpd_", i, "01"),
                                                           paste0("vpd_", i, "02"),
                                                           paste0("vpd_", i, "03"))])
  
  ppt.extr[, paste0("PPT_w_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i, "05"),
                                                          paste0("ppt_", i, "06"),
                                                          paste0("ppt_", i, "07"))])
  tmp.extr[, paste0("T_w_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i, "05"),
                                                         paste0("tmp_", i, "06"),
                                                         paste0("tmp_", i, "07"))])
  vpd.extr[, paste0("VPD_w_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i, "05"),
                                                           paste0("vpd_", i, "06"),
                                                           paste0("vpd_", i, "07"))])
}

# Add climate variables to data frame
for (i in 1:nrow(grData_remeas)) {
  print(i)
  grData_remeas[i, "PPT_c"] <- rowMeans(ppt.extr[i, paste0("PPT_c_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
  grData_remeas[i, "VPD_c"] <- rowMeans(vpd.extr[i, paste0("VPD_c_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
  grData_remeas[i, "PPT_w"] <- rowMeans(ppt.extr[i, paste0("PPT_w_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
  grData_remeas[i, "VPD_w"] <- rowMeans(vpd.extr[i, paste0("VPD_w_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
}

# Create output data frame
output <- grData_remeas[, c("CN", "PREV_TRE_CN", "PLT_CN", "PREV_PLT_CN", "LAT", "LON", "DRYBIO_AG", "PREV_DRYBIO_AG", "DRYBIO_AG_DIFF", "MEASYEAR", "PREV_MEASYEAR",
                        "CENSUS_INTERVAL", "baLive", "PPT_c", "VPD_c", "PPT_w", "VPD_w", "elev")]
write.csv(output, "D:/EvansLab/Final/Data/Processed/Growth/GrowthData.csv", row.names = F)