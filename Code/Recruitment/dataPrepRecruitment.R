library(raster)
library(car)
library(rgdal)

##PLOTS##

# Read plots
plots_subset <- read.csv("D:/EvansLab/Final/Data/FIA/PLOT_SUBSET.csv", header = T, stringsAsFactors = F)

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

##PIED##

# Read in tree data
trees <- read.csv("D:/EvansLab/Final/Data/FIA/TREE_COMBINED.csv")

# Loop through trees and code PIED presence/absence
for (i in 1:nrow(out)) {
  plot <- out[i, "plot"]
  print(i / nrow(out) * 100)
  tmp <- subset(trees, PLT_CN == plot)
  if (nrow(tmp) == 0) out[i, "PIED"] <- NA
  else if (106 %in% unique(tmp$SPCD)) out[i, "PIED"] <- 1
  else out[i, "PIED"] <- 0
  print(out[i, "PIED"])
}

# Get rid of all plots with no TREE records
out <- out[!is.na(out$PIED),]

##DATE##

# Add field for previous measurement year
for (i in 1:nrow(out)) {
  tmp <- plots_subset[out[i, "plotPrev"] == plots_subset$CN, "MEASYEAR"]
  print(tmp)
  out[i, "measYearPrev"] <- tmp
}
out$measInterval <- out$measYear - out$measYearPrev

# Create spatial data frame
outSp <- SpatialPointsDataFrame(coords = cbind(out$lon, out$lat), 
                                data = data.frame(plot = out$plot), 
                                proj4string = CRS("+proj=longlat +datum=NAD83"))

##ELEVATION##

# raster::extract elevation
elev <- raster("D:/EvansLab/Final/Data/SRTM/SRTM_all.tif")
out$elev <- raster::extract(elev, outSp)

##PRISM NORMALS##

# Read and raster::extract PRISM normals (PPT, T, VPD_max, VPD_min)
prism_files <- list.files(path = "D:/EvansLab/Final/Data/PRISM/Normals/", pattern = "bil$", full.names = T)
prism <- stack(prism_files)
outSp <- raster::extract(prism, outSp, sp = T)

# Aggregate PRISM normals (cool season: 10-03)
outSp$PPT_cool <- 
  outSp$PRISM_ppt_30yr_normal_800mM2_10_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_11_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_12_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_01_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_02_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_03_bil
outSp$T_cool <- (
  outSp$PRISM_tmean_30yr_normal_800mM2_10_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_11_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_12_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_01_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_02_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_03_bil) / 6
outSp$VPDmin_cool <- (
  outSp$PRISM_vpdmin_30yr_normal_800mM2_10_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_11_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_12_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_01_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_02_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_03_bil) / 6
outSp$VPDmax_cool <- (
  outSp$PRISM_vpdmax_30yr_normal_800mM2_10_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_11_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_12_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_01_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_02_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_03_bil) / 6
outSp$VPD_cool <- (outSp$VPDmax_cool + outSp$VPDmin_cool) / 2

# Aggregate PRISM normals (warm season: 05-07)

outSp$PPT_warm <- 
  outSp$PRISM_ppt_30yr_normal_800mM2_05_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_06_bil + 
  outSp$PRISM_ppt_30yr_normal_800mM2_07_bil
outSp$T_warm <- (
  outSp$PRISM_tmean_30yr_normal_800mM2_05_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_06_bil + 
  outSp$PRISM_tmean_30yr_normal_800mM2_07_bil) / 3
outSp$VPDmin_warm <- (
  outSp$PRISM_vpdmin_30yr_normal_800mM2_05_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_06_bil + 
  outSp$PRISM_vpdmin_30yr_normal_800mM2_07_bil) / 3
outSp$VPDmax_warm <- (
  outSp$PRISM_vpdmax_30yr_normal_800mM2_05_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_06_bil + 
  outSp$PRISM_vpdmax_30yr_normal_800mM2_07_bil) / 3
outSp$VPD_warm <- (outSp$VPDmax_warm + outSp$VPDmin_warm) / 2

# Copy normals to original data frame
out$PPT_w_normal <- outSp$PPT_warm
out$PPT_c_normal <- outSp$PPT_cool
out$T_w_normal <-   outSp$T_warm
out$T_c_normal <-   outSp$T_cool
out$VPD_w_normal <- outSp$VPD_warm
out$VPD_c_normal <- outSp$VPD_cool

##RECRUITS##

# Loop through plots
out$recruits <- 0
recruits <- subset(trees, RECONCILECD %in% c(1,2) & SPCD == 106)
for (i in unique(recruits$PLT_CN)) {
  tmp <- subset(recruits, PLT_CN == i)
  print(nrow(tmp))
  out[out$plot == i, "recruits"] <- nrow(tmp)
}

##COMPETITION##

# Read in conditions
conds <- read.csv("D:/EvansLab/Final/Data/FIA/COND_COMBINED.csv")
conds <- subset(conds, COND_STATUS_CD == 1)

for (i in 1:nrow(out)) {
  tmpPlot <- out[i, "plot"]
  condsMatch <- subset(conds, PLT_CN == as.character(tmpPlot))
  balive <- condsMatch$BALIVE
  print(i)
  props <- condsMatch$CONDPROP_UNADJ
  wmean <- weighted.mean(balive, props, na.rm = T)
  if (length(balive) == 0) out[i, "baLive"] <- NA
  else out[i, "baLive"] <- wmean
}

out[is.nan(out$baLive), "baLive"] <- NA
out <- out[!is.na(out$baLive),]

##PRISM conditions##

# Create spatial points layer
outPoints <- SpatialPoints(coords = cbind(out$lon, out$lat), 
                           proj4string = CRS("+proj=longlat +datum=NAD83"))

# Read in PRISM climate stacks
ppt <- stack("D:/EvansLab/Final/Data/PRISM/Historic/pptStack.tif")
tmp <- stack("D:/EvansLab/Final/Data/PRISM/Historic/tmpStack.tif")
vpd <- stack("D:/EvansLab/Final/Data/PRISM/Historic/vpdStack.tif")

# raster::extract PRISM data
ppt.extr <- raster::extract(ppt, outPoints)
tmp.extr <- raster::extract(tmp, outPoints)
vpd.extr <- raster::extract(vpd, outPoints)

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
write.csv(ppt.extr, "D:/EvansLab/Final/Data/Processed/Recruitment/ppt_extr.csv", row.names = F)
write.csv(tmp.extr, "D:/EvansLab/Final/Data/Processed/Recruitment/tmp_extr.csv", row.names = F)
write.csv(vpd.extr, "D:/EvansLab/Final/Data/Processed/Recruitment/vpd_extr.csv", row.names = F)

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

# Calculate seasonal averages
lags <- c(15, 20, 25)
for (j in lags) {
  print(j)
  lagLength <- j
  for (i in 1:nrow(ppt.extr)) {
    print(i)
    out[i, paste0("PPT_c_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_c_", (out[i, "measYear"]-lagLength):(out[i, "measYear"]-1))])
    out[i, paste0("T_c_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_c_", (out[i, "measYear"]-lagLength):(out[i, "measYear"]-1))])
    out[i, paste0("VPD_c_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_c_", (out[i, "measYear"]-lagLength):(out[i, "measYear"]-1))])
    out[i, paste0("PPT_w_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_w_", (out[i, "measYear"]-lagLength):(out[i, "measYear"]-1))])
    out[i, paste0("T_w_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_w_", (out[i, "measYear"]-lagLength):(out[i, "measYear"]-1))])
    out[i, paste0("VPD_w_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_w_", (out[i, "measYear"]-lagLength):(out[i, "measYear"]-1))])
  }
}

##COMPETITION##

# Calculate AGB (interspecific, intraspecific, all)
for (i in 1:nrow(out)) {
  print(i)
  tmpPlot <- out[i, "plot"]
  tmpTrees <- trees[trees$PLT_CN == tmpPlot,]
  tmpIntra <- sum(subset(tmpTrees, SPCD == 106)$DRYBIO_AG, na.rm = T)
  tmpInter <- sum(subset(tmpTrees, SPCD != 106)$DRYBIO_AG, na.rm = T)
  tmpSum <- tmpIntra + tmpInter
  out[i, "AGB_intra"] <- tmpIntra
  out[i, "AGB_inter"] <- tmpInter
  out[i, "AGB_all"] <- tmpSum
}

##EXPORT##
out <- out[complete.cases(out),]
write.csv(out, "D:/EvansLab/Final/Data/Processed/Recruitment/RecruitmentData.csv", row.names = F)