library(raster)
library(car)
library(rgdal)

##PLOTS##
data.path <- "./FIAdata/"

# Read plots
plots <- read.csv(paste(data.path,"PLOT_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
# Read in conditions
conds <- read.csv(paste(data.path,"COND_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
conds <- subset(conds, COND_STATUS_CD == 1) # "accessible forest land" by FIA classification

# Create data frame that will be used for zero-inflated Poisson (count) regression
# rData needs to be split by CONDID (I didn't implement this)
rData <- data.frame(plot = as.character(plots$CN),
                  lat = plots$LAT, 
                  lon = plots$LON, 
                  elev = plots$ELEV,
                  state = plots$STATECD,
                  county = plots$COUNTYCD,
                  plotID = plots$PLOT,
                  PApied = integer(length(plots$CN)), # will contain P/A information
                  measYear = plots$MEASYEAR, 
                  plotPrev = as.character(plots$PREV_PLT_CN))
rData$plot <- as.character(rData$plot) # strangely enough, this line is necessary for the PA loop to work
rData$plotPrev <- as.character(rData$plotPrev)

# Get rid of plots that were not remeasured
rData <- subset(rData, !is.na(rData$plotPrev))
# look up time 1 condition ID for each plot...but this isn't informative until plots are split by CONDID
rData$CONDID <- conds$CONDID[match(rData$plotPrev, conds$PLT_CN)]


## was PIED present or absent in each FIA plot? ##
# Read in tree data
trees <- read.csv(paste(data.path,"TREE_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
# make list of those tree records that fall in (time 2) plots that are found in "rData"
PA.trees <- subset(trees, PLT_CN %in% unique(rData$plot))
PA.list <- split(x=PA.trees, f=PA.trees$PLT_CN) 

# Loop through trees and code PIED presence/absence at time 2
# note that this includes dead trees, not just live trees
# needs to be split by CONDID within plots, see above...not implemented
for (i in 1:nrow(rData)) {
  plot <- rData[i, "plot"] # use concatenation of PLT and CONDID?
  print(i / nrow(rData) * 100)
  tmp <- PA.list[[plot]]
  if (is.null(tmp)) rData[i, "PApied"] <- NA # no tree records for this plot CN
  else if (106 %in% unique(tmp$SPCD)) rData[i, "PApied"] <- 1
  else rData[i, "PApied"] <- 0
  print(rData[i, "PApied"])
}

# Get rid of all plots with no TREE records, no CONDID
rData <- rData[!is.na(rData$PApied),]
#rData <- rData[!is.na(rData$CONDID),] # not necessary until the data are split by CONDID

#save(file = "regen_out.RData", list =c("out"))
# load saved RData
#load(file = "../../Processed/Recruitment/regen_out.RData")

# Add fields for previous measurement year, census interval
rData$PREV_MEASYEAR <- plots$MEASYEAR[match(rData$plotPrev, plots$CN)]
rData$CENSUS_INTERVAL <- rData$measYear - rData$PREV_MEASYEAR


## count PIED RECRUITS ##

# Loop through plots
# more generous definition of recruits
rData$recruits12 <- 0
recruits12 <- subset(trees, RECONCILECD %in% c(1,2) & SPCD == 106)
for (i in unique(recruits12$PLT_CN)) {
  tmp <- subset(recruits12, PLT_CN == i) # add CONDID to this subsetting
  print(nrow(tmp))
  rData[rData$plot == i, "recruits12"] <- nrow(tmp)
}

# more constrained definition of recruits
# this is what I ended up using in ZIP models
rData$recruits1 <- 0
recruits1 <- subset(trees, RECONCILECD == 1 & SPCD == 106)
for (i in unique(recruits1$PLT_CN)) {
  tmp <- subset(recruits1, PLT_CN == i) # add CONDID to this subsetting
  print(nrow(tmp))
  rData[rData$plot == i, "recruits1"] <- nrow(tmp)
}


# there is precious little difference between these two:
rData$diff <- rData$recruits12-rData$recruits1
sum(rData$diff) # returns the value 3


##create COMPETITION covariate data##

plot.conds <- subset(conds, PLT_CN %in% unique(rData$plotPrev))
for (i in 1:nrow(rData)) {
  tmpPlot <- rData[i, "plotPrev"] # time 1 plot control number
  condsMatch <- subset(plot.conds, PLT_CN == as.character(tmpPlot))
  balive <- condsMatch$BALIVE
  print(i)
  props <- condsMatch$CONDPROP_UNADJ
  wmean <- weighted.mean(balive, props, na.rm = T)
  if (length(balive) == 0) rData[i, "BALIVE"] <- NA
  else rData[i, "BALIVE"] <- wmean
}

rData[is.nan(rData$BALIVE), "BALIVE"] <- NA
rData <- subset(rData, !is.na(BALIVE))


##create OFFSET (SEED SOURCE) data##
## ...and additional competition/facilitation offset variables that might affect recruitment rate

# get rid of trees that are not in the time 1 plots in "rData" (just to make the search smaller)
plot.trees <- subset(trees, PLT_CN %in% unique(rData$plotPrev))
# create a list (indexed by plot) of trees (this also increases efficiency)
trees.list <- split(x=plot.trees, f=plot.trees$PLT_CN) #6723
missing <- setdiff(x=rData$plotPrev, y=trees$PLT_CN)
head(rData[rData$plotPrev %in% missing,])
# next line gets rid of those plots that had no trees at time 1 census
rData <- rData[!(rData$plotPrev %in% missing),]
# but trees could have colonized these 18 plots between time 1 and time 2

for (i in 1:nrow(rData)) {
  plot <- rData[i, "plotPrev"] # use concatenation of PLT and CONDID?
  print(i / nrow(rData) * 100)
#  tmp.PIED <- subset(trees.list[[plot]], SPCD == 106 & STATUSCD == 1)
#  tmp.notPIED <- subset(trees.list[[plot]], SPCD != 106 & STATUSCD == 1)
#  tmp.juniper <- subset(trees.list[[plot]], SPCD > 56 & SPCD < 70 & STATUSCD == 1)
#  tmp.PIPO <- subset(trees.list[[plot]], SPCD == 122 & STATUSCD == 1)
  plot.trees <- trees.list[[plot]]
  tmp.PIED <- plot.trees[plot.trees$SPCD == 106 & plot.trees$STATUSCD == 1, ]
  tmp.notPIED <- plot.trees[plot.trees$SPCD != 106 & plot.trees$STATUSCD == 1, ]
  tmp.juniper <- plot.trees[plot.trees$SPCD > 56 & plot.trees$SPCD < 70 & plot.trees$STATUSCD == 1, ]
  tmp.PIPO <- plot.trees[plot.trees$SPCD == 122 & plot.trees$STATUSCD == 1, ]
  rData$PIEDadults1[i] <- nrow(tmp.PIED) # all trees >1" DRC
  rData$otheradults1[i] <- nrow(tmp.notPIED) # all trees >1" DRC
  rData$PIEDadults4[i] <- nrow(subset(tmp.PIED, DIA > 4.0)) # number of trees >4" DRC
  rData$otheradults4[i] <- nrow(subset(tmp.notPIED, DIA > 4.0)) 
  rData$PIEDadults8[i] <- nrow(subset(tmp.PIED, DIA > 8.0)) # number of trees >8" DRC
  rData$otheradults8[i] <- nrow(subset(tmp.notPIED, DIA > 8.0))                               
  rData$cumDIA.PIED[i] <- sum(tmp.PIED$DIA)
  rData$cumDIA.others[i] <- sum(tmp.notPIED$DIA)
  rData$BA.PIED[i] <- sum(tmp.PIED$DIA^2*0.005454, na.rm =T) # basal area in units of sq feet
  # add BA.PIED of trees DIA > 4.0, i.e., reproductive
  rData$BA.notPIED[i] <- sum(tmp.notPIED$DIA^2*0.005454, na.rm =T)
  rData$BA.juniper[i] <- sum(tmp.juniper$DIA^2*0.005454, na.rm =T)
  rData$BA.PIPO[i] <- sum(tmp.PIPO$DIA^2*0.005454, na.rm =T)
  rData$AGB_intra[i] <- sum(tmp.PIED$DRYBIO_AG, na.rm = T)
  rData$AGB_inter[i] <- sum(tmp.notPIED$DRYBIO_AG, na.rm = T)
}

table(rData[, c("recruits1", "PApied")]) # PApied is presence/absence at time 2, including dead trees
table(rData[, c("PIEDadults1", "recruits1")])
table(rData[, c("PIEDadults4", "recruits1")])
# table(rData[, c("BA.juniper", "recruits1")])

##create PRISM covariate data##
### (most of the code in this entire script has to do with creating climate covariate data)

# Read and extract climate
# THE FOLLOWING LINES CAN BE SKIPPED if you've already generated the "ppt_extr.csv" etc. files
# Create spatial points layer
Points <- SpatialPoints(coords = cbind(rData$lon, rData$lat), 
                           proj4string = CRS("+proj=longlat +datum=NAD83"))

# Read in PRISM climate stacks
#clim.path <-  "F:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/"
clim.path <-  "./ClimateData/"
ppt <- stack(paste(clim.path,"pptStack.tif",sep=''))
tmp <- stack(paste(clim.path,"tmpStack.tif",sep=''))
vpd <- stack(paste(clim.path,"vpdStack.tif",sep=''))

# raster::extract PRISM data
ppt.extr <- raster::extract(ppt, Points) # takes ~9 minutes each
tmp.extr <- raster::extract(tmp, Points)
vpd.extr <- raster::extract(vpd, Points)

# Remove data after Oct, 2016 (because of different CRS Nov, 2016 vpdmax .bil)
# note that the work-around for this problem is to assign the CRS of another layer to Nov and Dec of 2016
# crs(vpdNov2016_raster) <- crs(vpdOct2016_raster)
#ppt.extr <- ppt.extr[, 1:430] 
#tmp.extr <- tmp.extr[, 1:430]
#vpd.extr <- vpd.extr[, 1:430]

# Add sensible column names for raster::extracted climate data
ppt.extr <- as.data.frame(ppt.extr)
tmp.extr <- as.data.frame(tmp.extr)
vpd.extr <- as.data.frame(vpd.extr)
#PRISM.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/"
PRISM.path <- "./ClimateData/PRISM"
pptFiles <- list.files(path = PRISM.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
#pptFiles <- pptFiles[1:430] # (hack to deal with CRS incompatibility, vpd .bil file Nov, 2016)
#tmpFiles <- list.files(path = climateDir, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
#vpdFiles <- list.files(path = climateDir, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
colNames <- lapply(strsplit(pptFiles, "4kmM._"), function (x) x[2])
colNames <- unlist(colNames)
colNames <- lapply(strsplit(colNames, "_"), function (x) x[1])
colNames <- unlist(colNames)
colnames(ppt.extr) <- paste0("ppt_", colNames)
colnames(tmp.extr) <- paste0("tmp_", colNames)
colnames(vpd.extr) <- paste0("vpd_", colNames)

# Export climate data
#recruit.path <- "C:/Users/mekevans/Documents/Cdrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PIED_IPM/MEKEvans/Processed/Recruitment/"
recruit.path <- "./Processed/Recruitment/"
write.csv(ppt.extr, paste0(recruit.path, "ppt_extr.csv"), row.names = F)
write.csv(tmp.extr, paste0(recruit.path, "tmp_extr.csv"), row.names = F)
write.csv(vpd.extr, paste0(recruit.path, "vpd_extr.csv"), row.names = F)

### CONTINUE HERE IF YOU ALREADY HAVE GENERATED THE ABOVE csv files
#recruit.path <- "C:/Users/mekevans/Documents/Cdrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PIED_IPM/MEKEvans/Processed/Recruitment/"
recruit.path <- "./Processed/Recruitment/"
ppt.extr <- read.csv(paste(recruit.path,"ppt_extr.csv",sep=''), header = T)
tmp.extr <- read.csv(paste(recruit.path,"tmp_extr.csv",sep=''), header = T)
vpd.extr <- read.csv(paste(recruit.path,"vpd_extr.csv",sep=''), header = T)


# Calculate seasonal climate for each year, 1982 through 2016
for (i in 1982:2017) {
  print(i)
  # cool season = pNov - Mar
  ppt.extr[, paste0("PPT_c_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i-1, "11"), 
                                                          paste0("ppt_", i-1, "12"), 
                                                          paste0("ppt_", i, "01"),
                                                          paste0("ppt_", i, "02"),
                                                          paste0("ppt_", i, "03"))])
  tmp.extr[, paste0("T_c_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i-1, "11"), 
                                                         paste0("tmp_", i-1, "12"), 
                                                         paste0("tmp_", i, "01"),
                                                         paste0("tmp_", i, "02"),
                                                         paste0("tmp_", i, "03"))])
  tmp.extr[, paste0("Tex_c_", i)] <- max(tmp.extr[, c(paste0("tmp_", i-1, "11"), 
                                                      paste0("tmp_", i-1, "12"), 
                                                      paste0("tmp_", i, "01"),
                                                      paste0("tmp_", i, "02"),
                                                      paste0("tmp_", i, "03"))])
  vpd.extr[, paste0("VPD_c_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i-1, "11"), 
                                                           paste0("vpd_", i-1, "12"), 
                                                           paste0("vpd_", i, "01"),
                                                           paste0("vpd_", i, "02"),
                                                           paste0("vpd_", i, "03"))])
  vpd.extr[, paste0("VPDex_c_", i)] <- max(vpd.extr[, c(paste0("vpd_", i-1, "11"), 
                                                        paste0("vpd_", i-1, "12"), 
                                                        paste0("vpd_", i, "01"),
                                                        paste0("vpd_", i, "02"),
                                                        paste0("vpd_", i, "03"))])
  # previous fall = pSep - pOct
  ppt.extr[, paste0("PPT_pf_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i-1, "09"),
                                                           paste0("ppt_", i-1, "10"))])
  tmp.extr[, paste0("T_pf_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i-1, "09"),
                                                          paste0("tmp_", i-1, "10"))])
  tmp.extr[, paste0("Tex_pf_", i)] <- max(tmp.extr[, c(paste0("tmp_", i-1, "09"),
                                                       paste0("tmp_", i-1, "10"))])
  vpd.extr[, paste0("VPD_pf_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i-1, "09"),
                                                            paste0("vpd_", i-1, "10"))])
  vpd.extr[, paste0("VPDex_pf_", i)] <- max(vpd.extr[, c(paste0("vpd_", i-1, "09"),
                                                         paste0("vpd_", i-1, "10"))])
  # foresummer = Apr - Jun
  ppt.extr[, paste0("PPT_fs_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i, "04"),
                                                           paste0("ppt_", i, "05"),
                                                           paste0("ppt_", i, "06"))])
  tmp.extr[, paste0("T_fs_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i, "04"),
                                                          paste0("tmp_", i, "05"),
                                                          paste0("tmp_", i, "06"))])
  tmp.extr[, paste0("Tex_fs_", i)] <- max(tmp.extr[, c(paste0("tmp_", i, "04"),
                                                       paste0("tmp_", i, "05"),
                                                       paste0("tmp_", i, "06"))])
  vpd.extr[, paste0("VPD_fs_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i, "04"),
                                                            paste0("vpd_", i, "05"),
                                                            paste0("vpd_", i, "06"))])
  vpd.extr[, paste0("VPDex_fs_", i)] <- max(vpd.extr[, c(paste0("vpd_", i, "04"),
                                                         paste0("vpd_", i, "05"),
                                                         paste0("vpd_", i, "06"))])
  
  # warm dry months = pSep - pOct + Apr - Jun
  ppt.extr[, paste0("PPT_wd_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i-1, "09"),
                                                           paste0("ppt_", i-1, "10"),
                                                           paste0("ppt_", i, "04"),
                                                           paste0("ppt_", i, "05"),
                                                           paste0("ppt_", i, "06"))])
  tmp.extr[, paste0("T_wd_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i-1, "09"),
                                                          paste0("tmp_", i-1, "10"),
                                                          paste0("tmp_", i, "04"),
                                                          paste0("tmp_", i, "05"),
                                                          paste0("tmp_", i, "06"))])
  tmp.extr[, paste0("Tex_wd_", i)] <- max(tmp.extr[, c(paste0("tmp_", i-1, "09"),
                                                       paste0("tmp_", i-1, "10"),
                                                       paste0("tmp_", i, "04"),
                                                       paste0("tmp_", i, "05"),
                                                       paste0("tmp_", i, "06"))])
  vpd.extr[, paste0("VPD_wd_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i-1, "09"),
                                                            paste0("vpd_", i-1, "10"),
                                                            paste0("vpd_", i, "04"),
                                                            paste0("vpd_", i, "05"),
                                                            paste0("vpd_", i, "06"))])
  vpd.extr[, paste0("VPDex_wd_", i)] <- max(vpd.extr[, c(paste0("vpd_", i-1, "09"),
                                                         paste0("vpd_", i-1, "10"),
                                                         paste0("vpd_", i, "04"),
                                                         paste0("vpd_", i, "05"),
                                                         paste0("vpd_", i, "06"))])
  
  # monsoon = Jul-Aug
  ppt.extr[, paste0("PPT_m_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i, "07"),
                                                          paste0("ppt_", i, "08"))])
  tmp.extr[, paste0("T_m_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i, "07"),
                                                         paste0("tmp_", i, "08"))])
  tmp.extr[, paste0("Tex_m_", i)] <- max(tmp.extr[, c(paste0("tmp_", i, "07"),
                                                      paste0("tmp_", i, "08"))])
  vpd.extr[, paste0("VPD_m_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i, "07"),
                                                           paste0("vpd_", i, "08"))])
  vpd.extr[, paste0("VPDex_m_", i)] <- max(vpd.extr[, c(paste0("vpd_", i, "07"),
                                                        paste0("vpd_", i, "08"))])
  # water year = pSept - Aug
  ppt.extr[, paste0("PPT_yr_", i)] <- rowSums(ppt.extr[, c(paste0("ppt_", i-1, "09"),
                                                           paste0("ppt_", i-1, "10"),
                                                           paste0("ppt_", i-1, "11"), 
                                                           paste0("ppt_", i-1, "12"), 
                                                           paste0("ppt_", i, "01"),
                                                           paste0("ppt_", i, "02"),
                                                           paste0("ppt_", i, "03"),
                                                           paste0("ppt_", i, "04"),
                                                           paste0("ppt_", i, "05"),
                                                           paste0("ppt_", i, "06"),
                                                           paste0("ppt_", i, "07"),
                                                           paste0("ppt_", i, "08"))])
  tmp.extr[, paste0("T_yr_", i)] <- rowMeans(tmp.extr[, c(paste0("tmp_", i-1, "09"),
                                                          paste0("tmp_", i-1, "10"),
                                                          paste0("tmp_", i-1, "11"), 
                                                          paste0("tmp_", i-1, "12"), 
                                                          paste0("tmp_", i, "01"),
                                                          paste0("tmp_", i, "02"),
                                                          paste0("tmp_", i, "03"),
                                                          paste0("tmp_", i, "04"),
                                                          paste0("tmp_", i, "05"),
                                                          paste0("tmp_", i, "06"),
                                                          paste0("tmp_", i, "07"),
                                                          paste0("tmp_", i, "08"))])
  vpd.extr[, paste0("VPD_yr_", i)] <- rowMeans(vpd.extr[, c(paste0("vpd_", i-1, "09"),
                                                            paste0("vpd_", i-1, "10"),
                                                            paste0("vpd_", i-1, "11"), 
                                                            paste0("vpd_", i-1, "12"), 
                                                            paste0("vpd_", i, "01"),
                                                            paste0("vpd_", i, "02"),
                                                            paste0("vpd_", i, "03"),
                                                            paste0("vpd_", i, "04"),
                                                            paste0("vpd_", i, "05"),
                                                            paste0("vpd_", i, "06"),
                                                            paste0("vpd_", i, "07"),
                                                            paste0("vpd_", i, "08"))])
}


##PRISM NORMALS##

# Read and raster::extract PRISM normals (PPT, T, VPD_max, VPD_min)
# import PRISM normals
PRISM.norm.path <-  "./ClimateData/"
#PRISM.norm.path <-  "F:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/"
PPT.norm <- stack(paste(clim.path,"pptNormals.tif",sep=''))
TMP.norm <- stack(paste(clim.path,"tmpNormals.tif",sep=''))
VPD.norm <- stack(paste(clim.path,"vpdNormals.tif",sep=''))

# raster::extract PRISM normals 1981-2010
ppt.norm.extr <- raster::extract(PPT.norm, Points)
tmp.norm.extr <- raster::extract(TMP.norm, Points)
vpd.norm.extr <- raster::extract(VPD.norm, Points)

ppt.norm.extr <- as.data.frame(ppt.norm.extr) # not sure this is necessary?
tmp.norm.extr <- as.data.frame(tmp.norm.extr)
vpd.norm.extr <- as.data.frame(vpd.norm.extr)

# reasonable column names (not actually necessary)
colnames(ppt.norm.extr) <- paste0("ppt_", 1:12) 
colnames(tmp.norm.extr) <- paste0("tmp_", 1:12)
colnames(vpd.norm.extr) <- paste0("vpd_", 1:12)

# make seasonal normals and add to growth data frame
# cool season = Nov-Mar
rData$PPT_c_norm <- rowSums(ppt.norm.extr[, c(1:3, 11:12)])
rData$T_c_norm <- rowMeans(tmp.norm.extr[, c(1:3, 11:12)])
rData$VPD_c_norm <- rowMeans(vpd.norm.extr[, c(1:3, 11:12)])
# previous fall = pSept-pOct
rData$PPT_pf_norm <- rowSums(ppt.norm.extr[, c(9:10)])
rData$T_pf_norm <- rowMeans(tmp.norm.extr[, c(9:10)])
rData$VPD_pf_norm <- rowMeans(vpd.norm.extr[, c(9:10)])
# foresummer = Apr-Jun
rData$PPT_fs_norm <- rowSums(ppt.norm.extr[, c(4:6)])
rData$T_fs_norm <- rowMeans(tmp.norm.extr[, c(4:6)])
rData$VPD_fs_norm <- rowMeans(vpd.norm.extr[, c(4:6)])
# warm, dry months = Apr-Jun + Sept-Oct
rData$PPT_wd_norm <- rowSums(ppt.norm.extr[, c(4:6, 9:10)])
rData$T_wd_norm <- rowMeans(tmp.norm.extr[, c(4:6, 9:10)])
rData$VPD_wd_norm <- rowMeans(vpd.norm.extr[, c(4:6, 9:10)])
# monsoon = Jul-Aug
rData$PPT_m_norm <- rowSums(ppt.norm.extr[, c(7:8)])
rData$T_m_norm <- rowMeans(tmp.norm.extr[, c(7:8)])
rData$VPD_m_norm <- rowMeans(vpd.norm.extr[, c(7:8)])
# water year
rData$PPT_yr_norm <- rowSums(ppt.norm.extr[, c(1:12)])
rData$T_yr_norm <- rowMeans(tmp.norm.extr[, c(1:12)])
rData$VPD_yr_norm <- rowMeans(vpd.norm.extr[, c(1:12)])

# Calculate seasonal averages
lags <- c(15, 20, 25) # if I try to add 30-yr lag, bump up against 1981 limit of PRISM data
for (j in lags) {
  print(j)
  lagLength <- j
  for (i in 1:nrow(ppt.extr)) {
    print(i)
# cool season (pNov-Mar)
    rData[i, paste0("PPT_c_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_c_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("PPTex_c_window_", lagLength)] <- min(ppt.extr[i, paste0("PPT_c_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("T_c_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_c_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("Tex_c_window_", lagLength)] <- max(tmp.extr[i, paste0("T_c_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])    
    rData[i, paste0("VPD_c_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_c_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPDex_c_window_", lagLength)] <- max(vpd.extr[i, paste0("VPD_c_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])    
# monsoon (Jul-Aug)    
    rData[i, paste0("PPT_m_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_m_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("PPTex_m_window_", lagLength)] <- min(ppt.extr[i, paste0("PPT_m_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("T_m_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_m_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("Tex_m_window_", lagLength)] <- max(tmp.extr[i, paste0("T_m_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPD_m_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_m_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPDex_m_window_", lagLength)] <- max(vpd.extr[i, paste0("VPD_m_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])    
# previous fall (pSep-pOct)    
    rData[i, paste0("PPT_pf_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_pf_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("PPTex_pf_window_", lagLength)] <- min(ppt.extr[i, paste0("PPT_pf_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("T_pf_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_pf_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("Tex_pf_window_", lagLength)] <- max(tmp.extr[i, paste0("T_pf_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPD_pf_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_pf_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPDex_pf_window_", lagLength)] <- max(vpd.extr[i, paste0("VPD_pf_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])    
# foresummer (Apr-Jun)
    rData[i, paste0("PPT_fs_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_fs_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("PPTex_fs_window_", lagLength)] <- min(ppt.extr[i, paste0("PPT_fs_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("T_fs_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_fs_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("Tex_fs_window_", lagLength)] <- max(tmp.extr[i, paste0("T_fs_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPD_fs_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_fs_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPDex_fs_window_", lagLength)] <- max(vpd.extr[i, paste0("VPD_fs_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])    
# warm dry months (pSep-pOct + Apr-Jun)    
    rData[i, paste0("PPT_wd_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_wd_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("PPTex_wd_window_", lagLength)] <- min(ppt.extr[i, paste0("PPT_wd_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("T_wd_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_wd_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("Tex_wd_window_", lagLength)] <- max(tmp.extr[i, paste0("T_wd_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPD_wd_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_wd_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPDex_wd_window_", lagLength)] <- max(vpd.extr[i, paste0("VPD_wd_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])    
# water year  
    rData[i, paste0("PPT_yr_window_", lagLength)] <- rowMeans(ppt.extr[i, paste0("PPT_yr_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("PPTex_yr_window_", lagLength)] <- min(ppt.extr[i, paste0("PPT_yr_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("T_yr_window_", lagLength)] <- rowMeans(tmp.extr[i, paste0("T_yr_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("Tex_yr_window_", lagLength)] <- max(tmp.extr[i, paste0("T_yr_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPD_yr_window_", lagLength)] <- rowMeans(vpd.extr[i, paste0("VPD_yr_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
    rData[i, paste0("VPDex_yr_window_", lagLength)] <- max(vpd.extr[i, paste0("VPD_yr_", (rData[i, "measYear"]-lagLength):(rData[i, "measYear"]-1))])
  }
}

# very little reason to expect that the climate variables above,
# specific to each plot's census interval, will differ from the 1981-2010 climate normals
# especially with the lags added

# create and add anomalies to growth data frame
# for just one lag...others would have to be calculated
# note that precip is lognormal (R-skewed), but precip anomalies are left-skewed (heavy tail in the negative)...because mean of a lognormal is right of center
# monsoon precip anomalies are centered on zero (no change)
# cool season, foresummer, and water year anomalies are left-skewed (heavy tail in negative values)
# temperature anomalies are positive (no surprise)
# cool season (pNov-Mar)
rData$PPT_c_anom <- rData$PPT_c_window_20 - rData$PPT_c_norm # negative
rData$T_c_anom <- rData$T_c_window_20 - rData$T_c_norm # positive
rData$VPD_c_anom <- rData$VPD_c_window_20 - rData$VPD_c_norm # positive
rData$PPTex_c_anom <- rData$PPTex_c_window_20 - rData$PPT_c_norm  # PPTex_c_window_20 does not exist
rData$Tex_c_anom <- rData$Tex_c_window_20 - rData$T_c_norm
rData$VPDex_c_anom <- rData$VPDex_c_window_20 - rData$VPD_c_norm

# previous fall (pSep-pOct)
rData$PPT_pf_anom <- rData$PPT_pf_window_20 - rData$PPT_pf_norm # no change
rData$T_pf_anom <- rData$T_pf_window_20 - rData$T_pf_norm # positive
rData$VPD_pf_anom <- rData$VPD_pf_window_20 - rData$VPD_pf_norm # positive
rData$PPTex_pf_anom <- rData$PPTex_pf_window_20 - rData$PPT_pf_norm # negative
rData$Tex_pf_anom <- rData$Tex_pf_window_20 - rData$T_pf_norm # positive
rData$VPDex_pf_anom <- rData$VPDex_pf_window_20 - rData$VPD_pf_norm # positive
# foresummer (Apr-Jun)
rData$PPT_fs_anom <- rData$PPT_fs_window_20 - rData$PPT_fs_norm # negative!
rData$T_fs_anom <- rData$T_fs_window_20 - rData$T_fs_norm # positive
rData$VPD_fs_anom <- rData$VPD_fs_window_20 - rData$VPD_fs_norm # positive
rData$PPTex_fs_anom <- rData$PPTex_fs_window_20 - rData$PPT_fs_norm # centered on zero
rData$Tex_fs_anom <- rData$Tex_fs_window_20 - rData$T_fs_norm # positive
rData$VPDex_fs_anom <- rData$VPDex_fs_window_20 - rData$VPD_fs_norm # positive
# warm dry months (pSep-pOct + Apr-Jun)
rData$PPT_wd_anom <- rData$PPT_wd_window_20 - rData$PPT_wd_norm # negative
rData$T_wd_anom <- rData$T_wd_window_20 - rData$T_wd_norm
rData$VPD_wd_anom <- rData$VPD_wd_window_20 - rData$VPD_wd_norm
rData$PPTex_wd_anom <- rData$PPTex_wd_window_20 - rData$PPT_wd_norm # positive!
rData$Tex_wd_anom <- rData$Tex_wd_window_20 - rData$T_wd_norm
rData$VPDex_wd_anom <- rData$VPDex_wd_window_20 - rData$VPD_wd_norm
# monsoon (Jul-Aug)
rData$PPT_m_anom <- rData$PPT_m_window_20 - rData$PPT_m_norm # negative
rData$T_m_anom <- rData$T_m_window_20 - rData$T_m_norm
rData$VPD_m_anom <- rData$VPD_m_window_20 - rData$VPD_m_norm
rData$PPTex_m_anom <- rData$PPTex_m_window_20 - rData$PPT_m_norm # centered on zero
rData$Tex_m_anom <- rData$Tex_m_window_20 - rData$T_m_norm
rData$VPDex_m_anom <- rData$VPDex_m_window_20 - rData$VPD_m_norm
# water year (pSept - Aug)
rData$PPT_yr_anom <- rData$PPT_yr_window_20 - rData$PPT_yr_norm # negative
rData$T_yr_anom <- rData$T_yr_window_20 - rData$T_yr_norm
rData$VPD_yr_anom <- rData$VPD_yr_window_20 - rData$VPD_yr_norm
rData$PPTex_yr_anom <- rData$PPTex_yr_window_20 - rData$PPT_yr_norm
rData$Tex_yr_anom <- rData$Tex_yr_window_20 - rData$T_yr_norm  # Tex_yr_window_20 does not exist
rData$VPDex_yr_anom <- rData$VPDex_yr_window_20 - rData$VPD_yr_norm  # VPDex_yr_window_20 does not exist

#### add climate variables specific to the drought period (previous fall 2000 through previous fall 2004)
### cumulative precip during 2000-2003 drought
rData$PPT_drought <- rowSums(ppt.extr[, c("PPT_pf_2000", "PPT_pf_2001", "PPT_pf_2002", "PPT_pf_2003", "PPT_pf_2004",
                                                  "PPT_c_2000", "PPT_c_2001", "PPT_c_2002", "PPT_c_2003",
                                                  "PPT_fs_2000", "PPT_fs_2001", "PPT_fs_2002", "PPT_fs_2003",
                                                  "PPT_m_2000", "PPT_m_2001", "PPT_m_2002", "PPT_m_2003")])
# expected cum precip for this 50-month period
rData$PPT_50mo_norm <- 4*rowSums(rData[, c("PPT_pf_norm", "PPT_c_norm", "PPT_fs_norm", "PPT_m_norm")]) + rData$PPT_pf_norm

# 2000-2003 precip anomaly
rData$PPT_dr_anom <- rData$PPT_drought - rData$PPT_50mo_norm

# 2000-2003 monsoon precip may have been critical?
rData$PPT_m_dr <- rowSums(ppt.extr[, c("PPT_m_2000", "PPT_m_2001", "PPT_m_2002", "PPT_m_2003")])

rData$PPT_4m_norm <- 4*rData$PPT_m_norm
rData$PPT_m_dr_anom <- rData$PPT_m_dr - rData$PPT_4m_norm

# compare against effect of 2000-2003 cool season precip (=conventional wisdom)
rData$PPT_c_dr <- rowSums(ppt.extr[, c("PPT_c_2000", "PPT_c_2001", "PPT_c_2002", "PPT_c_2003")])

rData$PPT_4c_norm <- 4*rData$PPT_c_norm
rData$PPT_c_dr_anom <- rData$PPT_c_dr - rData$PPT_4c_norm

# previous fall
rData$PPT_pf_dr <- rowSums(ppt.extr[, c("PPT_pf_2000", "PPT_pf_2001", "PPT_pf_2002", "PPT_pf_2003")])

rData$PPT_4pf_norm <- 4*rData$PPT_pf_norm
rData$PPT_pf_dr_anom <- rData$PPT_pf_dr - rData$PPT_4pf_norm

# foresummer
rData$PPT_fs_dr <- rowSums(ppt.extr[, c("PPT_fs_2000", "PPT_fs_2001", "PPT_fs_2002", "PPT_fs_2003")])

rData$PPT_4fs_norm <- 4*rData$PPT_fs_norm
rData$PPT_fs_dr_anom <- rData$PPT_fs_dr - rData$PPT_4fs_norm

# T mean during 50 months of 2000-2003 drought
rData$Tmean_drought <- (tmp.extr$tmp_199909 + tmp.extr$tmp_199910 + tmp.extr$tmp_199911 + tmp.extr$tmp_199912 +
                                  tmp.extr$tmp_200001 + tmp.extr$tmp_200002 + tmp.extr$tmp_200003 + tmp.extr$tmp_200004 + tmp.extr$tmp_200005 + tmp.extr$tmp_200006 +
                                  tmp.extr$tmp_200007 + tmp.extr$tmp_200008 + tmp.extr$tmp_200009 + tmp.extr$tmp_200010 + tmp.extr$tmp_200011 + tmp.extr$tmp_200012 +  
                                  tmp.extr$tmp_200101 + tmp.extr$tmp_200102 + tmp.extr$tmp_200103 + tmp.extr$tmp_200104 + tmp.extr$tmp_200105 + tmp.extr$tmp_200106 +
                                  tmp.extr$tmp_200107 + tmp.extr$tmp_200108 + tmp.extr$tmp_200109 + tmp.extr$tmp_200110 + tmp.extr$tmp_200111 + tmp.extr$tmp_200112 +  
                                  tmp.extr$tmp_200201 + tmp.extr$tmp_200202 + tmp.extr$tmp_200203 + tmp.extr$tmp_200204 + tmp.extr$tmp_200205 + tmp.extr$tmp_200206 +
                                  tmp.extr$tmp_200207 + tmp.extr$tmp_200208 + tmp.extr$tmp_200209 + tmp.extr$tmp_200210 + tmp.extr$tmp_200211 + tmp.extr$tmp_200212 +  
                                  tmp.extr$tmp_200301 + tmp.extr$tmp_200302 + tmp.extr$tmp_200303 + tmp.extr$tmp_200304 + tmp.extr$tmp_200305 + tmp.extr$tmp_200306 +
                                  tmp.extr$tmp_200307 + tmp.extr$tmp_200308 + tmp.extr$tmp_200309 + tmp.extr$tmp_200310)/50  

# expected T mean during an equivalent 50-month window
rData$T_50mo_norm <- (4*rowSums(tmp.norm.extr[, c("tmp_1", "tmp_2", "tmp_3", "tmp_4",
                                                          "tmp_5", "tmp_6", "tmp_7", "tmp_8",
                                                          "tmp_9", "tmp_10", "tmp_11", "tmp_12")]) +
                                tmp.norm.extr[,9] + tmp.norm.extr[,10])/50

# anomaly (different for each plot)                                              
rData$T_dr_anom <- rData$Tmean_drought - rData$T_50mo_norm

# cumulative Tmean anomaly during 2000-2003 drought
rData$cum_T_anom <- ((tmp.extr$tmp_199909 - tmp.norm.extr[,9]) + (tmp.extr$tmp_199910 - tmp.norm.extr[,10]) + (tmp.extr$tmp_199911 - tmp.norm.extr[,11]) + (tmp.extr$tmp_199912 - tmp.norm.extr[,12]) +
                               (tmp.extr$tmp_200001 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200002 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200003 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200004 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200005 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200006 - tmp.norm.extr[,6]) +
                               (tmp.extr$tmp_200007 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200008 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200009 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200010 - tmp.norm.extr[,10]) + (tmp.extr$tmp_200011 - tmp.norm.extr[,11]) + (tmp.extr$tmp_200012 - tmp.norm.extr[,12]) +  
                               (tmp.extr$tmp_200101 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200102 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200103 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200104 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200105 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200106 - tmp.norm.extr[,6]) +
                               (tmp.extr$tmp_200107 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200108 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200109 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200110 - tmp.norm.extr[,10]) + (tmp.extr$tmp_200111 - tmp.norm.extr[,11]) + (tmp.extr$tmp_200112 - tmp.norm.extr[,12]) +  
                               (tmp.extr$tmp_200201 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200202 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200203 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200204 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200205 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200206 - tmp.norm.extr[,6]) +
                               (tmp.extr$tmp_200207 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200208 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200209 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200210 - tmp.norm.extr[,10]) + (tmp.extr$tmp_200211 - tmp.norm.extr[,11]) + (tmp.extr$tmp_200212 - tmp.norm.extr[,12]) +  
                               (tmp.extr$tmp_200301 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200302 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200303 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200304 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200305 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200306 - tmp.norm.extr[,6]) +
                               (tmp.extr$tmp_200307 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200308 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200309 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200310 - tmp.norm.extr[,10]))  


### last thing that needs to be done for recruitment subkernel
### get size distribution of recruits (ingrowth)
## it would be better to make sure size distribution is taken from the exact same set of trees that get put into ZIP model
## but it probably doesn't make a difference
recruits <- subset(trees, RECONCILECD == 1 & SPCD == 106)
r.sizemean <- mean(log(recruits$DIA))
r.sizesd <- sd(log(recruits$DIA))
# other size (state) variables are possible, but have been rejected at this point


##EXPORT##
rData2 <- rData[complete.cases(rData),]
#write.csv(rData, "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Recruitment/RecruitData.csv", row.names = F)
write.csv(rData, "./Processed/Recruitment/RecruitData.csv", row.names = F)

save(r.sizemean, r.sizesd, file = "./Code/IPM/recrstats.rda")

