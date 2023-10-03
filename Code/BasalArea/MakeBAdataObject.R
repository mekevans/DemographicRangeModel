library(sp)
library(raster)
library(dplyr)
library(effects)
library(mgcv)

# path
climate.path <- "./ClimateData/"

### create data object with BALIVE and climate normals ------------------------------------------------------------------------------------

conds <- read.csv("./FIAdata/COND_COMBINED.csv", header = T, stringsAsFactors = F)
# Remove non-forest conditions (missing BALIVE, coded as 0 so causes underestimation)
# conds <- subset(conds, COND_STATUS_CD == 1) # "accessible forest land" by FIA classification

# try including nonforested land in model
conds.1 <- subset(conds, COND_STATUS_CD == 1) # "accessible forest land" by FIA classification
conds.2 <- subset(conds, COND_STATUS_CD == 2) # "nonforest land" by FIA classification
conds.2 <- conds.2[conds.2$PRESNFCD %in% c(20, 45), ] # only two categories of non-forest land that are not caused by (human) land use
conds.2$BALIVE <- 0 # populate with zeros

# join the two together
conds <- rbind(conds.1, conds.2)

hist(conds$BALIVE) # almost all of the data has values <=400; outliers above 400
summary(conds$BALIVE) # maximum value is 1363!
plot(conds$CONDPROP_UNADJ, conds$BALIVE) # reveals that one outlier for BALIVE is a cond with very low CONDPROP_UNADJ
conds <- subset(conds, BALIVE < 1000) # removes one outlier!
# another data point that perhaps should be removed is #2919 (see resids vs. leverage)

# Read in plot data and get coordinates and previous measurement year
plots <- read.csv("./FIAdata/PLOT_COMBINED.csv", header = T, stringsAsFactors = F)

conds$LAT <- plots$LAT[match(conds$PLT_CN, plots$CN)]
conds$LON <- plots$LON[match(conds$PLT_CN, plots$CN)]
conds$ELEV <- plots$ELEV[match(conds$PLT_CN, plots$CN)]
conds$MEASYEAR <- plots$MEASYEAR[match(conds$PLT_CN, plots$CN)]

# Make data spatial
CondsSpat <- SpatialPointsDataFrame(coords = cbind(conds$LON, conds$LAT), 
                                    data = conds, 
                                    proj4string = CRS("+proj=longlat +datum=NAD83"))

# import PRISM normals
PPT.norm <- stack(paste0(climate.path,"pptNormals.tif",sep=''))
TMP.norm <- stack(paste0(climate.path,"tmpNormals.tif",sep=''))
VPD.norm <- stack(paste0(climate.path,"vpdNormals.tif",sep=''))

# raster::extract PRISM normals 1981-2010
ppt.norm.extr <- raster::extract(PPT.norm, CondsSpat)
tmp.norm.extr <- raster::extract(TMP.norm, CondsSpat)
vpd.norm.extr <- raster::extract(VPD.norm, CondsSpat)

ppt.norm.extr <- as.data.frame(ppt.norm.extr) # not sure this is necessary?
tmp.norm.extr <- as.data.frame(tmp.norm.extr)
vpd.norm.extr <- as.data.frame(vpd.norm.extr)

# reasonable column names (not actually necessary)
colnames(ppt.norm.extr) <- paste0("ppt_", 1:12) 
colnames(tmp.norm.extr) <- paste0("tmp_", 1:12)
colnames(vpd.norm.extr) <- paste0("vpd_", 1:12)

# make seasonal normals and add to conditions data frame
# cool season = Nov-Mar
conds$PPT_c_norm <- rowSums(ppt.norm.extr[, c(1:3, 11:12)])
conds$T_c_norm <- rowMeans(tmp.norm.extr[, c(1:3, 11:12)])
conds$VPD_c_norm <- rowMeans(vpd.norm.extr[, c(1:3, 11:12)])
# previous fall = pSept-pOct
conds$PPT_pf_norm <- rowSums(ppt.norm.extr[, c(9:10)])
conds$T_pf_norm <- rowMeans(tmp.norm.extr[, c(9:10)])
conds$VPD_pf_norm <- rowMeans(vpd.norm.extr[, c(9:10)])
# foresummer = Apr-Jun
conds$PPT_fs_norm <- rowSums(ppt.norm.extr[, c(4:6)])
conds$T_fs_norm <- rowMeans(tmp.norm.extr[, c(4:6)])
conds$VPD_fs_norm <- rowMeans(vpd.norm.extr[, c(4:6)])
# warm, dry months = Apr-Jun + Sept-Oct
conds$PPT_wd_norm <- rowSums(ppt.norm.extr[, c(4:6, 9:10)])
conds$T_wd_norm <- rowMeans(tmp.norm.extr[, c(4:6, 9:10)])
conds$VPD_wd_norm <- rowMeans(vpd.norm.extr[, c(4:6, 9:10)])
# monsoon = Jul-Aug
conds$PPT_m_norm <- rowSums(ppt.norm.extr[, c(7:8)])
conds$T_m_norm <- rowMeans(tmp.norm.extr[, c(7:8)])
conds$VPD_m_norm <- rowMeans(vpd.norm.extr[, c(7:8)])
# water year
conds$PPT_yr_norm <- rowSums(ppt.norm.extr[, c(1:12)])
conds$T_yr_norm <- rowMeans(tmp.norm.extr[, c(1:12)])
conds$VPD_yr_norm <- rowMeans(vpd.norm.extr[, c(1:12)])

# Create output data frame
output <- conds[, c("PLT_CN", "CONDID", "CONDPROP_UNADJ",
                    "LAT", "LON", "ELEV",
                    "BALIVE", 
                    "PPT_c_norm", "T_c_norm", "VPD_c_norm",
                    "PPT_wd_norm", "T_wd_norm", "VPD_wd_norm",
                    "PPT_pf_norm", "T_pf_norm", "VPD_pf_norm",
                    "PPT_fs_norm", "T_fs_norm", "VPD_fs_norm",
                    "PPT_m_norm", "T_m_norm", "VPD_m_norm",
                    "PPT_yr_norm", "T_yr_norm", "VPD_yr_norm")]
# write.csv(output, paste0(path, "BA/BALIVEdata.csv"), row.names = F)
write.csv(output,"./BA/BALIVEdata2.csv", row.names = F)
