library(sp)
library(raster)
library(ggplot2)
library(wesanderson)

data.path <- "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/FIAdata/"

# Read data and subset PIED
grData <- read.csv(paste(data.path,"TREE_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
grData <- subset(grData, SPCD == 106)

# Only keep remeasured trees coded as dead or alive
    ### STATUSCD = 0 = no status (not in sample)
    ### STATUSCD = 1 = live tree
    ### STATUSCD = 2 = dead tree
    ### STATUSCD = 3 = harvested tree
grData_remeas <- subset(grData, !is.na(PREVDIA))
grData_remeas <- subset(grData_remeas, STATUSCD == 1 | STATUSCD == 2) 

# Look up previous AGB
### look up previous AGB, PLT_CN, and CONDID
grData_remeas$PREV_DRYBIO_AG <- grData$DRYBIO_AG[match(grData_remeas$PREV_TRE_CN, grData$CN)]
grData_remeas$PREV_PLT_CN <- grData$PLT_CN[match(grData_remeas$PREV_TRE_CN, grData$CN)]
grData_remeas$PREV_CONDID <- grData$CONDID[match(grData_remeas$PREV_TRE_CN, grData$CN)]

# Subset records with previous AGB found
grData_remeas <- subset(grData_remeas, !is.na(PREV_DRYBIO_AG))

# for growth analysis
grData_remeas$DRYBIO_AG_DIFF <- grData_remeas$DRYBIO_AG - grData_remeas$PREV_DRYBIO_AG
grData_remeas$DIA_DIFF <- grData_remeas$DIA - grData_remeas$PREVDIA

# basal area increment
grData_remeas$BAt1 <- ((grData_remeas$PREVDIA/2)^2)*3.14159
grData_remeas$BAt2 <- ((grData_remeas$DIA/2)^2)*3.14159
grData_remeas$BA_DIFF <- grData_remeas$BAt2 - grData_remeas$BAt1

# Read in plot data and get coordinates and previous measurement year
plots <- read.csv(paste(data.path,"PLOT_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)

grData_remeas$LAT <- plots$LAT[match(grData_remeas$PLT_CN, plots$CN)]
grData_remeas$LON <- plots$LON[match(grData_remeas$PLT_CN, plots$CN)]
grData_remeas$ELEV <- plots$ELEV[match(grData_remeas$PLT_CN, plots$CN)]
grData_remeas$MEASYEAR <- plots$MEASYEAR[match(grData_remeas$PLT_CN, plots$CN)]
grData_remeas$PREV_MEASYEAR <- plots$MEASYEAR[match(grData_remeas$PREV_PLT_CN, plots$CN)]

# Calculate census interval
grData_remeas$CENSUS_INTERVAL <- grData_remeas$MEASYEAR - grData_remeas$PREV_MEASYEAR

# look up previous (tree-specific) condition-level BALIVE
conds <- read.csv(paste(data.path,"COND_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
conds <- subset(conds, COND_STATUS_CD == 1) # "accessible forest land" by FIA classification
grData_remeas$BALIVE <- apply(X = grData_remeas[, c("PREV_PLT_CN", "PREV_CONDID")], 
                              MARGIN = 1, # applies function to each row in grData_remeas
                              FUN = function(x, conds.df) {
                                conds.df$BALIVE[conds.df$PLT_CN %in% x["PREV_PLT_CN"] &
                                                  conds.df$CONDID %in% x["PREV_CONDID"]]
                              },
                              conds.df = conds)
grData_remeas[is.nan(grData_remeas$BALIVE), "BALIVE"] <- NA
grData_remeas <- subset(grData_remeas, !is.na(BALIVE))

### look at mortality, by MEASYEAR and PREV_MEASYEAR
surv.tablet1 <- table(grData_remeas[, c("PREV_MEASYEAR", "STATUSCD")])
surv.tablet1[,1]/(surv.tablet1[,1]+surv.tablet1[,2])
# survival ranges from 0.738 to 0.810, no real temporal resolution

surv.tablet2 <- table(grData_remeas[, c("MEASYEAR", "STATUSCD")])
surv.tablet2[,1]/(surv.tablet2[,1]+surv.tablet2[,2])
# survival ranges from 0.751 to 0.811

# Make data spatial
grSpat <- SpatialPointsDataFrame(coords = cbind(grData_remeas$LON, grData_remeas$LAT), 
                                 data = grData_remeas, 
                                 proj4string = CRS("+proj=longlat +datum=NAD83"))

# Read and extract climate
# THE FOLLOWING LINES CAN BE SKIPPED if you've already generated the "ppt_extr.csv" etc. files
### should be moved to historic.R
### should be done once for both the survival and growth data, and dead trees subsetted out (one line) for growth analysis

# Read in PRISM climate stacks
clim.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/"
ppt <- stack(paste(clim.path,"pptStack.tif",sep=''))
tmp <- stack(paste(clim.path,"tmpStack.tif",sep=''))
vpd <- stack(paste(clim.path,"vpdStack.tif",sep=''))

# raster::extract PRISM data
ppt.extr <- raster::extract(ppt, grSpat)
tmp.extr <- raster::extract(tmp, grSpat)
vpd.extr <- raster::extract(vpd, grSpat)

# Remove data after Oct, 2016 (because of different CRS Nov, 2016 vpdmax .bil)
# note that the work-around for this problem is to assign the CRS of another layer to Nov and Dec of 2016
# crs(vpdNov2016_raster) <- crs(vpdOct2016_raster)
ppt.extr <- ppt.extr[, 1:430] 
tmp.extr <- tmp.extr[, 1:430]
vpd.extr <- vpd.extr[, 1:430]

# Add sensible column names for raster::extracted climate data
ppt.extr <- as.data.frame(ppt.extr)
tmp.extr <- as.data.frame(tmp.extr)
vpd.extr <- as.data.frame(vpd.extr)
PRISM.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/"
pptFiles <- list.files(path = PRISM.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
pptFiles <- pptFiles[1:430] # (hack to deal with CRS incompatibility, vpd .bil file Nov, 2016)
#tmpFiles <- list.files(path = PRISM.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
#vpdFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
colNames <- lapply(strsplit(pptFiles, "4kmM._"), function (x) x[2])
colNames <- unlist(colNames)
colNames <- lapply(strsplit(colNames, "_"), function (x) x[1])
colNames <- unlist(colNames)
colnames(ppt.extr) <- paste0("ppt_", colNames)
colnames(tmp.extr) <- paste0("tmp_", colNames)
colnames(vpd.extr) <- paste0("vpd_", colNames)

# Export climate data
surv.path <- "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Survival/"
write.csv(ppt.extr, paste0(surv.path, "ppt_extr.csv"), row.names = F)
write.csv(tmp.extr, paste0(surv.path, "tmp_extr.csv"), row.names = F)
write.csv(vpd.extr, paste0(surv.path, "vpd_extr.csv"), row.names = F)

### CONTINUE HERE IF YOU ALREADY HAVE GENERATED THE ABOVE csv files

ppt.extr <- read.csv(paste(surv.path,"ppt_extr.csv",sep=''), header = T)
tmp.extr <- read.csv(paste(surv.path,"tmp_extr.csv",sep=''), header = T)
vpd.extr <- read.csv(paste(surv.path,"vpd_extr.csv",sep=''), header = T)

# Calculate seasonal climate for each year
for (i in 1982:2016) {
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

### THE FOLLOWING ONLY NEEDS TO BE DONE ONCE
### should be moved to normals.R
# import PRISM normals
PRISM.norm.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/Normals/"
PPT.norm <- stack(paste(clim.path,"pptNormals.tif",sep='')) # should this actually read "PRISM.norm.path" not "clim.path"??
TMP.norm <- stack(paste(clim.path,"tmpNormals.tif",sep=''))
VPD.norm <- stack(paste(clim.path,"vpdNormals.tif",sep=''))

# raster::extract PRISM normals 1981-2010
ppt.norm.extr <- raster::extract(PPT.norm, grSpat)
tmp.norm.extr <- raster::extract(TMP.norm, grSpat)
vpd.norm.extr <- raster::extract(VPD.norm, grSpat)

ppt.norm.extr <- as.data.frame(ppt.norm.extr) # not sure this is necessary?
tmp.norm.extr <- as.data.frame(tmp.norm.extr)
vpd.norm.extr <- as.data.frame(vpd.norm.extr)

# reasonable column names (not actually necessary)
colnames(ppt.norm.extr) <- paste0("ppt_", 1:12) 
colnames(tmp.norm.extr) <- paste0("tmp_", 1:12)
colnames(vpd.norm.extr) <- paste0("vpd_", 1:12)

# make seasonal normals and add to growth data frame
# cool season = Nov-Mar
grData_remeas$PPT_c_norm <- rowSums(ppt.norm.extr[, c(1:3, 11:12)])
grData_remeas$T_c_norm <- rowMeans(tmp.norm.extr[, c(1:3, 11:12)])
grData_remeas$VPD_c_norm <- rowMeans(vpd.norm.extr[, c(1:3, 11:12)])
# previous fall = pSept-pOct
grData_remeas$PPT_pf_norm <- rowSums(ppt.norm.extr[, c(9:10)])
grData_remeas$T_pf_norm <- rowMeans(tmp.norm.extr[, c(9:10)])
grData_remeas$VPD_pf_norm <- rowMeans(vpd.norm.extr[, c(9:10)])
# foresummer = Apr-Jun
grData_remeas$PPT_fs_norm <- rowSums(ppt.norm.extr[, c(4:6)])
grData_remeas$T_fs_norm <- rowMeans(tmp.norm.extr[, c(4:6)])
grData_remeas$VPD_fs_norm <- rowMeans(vpd.norm.extr[, c(4:6)])
# warm, dry months = Apr-Jun + Sept-Oct
grData_remeas$PPT_wd_norm <- rowSums(ppt.norm.extr[, c(4:6, 9:10)])
grData_remeas$T_wd_norm <- rowMeans(tmp.norm.extr[, c(4:6, 9:10)])
grData_remeas$VPD_wd_norm <- rowMeans(vpd.norm.extr[, c(4:6, 9:10)])
# monsoon = Jul-Aug
grData_remeas$PPT_m_norm <- rowSums(ppt.norm.extr[, c(7:8)])
grData_remeas$T_m_norm <- rowMeans(tmp.norm.extr[, c(7:8)])
grData_remeas$VPD_m_norm <- rowMeans(vpd.norm.extr[, c(7:8)])
# water year
grData_remeas$PPT_yr_norm <- rowSums(ppt.norm.extr[, c(1:12)])
grData_remeas$T_yr_norm <- rowMeans(tmp.norm.extr[, c(1:12)])
grData_remeas$VPD_yr_norm <- rowMeans(vpd.norm.extr[, c(1:12)])


# Add seasonal climate variables (specific to each tree's census interval) to growth data frame
# could not get this to vectorize successfully
# close is:
# grData_remeas$PPT_c <- rowMeans(ppt.extr[, paste0("PPT_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPT_c <- rowMeans(ppt.extr[, paste0("PPT_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$T_c <- rowMeans(tmp.extr[, paste0("T_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$VPD_c <- rowMeans(vpd.extr[, paste0("VPD_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPTex_c <- apply(ppt.extr[, paste0("PPT_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, min)
grData_remeas$Tex_c <- apply(tmp.extr[, paste0("T_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)
grData_remeas$VPDex_c <- apply(vpd.extr[, paste0("VPD_c_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)

grData_remeas$PPT_pf <- rowMeans(ppt.extr[, paste0("PPT_pf_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$T_pf <- rowMeans(tmp.extr[, paste0("T_pf_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$VPD_pf <- rowMeans(vpd.extr[, paste0("VPD_pf_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPTex_pf <- apply(ppt.extr[, paste0("PPT_pf_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, min)
grData_remeas$Tex_pf <- apply(tmp.extr[, paste0("T_pf_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)
grData_remeas$VPDex_pf <- apply(vpd.extr[, paste0("VPD_pf_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)

grData_remeas$PPT_fs <- rowMeans(ppt.extr[, paste0("PPT_fs_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$T_fs <- rowMeans(tmp.extr[, paste0("T_fs_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$VPD_fs <- rowMeans(vpd.extr[, paste0("VPD_fs_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPTex_fs <- apply(ppt.extr[, paste0("PPT_fs_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, min)
grData_remeas$Tex_fs <- apply(tmp.extr[, paste0("T_fs_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)
grData_remeas$VPDex_fs <- apply(vpd.extr[, paste0("VPD_fs_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)

grData_remeas$PPT_wd <- rowMeans(ppt.extr[, paste0("PPT_wd_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$T_wd <- rowMeans(tmp.extr[, paste0("T_wd_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$VPD_wd <- rowMeans(vpd.extr[, paste0("VPD_wd_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPTex_wd <- apply(ppt.extr[, paste0("PPT_wd_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, min)
grData_remeas$Tex_wd <- apply(tmp.extr[, paste0("T_wd_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)
grData_remeas$VPDex_wd <- apply(vpd.extr[, paste0("VPD_wd_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)

grData_remeas$PPT_m <- rowMeans(ppt.extr[, paste0("PPT_m_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$T_m <- rowMeans(tmp.extr[, paste0("T_m_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$VPD_m <- rowMeans(vpd.extr[, paste0("VPD_m_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPTex_m <- apply(ppt.extr[, paste0("PPT_m_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]))], 1, min)
grData_remeas$Tex_m <- apply(tmp.extr[, paste0("T_m_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]))], 1, max)
grData_remeas$VPDex_m <- apply(vpd.extr[, paste0("VPD_m_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]))], 1, max)

grData_remeas$PPT_yr <- rowMeans(ppt.extr[, paste0("PPT_yr_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$T_yr <- rowMeans(tmp.extr[, paste0("T_yr_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$VPD_yr <- rowMeans(vpd.extr[, paste0("VPD_yr_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))])
grData_remeas$PPTex_yr <- apply(ppt.extr[, paste0("PPT_yr_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, min)
grData_remeas$Tex_yr <- apply(tmp.extr[, paste0("T_yr_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)
grData_remeas$VPDex_yr <- apply(vpd.extr[, paste0("VPD_yr_", (grData_remeas[i, "PREV_MEASYEAR"]+1):(grData_remeas[i, "MEASYEAR"]))], 1, max)

# create and add anomalies to growth data frame
# note that precip is lognormal (R-skewed), but precip anomalies are left-skewed (heavy tail in the negative)...because mean of a lognormal is right of center
# monsoon precip anomalies are centered on zero (no change)
# cool season, foresummer, and water year anomalies are left-skewed (heavy tail in negative values)
# temperature anomalies are positive (no surprise)
# cool season (pNov-Mar)
grData_remeas$PPT_c_anom <- grData_remeas$PPT_c - grData_remeas$PPT_c_norm # negative
grData_remeas$T_c_anom <- grData_remeas$T_c - grData_remeas$T_c_norm # positive
grData_remeas$VPD_c_anom <- grData_remeas$VPD_c - grData_remeas$VPD_c_norm # positive
grData_remeas$PPTex_c_anom <- grData_remeas$PPTex_c - grData_remeas$PPT_c_norm 
grData_remeas$Tex_c_anom <- grData_remeas$Tex_c - grData_remeas$T_c_norm
grData_remeas$VPDex_c_anom <- grData_remeas$VPDex_c - grData_remeas$VPD_c_norm

# previous fall (pSep-pOct)
grData_remeas$PPT_pf_anom <- grData_remeas$PPT_pf - grData_remeas$PPT_pf_norm # positive anomaly!
grData_remeas$T_pf_anom <- grData_remeas$T_pf - grData_remeas$T_pf_norm # positive
grData_remeas$VPD_pf_anom <- grData_remeas$VPD_pf - grData_remeas$VPD_pf_norm # positive
grData_remeas$PPTex_pf_anom <- grData_remeas$PPTex_pf - grData_remeas$PPT_pf_norm # negative
grData_remeas$Tex_pf_anom <- grData_remeas$Tex_pf - grData_remeas$T_pf_norm # positive
grData_remeas$VPDex_pf_anom <- grData_remeas$VPDex_pf - grData_remeas$VPD_pf_norm # positive
# foresummer (Apr-Jun)
grData_remeas$PPT_fs_anom <- grData_remeas$PPT_fs - grData_remeas$PPT_fs_norm # positive!
grData_remeas$T_fs_anom <- grData_remeas$T_fs - grData_remeas$T_fs_norm # positive
grData_remeas$VPD_fs_anom <- grData_remeas$VPD_fs - grData_remeas$VPD_fs_norm # positive
grData_remeas$PPTex_fs_anom <- grData_remeas$PPTex_fs - grData_remeas$PPT_fs_norm # centered on zero
grData_remeas$Tex_fs_anom <- grData_remeas$Tex_fs - grData_remeas$T_fs_norm # positive
grData_remeas$VPDex_fs_anom <- grData_remeas$VPDex_fs - grData_remeas$VPD_fs_norm # positive
# warm dry months (pSep-pOct + Apr-Jun)
grData_remeas$PPT_wd_anom <- grData_remeas$PPT_wd - grData_remeas$PPT_wd_norm # positive!
grData_remeas$T_wd_anom <- grData_remeas$T_wd - grData_remeas$T_wd_norm
grData_remeas$VPD_wd_anom <- grData_remeas$VPD_wd - grData_remeas$VPD_wd_norm
grData_remeas$PPTex_wd_anom <- grData_remeas$PPTex_wd - grData_remeas$PPT_wd_norm # positive!
grData_remeas$Tex_wd_anom <- grData_remeas$Tex_wd - grData_remeas$T_wd_norm
grData_remeas$VPDex_wd_anom <- grData_remeas$VPDex_wd - grData_remeas$VPD_wd_norm
# monsoon (Jul-Aug)
grData_remeas$PPT_m_anom <- grData_remeas$PPT_m - grData_remeas$PPT_m_norm # positive!
grData_remeas$T_m_anom <- grData_remeas$T_m - grData_remeas$T_m_norm
grData_remeas$VPD_m_anom <- grData_remeas$VPD_m - grData_remeas$VPD_m_norm
grData_remeas$PPTex_m_anom <- grData_remeas$PPTex_m - grData_remeas$PPT_m_norm # centered on zero
grData_remeas$Tex_m_anom <- grData_remeas$Tex_m - grData_remeas$T_m_norm
grData_remeas$VPDex_m_anom <- grData_remeas$VPDex_m - grData_remeas$VPD_m_norm
# water year (pSept - Aug)
grData_remeas$PPT_yr_anom <- grData_remeas$PPT_yr - grData_remeas$PPT_yr_norm # negative
grData_remeas$T_yr_anom <- grData_remeas$T_yr - grData_remeas$T_yr_norm
grData_remeas$VPD_yr_anom <- grData_remeas$VPD_yr - grData_remeas$VPD_yr_norm
grData_remeas$PPTex_yr_anom <- grData_remeas$PPTex_yr - grData_remeas$PPT_yr_norm
grData_remeas$Tex_yr_anom <- grData_remeas$Tex_yr - grData_remeas$T_yr_norm
grData_remeas$VPDex_yr_anom <- grData_remeas$VPDex_yr - grData_remeas$VPD_yr_norm

#### add climate variables specific to the drought period (previous fall 2000 through previous fall 2004)
### cumulative precip during 2000-2003 drought
grData_remeas$PPT_drought <- rowSums(ppt.extr[, c("PPT_pf_2000", "PPT_pf_2001", "PPT_pf_2002", "PPT_pf_2003", "PPT_pf_2004",
                                                  "PPT_c_2000", "PPT_c_2001", "PPT_c_2002", "PPT_c_2003",
                                                  "PPT_fs_2000", "PPT_fs_2001", "PPT_fs_2002", "PPT_fs_2003",
                                                  "PPT_m_2000", "PPT_m_2001", "PPT_m_2002", "PPT_m_2003")])
# expected cum precip for this 50-month period
grData_remeas$PPT_50mo_norm <- 4*rowSums(grData_remeas[, c("PPT_pf_norm", "PPT_c_norm", "PPT_fs_norm", "PPT_m_norm")]) + grData_remeas$PPT_pf_norm

# 2000-2003 precip anomaly
grData_remeas$PPT_dr_anom <- grData_remeas$PPT_drought - grData_remeas$PPT_50mo_norm

# 2000-2003 monsoon precip may have been critical?
grData_remeas$PPT_m_dr <- rowSums(ppt.extr[, c("PPT_m_2000", "PPT_m_2001", "PPT_m_2002", "PPT_m_2003")])

grData_remeas$PPT_4m_norm <- 4*grData_remeas$PPT_m_norm
grData_remeas$PPT_m_dr_anom <- grData_remeas$PPT_m_dr - grData_remeas$PPT_4m_norm

# compare against effect of 2000-2003 cool season precip (=conventional wisdom)
grData_remeas$PPT_c_dr <- rowSums(ppt.extr[, c("PPT_c_2000", "PPT_c_2001", "PPT_c_2002", "PPT_c_2003")])

grData_remeas$PPT_4c_norm <- 4*grData_remeas$PPT_c_norm
grData_remeas$PPT_c_dr_anom <- grData_remeas$PPT_c_dr - grData_remeas$PPT_4c_norm

# previous fall
grData_remeas$PPT_pf_dr <- rowSums(ppt.extr[, c("PPT_pf_2000", "PPT_pf_2001", "PPT_pf_2002", "PPT_pf_2003")])

grData_remeas$PPT_4pf_norm <- 4*grData_remeas$PPT_pf_norm
grData_remeas$PPT_pf_dr_anom <- grData_remeas$PPT_pf_dr - grData_remeas$PPT_4pf_norm

# foresummer
grData_remeas$PPT_fs_dr <- rowSums(ppt.extr[, c("PPT_fs_2000", "PPT_fs_2001", "PPT_fs_2002", "PPT_fs_2003")])

grData_remeas$PPT_4fs_norm <- 4*grData_remeas$PPT_fs_norm
grData_remeas$PPT_fs_dr_anom <- grData_remeas$PPT_fs_dr - grData_remeas$PPT_4fs_norm

# T mean during 50 months of 2000-2003 drought
grData_remeas$Tmean_drought <- (tmp.extr$tmp_199909 + tmp.extr$tmp_199910 + tmp.extr$tmp_199911 + tmp.extr$tmp_199912 +
                  tmp.extr$tmp_200001 + tmp.extr$tmp_200002 + tmp.extr$tmp_200003 + tmp.extr$tmp_200004 + tmp.extr$tmp_200005 + tmp.extr$tmp_200006 +
                  tmp.extr$tmp_200007 + tmp.extr$tmp_200008 + tmp.extr$tmp_200009 + tmp.extr$tmp_200010 + tmp.extr$tmp_200011 + tmp.extr$tmp_200012 +  
                  tmp.extr$tmp_200101 + tmp.extr$tmp_200102 + tmp.extr$tmp_200103 + tmp.extr$tmp_200104 + tmp.extr$tmp_200105 + tmp.extr$tmp_200106 +
                  tmp.extr$tmp_200107 + tmp.extr$tmp_200108 + tmp.extr$tmp_200109 + tmp.extr$tmp_200110 + tmp.extr$tmp_200111 + tmp.extr$tmp_200112 +  
                  tmp.extr$tmp_200201 + tmp.extr$tmp_200202 + tmp.extr$tmp_200203 + tmp.extr$tmp_200204 + tmp.extr$tmp_200205 + tmp.extr$tmp_200206 +
                  tmp.extr$tmp_200207 + tmp.extr$tmp_200208 + tmp.extr$tmp_200209 + tmp.extr$tmp_200210 + tmp.extr$tmp_200211 + tmp.extr$tmp_200212 +  
                  tmp.extr$tmp_200301 + tmp.extr$tmp_200302 + tmp.extr$tmp_200303 + tmp.extr$tmp_200304 + tmp.extr$tmp_200305 + tmp.extr$tmp_200306 +
                  tmp.extr$tmp_200307 + tmp.extr$tmp_200308 + tmp.extr$tmp_200309 + tmp.extr$tmp_200310)/50  
                    
# expected T mean during an equivalent 50-month window
grData_remeas$T_50mo_norm <- (4*rowSums(tmp.norm.extr[, c("tmp_1", "tmp_2", "tmp_3", "tmp_4",
                                         "tmp_5", "tmp_6", "tmp_7", "tmp_8",
                                         "tmp_9", "tmp_10", "tmp_11", "tmp_12")]) +
              tmp.norm.extr[,9] + tmp.norm.extr[,10])/50

# anomaly (different for each plot)                                              
grData_remeas$T_dr_anom <- grData_remeas$Tmean_drought - grData_remeas$T_50mo_norm

# cumulative Tmean anomaly during 2000-2003 drought
grData_remeas$cum_T_anom <- ((tmp.extr$tmp_199909 - tmp.norm.extr[,9]) + (tmp.extr$tmp_199910 - tmp.norm.extr[,10]) + (tmp.extr$tmp_199911 - tmp.norm.extr[,11]) + (tmp.extr$tmp_199912 - tmp.norm.extr[,12]) +
              (tmp.extr$tmp_200001 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200002 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200003 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200004 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200005 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200006 - tmp.norm.extr[,6]) +
              (tmp.extr$tmp_200007 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200008 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200009 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200010 - tmp.norm.extr[,10]) + (tmp.extr$tmp_200011 - tmp.norm.extr[,11]) + (tmp.extr$tmp_200012 - tmp.norm.extr[,12]) +  
              (tmp.extr$tmp_200101 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200102 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200103 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200104 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200105 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200106 - tmp.norm.extr[,6]) +
              (tmp.extr$tmp_200107 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200108 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200109 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200110 - tmp.norm.extr[,10]) + (tmp.extr$tmp_200111 - tmp.norm.extr[,11]) + (tmp.extr$tmp_200112 - tmp.norm.extr[,12]) +  
              (tmp.extr$tmp_200201 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200202 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200203 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200204 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200205 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200206 - tmp.norm.extr[,6]) +
              (tmp.extr$tmp_200207 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200208 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200209 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200210 - tmp.norm.extr[,10]) + (tmp.extr$tmp_200211 - tmp.norm.extr[,11]) + (tmp.extr$tmp_200212 - tmp.norm.extr[,12]) +  
              (tmp.extr$tmp_200301 - tmp.norm.extr[,1]) + (tmp.extr$tmp_200302 - tmp.norm.extr[,2]) + (tmp.extr$tmp_200303 - tmp.norm.extr[,3]) + (tmp.extr$tmp_200304 - tmp.norm.extr[,4]) +  (tmp.extr$tmp_200305 - tmp.norm.extr[,5]) + (tmp.extr$tmp_200306 - tmp.norm.extr[,6]) +
              (tmp.extr$tmp_200307 - tmp.norm.extr[,7]) + (tmp.extr$tmp_200308 - tmp.norm.extr[,8]) + (tmp.extr$tmp_200309 - tmp.norm.extr[,9]) + (tmp.extr$tmp_200310 - tmp.norm.extr[,10]))  

# make plot of survival rate as a function of cumulative PPT during drought (and anomaly)...like Clifford et al
# make plot of survival rate as a function of cumulative T anomaly during drought...like Clifford et al
# thanks to Jeff Oliver for help
surv.table <- grData_remeas %>%
  group_by(PLT_CN, CONDID) %>%
  summarise(survival.rate = sum(STATUSCD ==1)/n())
cols.to.include <- c("PLT_CN", "CONDID", "BALIVE", "ELEV", "LAT", 
                     "PPT_pf", "T_pf", "VPD_pf", "PPT_pf_anom", "PPT_pf_dr", "PPT_pf_dr_anom", 
                     "PPT_c", "T_c", "VPD_c", "PPT_c_anom", "PPT_c_dr", "PPT_c_dr_anom", 
                     "PPT_fs", "T_fs", "VPD_fs", "PPT_fs_anom", "PPT_fs_dr", "PPT_fs_dr_anom", 
                     "PPT_m", "T_m", "VPD_m", "PPT_m_anom", "PPT_m_dr", "PPT_m_dr_anom",
                     "PPT_yr_norm", "T_yr_norm", "VPD_yr_norm", # normals
                     "PPT_drought", "PPT_dr_anom", "Tmean_drought", "T_dr_anom", "cum_T_anom") # 2000-2003 drought
surv.table <- merge(x = surv.table,
                       y = grData_remeas[!duplicated(grData_remeas[c("PLT_CN", "CONDID")], ), cols.to.include],
                       by = c("PLT_CN", "CONDID"))

surv_competition <- ggplot(data = surv.table, aes(y = survival.rate, x = BALIVE)) + geom_point(aes(col = ELEV)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_competition + scale_color_gradient(low="brown", high="chartreuse1")

surv_elevation <- ggplot(data = surv.table, aes(y = survival.rate, x = ELEV)) + geom_point(aes(col = BALIVE)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_elevation + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptyr <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_yr_norm)) + geom_point(aes(col = ELEV)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptyr + scale_color_gradient(low="brown", high="chartreuse1")

surv_tyr <- ggplot(data = surv.table, aes(y = survival.rate, x = T_yr_norm)) + geom_point(aes(col = ELEV)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_tyr + scale_color_gradient(low="brown", high="chartreuse1")

surv_tyr <- ggplot(data = surv.table, aes(y = survival.rate, x = T_yr_norm)) + geom_point(aes(col = PPT_yr_norm)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_tyr + scale_color_gradient(low="brown", high="chartreuse1")

surv_tyr <- ggplot(data = surv.table, aes(y = survival.rate, x = T_yr_norm)) + geom_point(aes(col = PPT_drought)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_tyr + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptdr <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_drought)) + geom_point(aes(col = ELEV)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptdr + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptdranom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_dr_anom)) + geom_point(aes(col = cum_T_anom)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptdranom + scale_color_gradient(low="chartreuse1", high="brown")

surv_pptm <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_m)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptm + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptmanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_m_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptmanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptmanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_m_dr_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptmanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptmdr <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_m_dr)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptmdr + scale_color_gradient(low="brown", high="chartreuse1")

### check on effect of cool season precip
surv_pptc <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_c)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptc + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptcanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_c_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptcanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptcanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_c_dr_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptcanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptcdr <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_c_dr)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptcdr + scale_color_gradient(low="brown", high="chartreuse1")

# effect of previous fall precip?
surv_pptpf <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_pf)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptpf + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptpfanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_pf_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptpfanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptpfanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_pf_dr_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptpfanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptpfdr <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_pf_dr)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptpfdr + scale_color_gradient(low="brown", high="chartreuse1")

# effect of foresummer precip?
surv_pptfs <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_fs)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptfs + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptfsanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_fs_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptfsanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptfsanom <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_fs_dr_anom)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptfsanom + scale_color_gradient(low="brown", high="chartreuse1")

surv_pptfsdr <- ggplot(data = surv.table, aes(y = survival.rate, x = PPT_fs_dr)) + geom_point(aes(col = LAT)) + stat_smooth(method = "glm", method.args = list(family = "binomial"))
surv_pptfsdr + scale_color_gradient(low="brown", high="chartreuse1")


# write code to search for min cum precip within census interval + preceeding 10 years (1-yr min, 2-yr min, 3-yr min, 4-yr min, 5-yr min, etc.)
# and cumulative T anomaly in that same window

# examine anomalies - how intense was the drought, and what spatial variation (distribution)?
PPT_pf_1999_anom <- ppt.extr$PPT_pf_1999-grData_remeas$PPT_pf_norm # 98 = wet
PPT_pf_2000_anom <- ppt.extr$PPT_pf_2000-grData_remeas$PPT_pf_norm # 99 = dry
PPT_pf_2001_anom <- ppt.extr$PPT_pf_2001-grData_remeas$PPT_pf_norm # 00 = wet
PPT_pf_2002_anom <- ppt.extr$PPT_pf_2002-grData_remeas$PPT_pf_norm # 01 = dry
PPT_pf_2003_anom <- ppt.extr$PPT_pf_2003-grData_remeas$PPT_pf_norm # 02 = wet
PPT_pf_2004_anom <- ppt.extr$PPT_pf_2004-grData_remeas$PPT_pf_norm # 03 = dry

PPT_c_1999_anom <- ppt.extr$PPT_c_1999-grData_remeas$PPT_c_norm # 98-99 = dry
PPT_c_2000_anom <- ppt.extr$PPT_c_2000-grData_remeas$PPT_c_norm # 99-00 = dry
PPT_c_2001_anom <- ppt.extr$PPT_c_2001-grData_remeas$PPT_c_norm # 00-01 = normal to slightly dry
PPT_c_2002_anom <- ppt.extr$PPT_c_2002-grData_remeas$PPT_c_norm # 01-02 = dry
PPT_c_2003_anom <- ppt.extr$PPT_c_2003-grData_remeas$PPT_c_norm # 02-03 = normal
PPT_c_2004_anom <- ppt.extr$PPT_c_2004-grData_remeas$PPT_c_norm # 03 = dry

PPT_fs_1999_anom <- ppt.extr$PPT_fs_1999-grData_remeas$PPT_fs_norm # 99 = wet
PPT_fs_2000_anom <- ppt.extr$PPT_fs_2000-grData_remeas$PPT_fs_norm # 00 = dry
PPT_fs_2001_anom <- ppt.extr$PPT_fs_2001-grData_remeas$PPT_fs_norm # 01 = dry to normal
PPT_fs_2002_anom <- ppt.extr$PPT_fs_2002-grData_remeas$PPT_fs_norm # 02 = dry
PPT_fs_2003_anom <- ppt.extr$PPT_fs_2003-grData_remeas$PPT_fs_norm # 03 = dry
PPT_fs_2004_anom <- ppt.extr$PPT_fs_2004-grData_remeas$PPT_fs_norm # 04 = normal

PPT_m_1999_anom <- ppt.extr$PPT_m_1999-grData_remeas$PPT_m_norm # 99 = wet
PPT_m_2000_anom <- ppt.extr$PPT_m_2000-grData_remeas$PPT_m_norm # 00 = dry
PPT_m_2001_anom <- ppt.extr$PPT_m_2001-grData_remeas$PPT_m_norm # 01 = dry to normal
PPT_m_2002_anom <- ppt.extr$PPT_m_2002-grData_remeas$PPT_m_norm # 02 = dry
PPT_m_2003_anom <- ppt.extr$PPT_m_2003-grData_remeas$PPT_m_norm # 03 = dry
PPT_m_2004_anom <- ppt.extr$PPT_m_2004-grData_remeas$PPT_m_norm # 04 = normal

# write csv with the PPT anomaly variables (for figure)
PPT.dr.anom <- cbind(PPT_pf_1999_anom, PPT_c_1999_anom, PPT_fs_1999_anom, PPT_m_1999_anom,
                     PPT_pf_2000_anom, PPT_c_2000_anom, PPT_fs_2000_anom, PPT_m_2000_anom,
                     PPT_pf_2001_anom, PPT_c_2001_anom, PPT_fs_2001_anom, PPT_m_2001_anom,
                     PPT_pf_2002_anom, PPT_c_2002_anom, PPT_fs_2002_anom, PPT_m_2002_anom,
                     PPT_pf_2003_anom, PPT_c_2003_anom, PPT_fs_2003_anom, PPT_m_2003_anom,
                     PPT_pf_2004_anom, PPT_c_2004_anom, PPT_fs_2004_anom, PPT_m_2004_anom)
write.csv(PPT.dr.anom, "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/PPTanomalies.csv", row.names = F)
PPT.normals <- cbind(grData_remeas$PPT_pf_norm, grData_remeas$PPT_c_norm, grData_remeas$PPT_fs_norm, grData_remeas$PPT_m_norm)
write.csv(PPT.dr.anom, "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/PPTnormals.csv", row.names = F)

# Create output data frame
output <- grData_remeas[, c("CN", "PREV_TRE_CN", "PLT_CN", "PREV_PLT_CN", "CONDID", 
                            "LAT", "LON", "ELEV",
                            "PREV_DRYBIO_AG", "DRYBIO_AG", "DRYBIO_AG_DIFF", #state variable
                            "PREVDIA", "DIA", "DIA_DIFF", #alternative state variable
                            "BAt1", "BAt2", "BA_DIFF", # another state variable
                            "STATUSCD", # live/dead
                            "MEASYEAR", "PREV_MEASYEAR", "CENSUS_INTERVAL", "BALIVE", 
                            "PPT_c", "T_c", "VPD_c", "PPTex_c", "Tex_c", "VPDex_c",
                            "PPT_wd", "T_wd", "VPD_wd", "PPTex_wd", "Tex_wd", "VPDex_wd",
                            "PPT_pf", "T_pf", "VPD_pf", "PPTex_pf", "Tex_pf", "VPDex_pf",
                            "PPT_fs", "T_fs", "VPD_fs", "PPTex_fs", "Tex_fs", "VPDex_fs",
                            "PPT_m", "T_m", "VPD_m", "PPTex_m", "Tex_m", "VPDex_m",
                            "PPT_yr", "T_yr", "VPD_yr", "PPTex_yr", "Tex_yr", "VPDex_yr",
                            "PPT_c_norm", "T_c_norm", "VPD_c_norm",
                            "PPT_wd_norm", "T_wd_norm", "VPD_wd_norm",
                            "PPT_pf_norm", "T_pf_norm", "VPD_pf_norm",
                            "PPT_fs_norm", "T_fs_norm", "VPD_fs_norm",
                            "PPT_m_norm", "T_m_norm", "VPD_m_norm",
                            "PPT_yr_norm", "T_yr_norm", "VPD_yr_norm",
                            "PPT_c_anom", "T_c_anom", "VPD_c_anom", "PPTex_c_anom", "Tex_c_anom", "VPDex_c_anom",
                            "PPT_wd_anom", "T_wd_anom", "VPD_wd_anom", "PPTex_wd_anom", "Tex_wd_anom", "VPDex_wd_anom",
                            "PPT_pf_anom", "T_pf_anom", "VPD_pf_anom", "PPTex_pf_anom", "Tex_pf_anom", "VPDex_pf_anom",
                            "PPT_fs_anom", "T_fs_anom", "VPD_fs_anom", "PPTex_fs_anom", "Tex_fs_anom", "VPDex_fs_anom",
                            "PPT_m_anom", "T_m_anom", "VPD_m_anom", "PPTex_m_anom", "Tex_m_anom", "VPDex_m_anom",
                            "PPT_yr_anom", "T_yr_anom", "VPD_yr_anom", "PPTex_yr_anom", "Tex_yr_anom", "VPDex_yr_anom",
                            "PPT_drought", "PPT_50mo_norm", "PPT_dr_anom",
                            "Tmean_drought", "T_50mo_norm", "T_dr_anom", "cum_T_anom",
                            "PPT_pf_dr", "PPT_c_dr", "PPT_fs_dr", "PPT_m_dr",
                            "PPT_pf_dr_anom", "PPT_c_dr_anom", "PPT_fs_dr_anom", "PPT_m_dr_anom")]
write.csv(output, "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Survival/SurvivalData.csv", row.names = F)


#### OLD CODE (FOR LOOP LOOK-UPS)
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


for (i in 1:nrow(grData_remeas)) {
  print(i)
  grData_remeas[i, "PPT_c"] <- rowMeans(ppt.extr[i, paste0("PPT_c_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
  grData_remeas[i, "VPD_c"] <- rowMeans(vpd.extr[i, paste0("VPD_c_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
  grData_remeas[i, "PPT_w"] <- rowMeans(ppt.extr[i, paste0("PPT_w_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
  grData_remeas[i, "VPD_w"] <- rowMeans(vpd.extr[i, paste0("VPD_w_", grData_remeas[i, "PREV_MEASYEAR"]:(grData_remeas[i, "MEASYEAR"]-1))])
}
