### create model of BALIVE
### to make a map of predicted BALIVE as a function of MAT and MAP
### for use in building IPM

library(sp)
library(raster)
library(dplyr)
library(effects)
library(mgcv)
library(glmmTMB)
library(DHARMa)

# path
# path = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/"

path = ""

# climate.path <- "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/"

climate.path = ""



### create data object with BALIVE and climate normals ------------------------------------------------------------------------------------

# conds <- read.csv(paste0(path,"FIAdata/COND_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
# # Remove non-forest conditions (missing BALIVE, coded as 0 so causes underestimation)
# # conds <- subset(conds, COND_STATUS_CD == 1) # "accessible forest land" by FIA classification
# 
# # try including nonforested land in model
# conds.1 <- subset(conds, COND_STATUS_CD == 1) # "accessible forest land" by FIA classification
# conds.2 <- subset(conds, COND_STATUS_CD == 2) # "nonforest land" by FIA classification
# conds.2 <- conds.2[conds.2$PRESNFCD %in% c(20, 45), ] # only two categories of non-forest land that are not caused by (human) land use
# conds.2$BALIVE <- 0 # populate with zeros
# 
# # join the two together
# conds <- rbind(conds.1, conds.2)
# 
# hist(conds$BALIVE) # almost all of the data has values <=400; outliers above 400
# summary(conds$BALIVE) # maximum value is 1363!
# plot(conds$CONDPROP_UNADJ, conds$BALIVE) # reveals that one outlier for BALIVE is a cond with very low CONDPROP_UNADJ
# conds <- subset(conds, BALIVE < 1000) # removes one outlier!
# # another data point that perhaps should be removed is #2919 (see resids vs. leverage)
# 
# # Read in plot data and get coordinates and previous measurement year
# plots <- read.csv(paste0(path,"FIAdata/PLOT_COMBINED.csv",sep=''), header = T, stringsAsFactors = F)
# 
# conds$LAT <- plots$LAT[match(conds$PLT_CN, plots$CN)]
# conds$LON <- plots$LON[match(conds$PLT_CN, plots$CN)]
# conds$ELEV <- plots$ELEV[match(conds$PLT_CN, plots$CN)]
# conds$MEASYEAR <- plots$MEASYEAR[match(conds$PLT_CN, plots$CN)]
# 
# # Make data spatial
# CondsSpat <- SpatialPointsDataFrame(coords = cbind(conds$LON, conds$LAT), 
#                                  data = conds, 
#                                  proj4string = CRS("+proj=longlat +datum=NAD83"))
# 
# # import PRISM normals
# PPT.norm <- stack(paste0(climate.path,"pptNormals.tif",sep=''))
# TMP.norm <- stack(paste0(climate.path,"tmpNormals.tif",sep=''))
# VPD.norm <- stack(paste0(climate.path,"vpdNormals.tif",sep=''))
# 
# # raster::extract PRISM normals 1981-2010
# ppt.norm.extr <- raster::extract(PPT.norm, CondsSpat)
# tmp.norm.extr <- raster::extract(TMP.norm, CondsSpat)
# vpd.norm.extr <- raster::extract(VPD.norm, CondsSpat)
# 
# ppt.norm.extr <- as.data.frame(ppt.norm.extr) # not sure this is necessary?
# tmp.norm.extr <- as.data.frame(tmp.norm.extr)
# vpd.norm.extr <- as.data.frame(vpd.norm.extr)
# 
# # reasonable column names (not actually necessary)
# colnames(ppt.norm.extr) <- paste0("ppt_", 1:12) 
# colnames(tmp.norm.extr) <- paste0("tmp_", 1:12)
# colnames(vpd.norm.extr) <- paste0("vpd_", 1:12)
# 
# # make seasonal normals and add to conditions data frame
# # cool season = Nov-Mar
# conds$PPT_c_norm <- rowSums(ppt.norm.extr[, c(1:3, 11:12)])
# conds$T_c_norm <- rowMeans(tmp.norm.extr[, c(1:3, 11:12)])
# conds$VPD_c_norm <- rowMeans(vpd.norm.extr[, c(1:3, 11:12)])
# # previous fall = pSept-pOct
# conds$PPT_pf_norm <- rowSums(ppt.norm.extr[, c(9:10)])
# conds$T_pf_norm <- rowMeans(tmp.norm.extr[, c(9:10)])
# conds$VPD_pf_norm <- rowMeans(vpd.norm.extr[, c(9:10)])
# # foresummer = Apr-Jun
# conds$PPT_fs_norm <- rowSums(ppt.norm.extr[, c(4:6)])
# conds$T_fs_norm <- rowMeans(tmp.norm.extr[, c(4:6)])
# conds$VPD_fs_norm <- rowMeans(vpd.norm.extr[, c(4:6)])
# # warm, dry months = Apr-Jun + Sept-Oct
# conds$PPT_wd_norm <- rowSums(ppt.norm.extr[, c(4:6, 9:10)])
# conds$T_wd_norm <- rowMeans(tmp.norm.extr[, c(4:6, 9:10)])
# conds$VPD_wd_norm <- rowMeans(vpd.norm.extr[, c(4:6, 9:10)])
# # monsoon = Jul-Aug
# conds$PPT_m_norm <- rowSums(ppt.norm.extr[, c(7:8)])
# conds$T_m_norm <- rowMeans(tmp.norm.extr[, c(7:8)])
# conds$VPD_m_norm <- rowMeans(vpd.norm.extr[, c(7:8)])
# # water year
# conds$PPT_yr_norm <- rowSums(ppt.norm.extr[, c(1:12)])
# conds$T_yr_norm <- rowMeans(tmp.norm.extr[, c(1:12)])
# conds$VPD_yr_norm <- rowMeans(vpd.norm.extr[, c(1:12)])
# 
# # Create output data frame
# output <- conds[, c("PLT_CN", "CONDID", "CONDPROP_UNADJ",
#                             "LAT", "LON", "ELEV",
#                             "BALIVE", 
#                             "PPT_c_norm", "T_c_norm", "VPD_c_norm",
#                             "PPT_wd_norm", "T_wd_norm", "VPD_wd_norm",
#                             "PPT_pf_norm", "T_pf_norm", "VPD_pf_norm",
#                             "PPT_fs_norm", "T_fs_norm", "VPD_fs_norm",
#                             "PPT_m_norm", "T_m_norm", "VPD_m_norm",
#                             "PPT_yr_norm", "T_yr_norm", "VPD_yr_norm")]
# # write.csv(output, paste0(path, "BA/BALIVEdata.csv"), row.names = F)
# write.csv(output, paste0(path, "BA/BALIVEdata2.csv"), row.names = F)
# 



# load data object for model of BALIVE ---------------------------------------------------------------
BAdata <- read.csv(paste0(path, "BA/BALIVEdata2.csv"), header = T, stringsAsFactors = F)
BAdataold <- read.csv(paste0(path, "BA/BALIVEdata.csv"), header = T, stringsAsFactors = F)

# standardize covariates
BAdata.scaled <- BAdata %>% mutate_at(scale, .vars = vars(-PLT_CN, 
                                                          -CONDPROP_UNADJ, 
                                                          -CONDID,
                                                          -BALIVE,
                                                          -LAT, -LON))

# make scaling objects ------------------------------------------------------------------------


ba.predictors <- c("T_yr_norm", "PPT_yr_norm")

get_scale = function(data, predictors) {
  sc = list("scale" = NULL, "center"  = NULL)
  for (i in predictors) {
    sc$scale[i] = attributes(data[, i])$"scaled:scale"
    sc$center[i] = attributes(data[, i])$"scaled:center"
  }
  return(sc)
}

ba.scaling = get_scale(BAdata.scaled, ba.predictors)

# remove scaling information from the dataset so that the model doesnt expect scaled data in predict()
for (i in ba.predictors) {
  attributes(BAdata.scaled[, i]) = NULL
}


# simple models -----------------------------------------------------------

# # lm with sqrt transformation
# BAmodel1 <- lm(sqrt(BALIVE) ~ (PPT_yr_norm + T_yr_norm)^2 + 
#                   I(T_yr_norm^2) + I(PPT_yr_norm^2), 
#                   data = BAdata.scaled)
# plot(BAmodel1)
# plot(allEffects(BAmodel1))
# # res = resid(BAmodel1)
# # plot(density(res))
# plot(qqnorm(res))
# qqline(res)
# #plot(conds$BALIVE, res, ylab="residuals", xlab="BALIVE"); abline(0, 0)
# 
# 
# # lm with sqrt and unit transformation
# BAmodel1 <- lm(sqrt(BALIVE/10) ~ (PPT_yr_norm + T_yr_norm)^2 + 
#                  I(T_yr_norm^2) + I(PPT_yr_norm^2), 
#                data = BAdata.scaled)
# plot(BAmodel1)
# plot(allEffects(BAmodel1))
# # res = resid(BAmodel1)
# # plot(density(res))
# plot(qqnorm(res))
# qqline(res)
# #plot(conds$BALIVE, res, ylab="residuals", xlab="BALIVE"); abline(0, 0)
# 
# 
# 
# # gamma model
# BAmodel.g <- glm((BALIVE+1) ~ (PPT_yr_norm + T_yr_norm)^2 + 
#                 I(T_yr_norm^2) + I(PPT_yr_norm^2), 
#                  family=Gamma(link="log"),
#                  data = BAdata.scaled)
# plot(BAmodel.g)
# plot(allEffects(BAmodel.g))
# #res = simulateResiduals(BAmodel1)
# #plot(res)
# #plotResiduals(BAdata.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm") #
# #plotResiduals(BAdata.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm") # 



# gam model ---------------------------------------------------------------


BAmodel.gam <- gam(sqrt(BALIVE) ~ s(PPT_yr_norm, T_yr_norm) 
                   # + s(T_yr_norm)
                   + s(LON, LAT, k = 200),
                   # family=nb,
                   # family = Tweedie(p = 1.5, link = "identity"),
                   data = BAdata.scaled)
summary(BAmodel.gam)
mgcv::plot.gam(BAmodel.gam, pages = 1, scheme = 2, n2 = 120)
gam.check(BAmodel.gam)


library(mgcViz)
b <- getViz(BAmodel.gam)
# qq(b, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(b, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), 
   a.replin = list("alpha" = 0.2))

smoothScatter(predict(BAmodel.gam)^2, BAdata.scaled$BALIVE)
abline(c(0, 0), c(1, 1), col = "gold")





# glmmTMB -----------------------------------------------------------------


BAmodel.glmmTMB <- glmmTMB(sqrt(BALIVE) ~ PPT_yr_norm * T_yr_norm
                       + I(T_yr_norm^2)
                       + I(PPT_yr_norm^2),
                       ziformula = ~ 1 + PPT_yr_norm * T_yr_norm,
                       family=gaussian,
                       data = BAdata.scaled)
summary(BAmodel.glmmTMB)
plot(simulateResiduals(BAmodel.glmmTMB, n = 500))
plot(allEffects(BAmodel.glmmTMB))





# # gamlss model ---------------------------------------------------------------
# 
# library(gamlss)
# library(gamlss.add)
# BAmodel.gam <- gamlss::gamlss(sqrt(BALIVE) ~ ga(~ s(PPT_yr_norm)
#                                             + s(T_yr_norm)
#                                             + s(LON, LAT)),
#                               family = gaussian,
#                               sigma.formula = ~ 1,
#                               data = na.omit(BAdata.scaled))
# rqres.plot(BAmodel.gam)
# summary(BAmodel.gam)
# plot(BAmodel.gam)
# term.plot(BAmodel.gam, pages = 1)
# 
# 
# smoothScatter(predict(BAmodel.gam, type = "response")-1, BAdata.scaled$BALIVE)
# abline(c(0, 0), c(1, 1), col = "gold")






# export model for coefficients and scaling information -------------------


save(BAmodel1, ba.scaling, file = paste0(path, "BA/BArescaling.Rdata"))
save(BAmodel.g, ba.scaling, file = paste0(path, "BA/BArescalingGamma.Rdata"))




# prepare to make map --------------------------------------------------------------

# Load climate layers
# these created using the script "current.R"
ppt_yr_raster <- raster(paste0(climate.path, "ClimateData/PPT_year.tif"))
t_yr_raster <- raster(paste0(climate.path, "ClimateData/T_year.tif"))

# bring in scaling
load(paste0(path, "BA/BArescaling.Rdata"))




# function to predict balive using (scaled) climate data -----------------------------------
ba.pred <- function(PPTann, Tann, LON, LAT) { 
  
  # raw (unscaled) data
  badata <- data.frame(PPT_yr_norm = PPTann,
                      T_yr_norm = Tann)
  
  # rescaled data
  scaled.badata = data.frame(scale(badata, 
                                  scale = ba.scaling$scale[match(names(badata), names(ba.scaling$scale))], 
                                  center = ba.scaling$center[match(names(badata), names(ba.scaling$center))]))  
  scaled.badata$LON = LON
  scaled.badata$LAT = LAT
  
  # apply the fitted model to the scaled input data
  # bapred <- predict(BAmodel1, newdata = scaled.badata) 
  # bapred <- predict(BAmodel.gam, newdata = scaled.badata)
  bapred <- predict(BAmodel.glmmTMB, newdata = scaled.badata, type = "response")
  return(bapred^2) 

}


# use model and MAT, MAP to predict balive (make map) -----------------------------------

# Extract climate for cell
ppt_yr_val <- raster::getValues(ppt_yr_raster)
t_yr_val <- raster::getValues(t_yr_raster)
LON = raster::xFromCell(ppt_yr_raster, 1:ncell(ppt_yr_raster))
LAT = raster::yFromCell(ppt_yr_raster, 1:ncell(ppt_yr_raster))

# use model to predict balive
BA <- ba.pred(PPTann = ppt_yr_val, Tann = t_yr_val, 
              LON = LON, LAT = LAT)

# Prepare empty raster
balive <- ppt_yr_raster
balive <- setValues(balive, as.vector(BA))


par(mfrow = c(1, 1))
plot(balive)
points(LAT ~ LON, BAdataold, pch = 16, cex= 0.01, col = rgb(0,0,0,0.1))
hist(balive)






# explore observed BA in space ----------------------------------------------------

plot(LAT ~ LON, BAdata, pch = 16, cex = 0.1,
     col = rev(terrain.colors(10))[cut(log1p(BAdata$BALIVE), 10, ordered_result = T)], 
     asp = 1)

hist(BAdata$BALIVE)
hist(getValues(balive))



# random forest -----------------------------------------------------------



library(ranger)
library(GSIF)
library(rgdal)
library(sp)

BAdata.sp = BAdata[sample(1:nrow(BAdata), ceiling(0.05*nrow(BAdata))), ]
coordinates(BAdata.sp) <- ~ LON + LAT
proj4string(BAdata.sp) <- crs(ppt_yr_raster)


dat.grid = as(stack(ppt_yr_raster, t_yr_raster), "SpatialPixelsDataFrame")
dat.grid$LON = dat.grid@coords[, 1]
dat.grid$LAT = dat.grid@coords[, 2]


# with buffers

# grid.dist0 <- GSIF::buffer.dist(BAdata.sp["BALIVE"], 
#                                 predictionDomain = dat.grid[], 
#                                 classes = as.factor(1:nrow(BAdata.sp)))
# # summary(grid.dist0)
# 
# dn0 <- paste(names(grid.dist0), collapse="+")
# fm0 <- as.formula(paste("BALIVE ~ ", dn0))
# 
# ov.ba <- over(BAdata.sp["BALIVE"], grid.dist0)
# rm.ba <- cbind(BAdata.sp@data["BALIVE"], ov.ba)
# 
# m.ba <- ranger(fm0, rm.ba, quantreg=TRUE, num.trees=150, seed=1)
# m.ba
# 
# ba.rfd <- predict(m.ba, grid.dist0@data, type="quantiles")$predictions
# str(ba.rfd)
# 
# dat.grid$ba_rfd = ba.rfd[,2]
# dat.grid$ba_rfd_range = (ba.rfd[,3]-ba.rfd[,1])/2
# 
# sp::plot(dat.grid["ba_rfd"])
# 
# ba_rf_raster = as(dat.grid["ba_rfd"], "RasterLayer")
# plot(ba_rf_raster)
# 
# 
# 
# 
# # with covariates
# 
# grids.spc = GSIF::spc(dat.grid, as.formula("~ PPT_year + T_year"))
# 
# fm1 <- as.formula(paste("BALIVE ~ ", dn0, " + ", paste(names(grids.spc@predicted), collapse = "+")))
# fm1
# 
# ov.ba1 <- over(BAdata.sp["BALIVE"], grids.spc@predicted)
# rm.ba1 <- do.call(cbind, list(BAdata.sp@data["BALIVE"], ov.ba, ov.ba1))
# 
# m1.ba <- ranger(fm1, rm.ba1, importance="impurity", quantreg=TRUE, num.trees=150, seed=1)
# m1.ba
# 
# barplot(sort(ranger::importance(m1.ba), decreasing = T)[1:10])
# 
# ba.rfd1 <- predict(m1.ba, cbind(grids.spc@predicted@data, grid.dist0@data), type="quantiles")$predictions
# # str(ba.rfd1)
# 
# dat.grid$ba_rfd1 = ba.rfd1[,2]
# dat.grid$ba_rfd1_range = (ba.rfd1[,3]-ba.rfd1[,1])/2
# 
# sp::plot(dat.grid["ba_rfd1"])
# 
# ba_rf_raster = as(dat.grid["ba_rfd1"], "RasterLayer")
# plot(ba_rf_raster)



# with covariates only (PCA)

BAdata.sp = BAdata
coordinates(BAdata.sp) <- ~ LON + LAT
proj4string(BAdata.sp) <- crs(ppt_yr_raster)


grids.spc = GSIF::spc(dat.grid, as.formula("~ PPT_year + T_year"))
grids.spc@predicted$LON = dat.grid@coords[, 1]
grids.spc@predicted$LAT = dat.grid@coords[, 2]

fm2 <- as.formula(paste("BALIVE ~ ", " + ", paste(names(grids.spc@predicted), collapse = "+")))
fm2

ov.ba2 <- over(BAdata.sp["BALIVE"], grids.spc@predicted)
rm.ba2 <- do.call(cbind, list(BAdata.sp@data["BALIVE"], ov.ba2))

m2.ba <- ranger(fm2, na.omit(rm.ba2), importance="impurity", quantreg=TRUE, num.trees=300, seed=1)
m2.ba

# # with validation in prediction 
# sel = sample(1:nrow(na.omit(rm.ba2)), ceiling(0.5*nrow(na.omit(rm.ba2))))
# m2.ba <- ranger(fm2, na.omit(rm.ba2), importance="impurity", quantreg=TRUE, num.trees=300, seed=1, 
#                 case.weights = ifelse(1:nrow(na.omit(rm.ba2)) %in% sel, 1, 0), holdout = F)
# m2.ba <- ranger(fm2, na.omit(rm.ba2)[sel, ], importance="impurity", quantreg=TRUE, num.trees=300, seed=1)
#                 # case.weights = ifelse(1:nrow(na.omit(rm.ba2)) %in% sel, 1, 0), holdout = T)
# m2.ba


# barplot(sort(ranger::importance(m2.ba), decreasing = T)[1:10])

ba.rfd2 <- predict(m2.ba, cbind(grids.spc@predicted@data), type="quantiles")$predictions
# str(ba.rfd2)

dat.grid$ba_rfd2 = ba.rfd2[,2]
dat.grid$ba_rfd2_range = (ba.rfd2[,3]-ba.rfd2[,1])/2

# sp::plot(dat.grid["ba_rfd2"])

ba_rf_raster = as(dat.grid["ba_rfd2"], "RasterLayer")
plot(ba_rf_raster)

hist(getValues(ba_rf_raster))
hist(BAdata.scaled$BALIVE)





# with covariates only (original)


fm2 <- as.formula(paste("BALIVE ~ ", " + ", paste(names(dat.grid@data)[1:4], collapse = "+")))
fm2

ov.ba2 <- over(BAdata.sp["BALIVE"], dat.grid)
rm.ba2 <- do.call(cbind, list(BAdata.sp@data["BALIVE"], ov.ba2))

m2.ba <- ranger(fm2, na.omit(rm.ba2), importance="impurity", quantreg=TRUE, num.trees=150, seed=1)
m2.ba
# barplot(sort(ranger::importance(m2.ba), decreasing = T)[1:10])

ba.rfd2 <- predict(m2.ba, cbind(dat.grid@data), type="quantiles")$predictions
# str(ba.rfd2)

dat.grid$ba_rfd2 = ba.rfd2[,2]
dat.grid$ba_rfd2_range = (ba.rfd2[,3]-ba.rfd2[,1])/2

# sp::plot(dat.grid["ba_rfd2"])

ba_rf_raster = as(dat.grid["ba_rfd2"], "RasterLayer")
plot(ba_rf_raster)



# RF with tuning ----------------------------------------------------------

dat.grid@data)






# # Export ----------------------------------------------------------------------
# as pdf
pdf(paste0(path, "BA/BAmap.pdf"))
plot(balive, main = "BALIVE")
dev.off()

# as raster (for use in building IPM)
writeRaster(balive, paste0(path, "BA/balive.tif"), overwrite = T)
