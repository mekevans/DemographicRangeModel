##Evaluate PIED integral projection model

library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(glmmTMB)

# Load vital rate and IPM functions
source("./Code/IPM/BuildIPM.R")

# define path
#path = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/"
#path = "C:/Users/mekevans/Documents/Cdrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PIED_IPM/MEKEvans/"

#PRISM.norm.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/Normals/"
PRISM.norm.path <-  "./ClimateData/PRISM/Normals/"

# load models and scaling --------------------------------------------------

# growth model + scaling
# from modelSelection_Growth.R
load("./Code/IPM/GrRescaling.Rdata")

# survival model + scaling
# from modelSelection_Survival.R
# load(paste0(path, "Code/IPM/SurvRescaling.Rdata"))
load("./Code/IPM/SurvRescalingNoFire.Rdata")
#load(paste0(path, "Code/IPM/SurvRescalingBA.Rdata"))

# recruitment model + scaling
# from modelSelection_Recruit.R
load("./Code/IPM/RecruitRescaling.Rdata")

# information on the size distribution of recruits (ingrowth)
# from dataPrepRecruitment.R
load("./Code/IPM/recrstats.rda")

# Load FIA survival, growth data
FIA <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)
# when running IPM without trees killed by fire, technically should filter out those trees
# prolly makes ~0 difference
FIA <- FIA[!(FIA$DSTRBCD1 %in% c(30, 31, 32, 80)), ]

xy = FIA[, c("LON", "LAT")]
spFIA = sp::SpatialPoints(xy)
spFIA = SpatialPointsDataFrame(spFIA, FIA)

## Test vital rate functions
# try out s.x on some realistic combinations of values
#s.x(size.x = c(8.4, 8.4, 8.4), Tann = c(2.0, 8.9, 15.0), PPTann = c(395, 395, 395), interval = 10)
# look at survival as a function of PPTann (ave-sized tree, ave Tann...)
data_test<-data.frame(PPT_yr=seq(1, 1000, length.out = 100),BALIVE=110,T_yr=8.9)
test = s.x(model = smodel.clim, size.x = 8.4, data = data_test, interval = 10)
plot(seq(1, 1000, length.out = 100), test)
# look at survival as a function of Tann
data_test<-data.frame(PPT_yr=395,BALIVE=110,T_yr=seq(-4, 30, length.out = 100))
test = s.x(model = smodel.clim, size.x = 8.4, data = data_test, interval = 10)
plot(seq(-4, 30, length.out = 100), test)
# look at survival as a function of balive, if balive is in mort model
data_test<-data.frame(PPT_yr=395,BALIVE=seq(0, 360, length.out = 100),T_yr=8.9)
test = s.x(model = sbase_int, size.x = 8.4, data = data_test, interval = 10)
plot(seq(0, 360, length.out = 100), test, ylab = "10-yr survival probability", xlab = "basal area live trees")
#abline(v = min(FIA$BALIVE)); abline(v = max(FIA$BALIVE)) # extrapolation lines
abline(v = 239) # 239 is the maximum value in Michiel's interpolated BALIVE; 360 is the largest value observed at PIED plots

# try out g.yx
#d_growth <- g.yx(y, 10, balive = 110, PPTann = 348, Tann = 17) # time step = 1 year
dx<-0.1
y=seq(min(FIA$PREVDIA,na.rm=T),(max(FIA$PREVDIA,na.rm=T)+5),dx)
data_test<-data.frame(PPT_yr=395,BALIVE=110,T_yr=8.9)
d_growth <- g.yx(model = gmodel.clim, growSD = growSD.clim, y, 39, data = data_test) # time step = 1 year
plot(y, d_growth, type = "l", ylab = "density", xlim = c(30,40))
sum(d_growth)*dx
# sum of d_growth should = 1.0 (but isn't)
# sum of d_growth*dx = 1 for small enough dx

# look at mean growth as a function of PPTann (ave-sized tree, ave Tann...)
data_test<-data.frame(BALIVE = 110, T_yr = 8.9, PPT_yr = seq(1, 1500, length.out = 100))
test = g.mean(model = gmodel.clim, size.x = 8.4, data = data_test, interval = 10)
plot(seq(1, 1000, length.out = 100), test, ylab = "10-yr growth incr (in)", xlab = "mean annual precipitation (mm)")
abline(v = min(FIA$PPT_yr_norm)); abline(v = max(FIA$PPT_yr_norm)) # extrapolation lines
# look at mean growth as a function of Tann
data_test<-data.frame(BALIVE = 110, T_yr = seq(-1, 24, length.out = 100), PPT_yr = 395)
test = g.mean(model = gmodel_int, size.x = 8.4, data = data_test, interval = 10)
plot(seq(-4, 30, length.out = 100), test, ylab = "10-yr growth incr (in)", xlab = "mean annual temperature (C)")
abline(v = min(FIA$T_yr_norm)); abline(v = max(FIA$T_yr_norm)) # extrapolation lines
# look at mean growth as a function of balive
data_test<-data.frame(BALIVE = seq(0, 240, length.out = 100), T_yr = 8.9, PPT_yr = 395)
test = g.mean(model = gmodel_int, size.x = 8.4, data = data_test, interval = 10)
plot(seq(0, 240, length.out = 100), test, ylab = "10-yr growth incr (in)", xlab = "basal area live trees")
abline(v = min(FIA$BALIVE)); abline(v = max(FIA$BALIVE)) # extrapolation lines
abline(v = 190, col = "blue", lty = 2, lwd = 3) # clamping line

# try out fec (y's are the 500 mid-point sizes)
data_test<-data.frame(BALIVE=110,PPT_yr=348,T_yr=17)
d_recruit <- fec(model = r_q, y, 10, data = data_test) # time step = 1 year
data_test<-data.frame(BALIVE=88.6,PPT_yr=244,T_yr=7.5)
d_recruit <- fec(model = r_balive_clim, y, 1.0581, data = data_test) # smallest tree
d_recruit <- fec(model = r_balive_clim, y, 35, data = data_test) # largest tree
# the probability density looks the same, no matter what the size of the "parent" tree
# because our model of recruitment doesn't account for the influence of tree size
plot(y, d_recruit, type = "l", ylab = "density")

# Load climate layers
# these created using the script "current.R"
ppt_yr_raster <- raster(paste0(PRISM.norm.path, "PPT_year.tif"))
t_yr_raster <- raster(paste0(PRISM.norm.path, "T_year.tif"))
t_wd_raster <- raster(paste0(PRISM.norm.path, "T_wd.tif"))
t_c_raster <- raster(paste0(PRISM.norm.path, "T_c.tif"))
t_m_raster <- raster(paste0(PRISM.norm.path, "T_m.tif"))
# stand-level basal area raster
ba_raster <- raster("./BA/BA.tif")
# ba_raster <- raster("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/BA/BA.tif")

ppt_yr_raster <- resample(ppt_yr_raster, ba_raster)
t_yr_raster <- resample(t_yr_raster, ba_raster)
t_wd_raster <- resample(t_wd_raster, ba_raster)
t_c_raster <- resample(t_c_raster, ba_raster)
t_m_raster <- resample(t_m_raster, ba_raster)

# aggregate for now, just to make it faster
#ppt_yr_raster <- aggregate(ppt_yr_raster, 2)
#t_yr_raster <- aggregate(t_yr_raster, 2)
#ba_raster <- aggregate(ba_raster, 2)


# Prepare empty rasters
lambda <- ppt_yr_raster
lambda <- setValues(lambda, NA)
growth <- ppt_yr_raster
growth <- setValues(growth, 0)
survival <- ppt_yr_raster
survival <- setValues(survival, 0)
reproduction <- ppt_yr_raster
reproduction <- setValues(reproduction, 0)

min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500
# I think max.size should be smaller...59.1 inches DRC is totally unrealistic

# Build IPMs and calculate lambda
for (i in 1:nrow(ppt_yr_raster)) {
  for (j in 1:ncol(ppt_yr_raster)) {
    # Extract climate for cell
    pred_data <- data.frame(PPT_yr=as.numeric(ppt_yr_raster[i,j]),
                            T_yr=as.numeric(t_yr_raster[i,j]),
                            BALIVE=as.numeric(ba_raster[i,j]),
                            T_wd_norm=as.numeric(t_wd_raster[i,j]),
                            T_c_norm=as.numeric(t_c_raster[i,j]),
                            T_m_norm=as.numeric(t_m_raster[i,j]))
    # Check for missing value
    if (is.na(pred_data$PPT_yr) | is.na(pred_data$T_yr) | is.na(pred_data$BALIVE)) {
      lambda[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    # Calculate lambda
    K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.clim, smodel=smodel.clim, 
               rmodel=rmodel.clim, gSD=growSD.clim,
               data=pred_data)
    lambda_val <- Re(eigen(K)$values[1])
    print(lambda_val)
    lambda[i,j] <- lambda_val
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(lambda, main = "Lambda"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

# Fill growth, survival, reproduction rasters
for (i in 1:nrow(ppt_yr_raster)) {
  for (j in 1:ncol(ppt_yr_raster)) {
    # Extract climate for cell
    pred_data <- data.frame(PPT_yr=as.numeric(ppt_yr_raster[i,j]),
                            T_yr=as.numeric(t_yr_raster[i,j]),
                            BALIVE=as.numeric(ba_raster[i,j]),
                            T_wd_norm=as.numeric(t_wd_raster[i,j]),
                            T_c_norm=as.numeric(t_c_raster[i,j]),
                            T_m_norm=as.numeric(t_m_raster[i,j]))
    # Check for missing value
    if (is.na(pred_data$PPT_yr) | is.na(pred_data$T_yr) | is.na(pred_data$BALIVE)) {
      lambda[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    
    G <- g.mean(size.x = median(FIA$DIA, na.rm = T), model = gmodel_q, data=pred_data) # interval = 1 in g.mean function
    # S <- s.x(size.x = median(FIA$DIA, na.rm = T), PPTann = ppt_yr_val, Tann = t_yr_val)
    S <- s.x(size.x = median(FIA$DIA, na.rm = T), model = sbase_q, data=pred_data)
    R <- f.mean(model = r_q, data=pred_data)
    growth[i,j] <- G
    survival[i,j] <- S
    reproduction[i,j] <- R
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(growth); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

pdf("./Output/PIED.q.pdf")
plot(lambda, main = "Lambda"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(growth, main = "Growth"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(survival, main = "Survival", zlim = c(0.95, 1)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(reproduction, main = "Reproduction"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/Climate_maps.pdf")
plot(ba_raster, main = "Live Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(ppt_yr_raster, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(t_yr_raster, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

# Crop rasters
fourState <- readOGR("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans", "4state")
fourState <- gUnaryUnion(fourState)
lambda <- mask(lambda, fourState)
growth <- mask(growth, fourState)
survival <- mask(survival, fourState)
reproduction <- mask(reproduction, fourState)

# Export
writeRaster(lambda, "./Output/BC/PIED.q_lambda.tif", overwrite = T)
writeRaster(growth, "./Output/BC/PIED.q_growth.tif", overwrite = T)
writeRaster(survival, "./Output/BC/PIED.q_survival.tif", overwrite = T)
writeRaster(reproduction, "./Output/BC/PIED.q_reproduction.tif", overwrite = T)

lambda<- readRaster("./Output/BC/PIED.int.q_lambda.tif")
