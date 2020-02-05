library(ggplot2)
library(cowplot)
library(glmmTMB)
library(raster)

# Load IPM functions
source("./Code/IPM/BuildIPM.R")

PRISM.norm.path <-  "./ClimateData/PRISM/Normals/"

# Load climate layers
# these created using the script "current.R"
ppt_yr_raster <- raster(paste0(PRISM.norm.path, "PPT_year.tif"))
t_yr_raster <- raster(paste0(PRISM.norm.path, "T_year.tif"))
# stand-level basal area raster
ba_raster <- raster("./BA/balive_RF.tif")

# load models and scaling --------------------------------------------------

# growth model + scaling
# from modelSelection_Growth.R
load("./Code/IPM/GrRescaling.Rdata")

# survival model + scaling
# from modelSelection_Survival.R
# load(paste0(path, "Code/IPM/SurvRescaling.Rdata"))
load("./Code/IPM/SurvRescaling.Rdata")
#load(paste0(path, "Code/IPM/SurvRescalingBA.Rdata"))

# recruitment model + scaling
# from modelSelection_Recruit.R
load("./Code/IPM/RecruitRescaling.Rdata")

# information on the size distribution of recruits (ingrowth)
# from dataPrepRecruitment.R
load("./Code/IPM/recrstats.rda")

# Alternative models - gams
load("./Code/IPM/GrRescaling_gam.Rdata")
load("./Code/IPM/SurvRescaling_gam.Rdata")
load("./Code/IPM/RecruitRescaling_gam.Rdata")
load("./Code/IPM/recrstats.rda")

# Load FIA survival, growth data
FIA <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)
# when running IPM without trees killed by fire, technically should filter out those trees
# prolly makes ~0 difference
FIA <- FIA[!(FIA$DSTRBCD1 %in% c(30, 31, 32, 80)), ]

min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500

mins<-c(cellStats(ba_raster,stat="min"),cellStats(ppt_yr_raster,stat="min"),
        cellStats(t_yr_raster,stat="min"))
maxs<-c(cellStats(ba_raster,stat="max"),cellStats(ppt_yr_raster,stat="max"),
        cellStats(t_yr_raster,stat="max"))

# Prepare prediction data frame
noPoints <- 50
predictorNames <- c("BALIVE", "PPT_yr_norm", "T_yr_norm") #, "T_wd_norm", "T_c_norm", "T_m_norm")
predictors <- FIA[, predictorNames]
predictorDFs_clim <- list()
for (i in 1:ncol(predictors)) {
  currentPredictor <- predictorNames[i]
  predictorVals <- data.frame(BALIVE = rep(median(FIA$BALIVE), noPoints), 
                              PPT_yr = rep(median(FIA$PPT_yr_norm), noPoints), 
                              T_yr = rep(median(FIA$T_yr_norm), noPoints), 
                              #T_wd_norm = rep(median(FIA$T_wd_norm), noPoints), 
                              #T_c_norm = rep(median(FIA$T_c_norm), noPoints), 
                              #T_m_norm = rep(median(FIA$T_m_norm), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(mins[i], maxs[i], length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- Re(eigen(ipm_fun(min=min.size, max=max.size, n=n_dim, 
                                                   gmodel=gmodel.clim.gam, 
                                                   smodel=smodel.clim.gam, 
                                                   rmodel=rmodel.clim.gam, 
                                                   gSD=growSD.clim.gam,
                                                   data=predictorVals[j,],
                                                   s.t.clamp=F,g.t.clamp=F,g.ba.clamp=F,r.ba.clamp=F))$values[1])
  }
  predictorDFs_clim[[i]] <- predictorVals
}

predictorNames <- c("BALIVE", "PPT_yr_norm", "T_yr_norm") #, "T_wd_norm", "T_c_norm", "T_m_norm")
predictors <- FIA[, predictorNames]
predictorDFs_climcomp <- list()
for (i in 1:ncol(predictors)) {
  currentPredictor <- predictorNames[i]
  predictorVals <- data.frame(BALIVE = rep(median(FIA$BALIVE), noPoints), 
                              PPT_yr = rep(median(FIA$PPT_yr_norm), noPoints), 
                              T_yr = rep(median(FIA$T_yr_norm), noPoints), 
                              #T_wd_norm = rep(median(FIA$T_wd_norm), noPoints), 
                              #T_c_norm = rep(median(FIA$T_c_norm), noPoints), 
                              #T_m_norm = rep(median(FIA$T_m_norm), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(mins[i], maxs[i], length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- Re(eigen(ipm_fun(min=min.size, max=max.size, n=n_dim, 
                                                   gmodel=gmodel.clim.comp.gam, 
                                                   smodel=smodel.clim.comp.gam, 
                                                   rmodel=rmodel.clim.comp.gam,
                                                   gSD=growSD.clim.comp.gam,
                                                   data=predictorVals[j,],
                                                   s.t.clamp=F,g.t.clamp=F,g.ba.clamp=F,r.ba.clamp=F))$values[1])
  }
  predictorDFs_climcomp[[i]] <- predictorVals
}

predictorNames <- c("BALIVE", "PPT_yr_norm", "T_yr_norm") #, "T_wd_norm", "T_c_norm", "T_m_norm")
predictors <- FIA[, predictorNames]
predictorDFs_climcompfire <- list()
for (i in 1:ncol(predictors)) {
  currentPredictor <- predictorNames[i]
  predictorVals <- data.frame(BALIVE = rep(median(FIA$BALIVE), noPoints), 
                              PPT_yr = rep(median(FIA$PPT_yr_norm), noPoints), 
                              T_yr = rep(median(FIA$T_yr_norm), noPoints), 
                              #T_wd_norm = rep(median(FIA$T_wd_norm), noPoints), 
                              #T_c_norm = rep(median(FIA$T_c_norm), noPoints), 
                              #T_m_norm = rep(median(FIA$T_m_norm), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(mins[i], maxs[i], length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- Re(eigen(ipm_fun(min=min.size, max=max.size, n=n_dim, 
                                                   gmodel=gmodel.clim.comp, 
                                                   smodel=smodel.clim.comp.fire.lin, 
                                                   rmodel=rmodel.clim.comp.lin,
                                                   gSD=growSD.clim.comp,
                                                   data=predictorVals[j,],s.t.clamp=F))$values[1])
  }
  predictorDFs_climcompfire[[i]] <- predictorVals
}

predictorNames <- c("BALIVE", "PPT_yr_norm", "T_yr_norm") #, "T_wd_norm", "T_c_norm", "T_m_norm")
predictors <- FIA[, predictorNames]
predictorDFs_int <- list()
for (i in 1:ncol(predictors)) {
  currentPredictor <- predictorNames[i]
  predictorVals <- data.frame(BALIVE = rep(median(FIA$BALIVE), noPoints), 
                              PPT_yr = rep(median(FIA$PPT_yr_norm), noPoints), 
                              T_yr = rep(median(FIA$T_yr_norm), noPoints), 
                              #T_wd_norm = rep(median(FIA$T_wd_norm), noPoints), 
                              #T_c_norm = rep(median(FIA$T_c_norm), noPoints), 
                              #T_m_norm = rep(median(FIA$T_m_norm), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(mins[i], maxs[i], length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- Re(eigen(ipm_fun(min=min.size, max=max.size, n=n_dim, 
                                                   gmodel=gmodel.int.gam, 
                                                   smodel=smodel.int.gam, 
                                                   rmodel=rmodel.int.gam, 
                                                   gSD=growSD.int.gam,
                                                   data=predictorVals[j,],
                                                   s.t.clamp=F,g.t.clamp=F,g.ba.clamp=F,r.ba.clamp=F))$values[1])
  }
  predictorDFs_int[[i]] <- predictorVals
}

save(predictorDFs_clim, predictorDFs_climcomp, 
     predictorDFs_int, file="./Output/lambda_effects_gam.rda")

