library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(glmmTMB)

# Load vital rate and IPM functions
source("./Code/IPM/BuildIPM.R")

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


# Perturb coefficient
#coefN <- 5
#oldCoef <- as.numeric(smodel$coefficients[coefN])
#smodel$coefficients[coefN] <- smodel$coefficients[coefN] * 1.01

# Set IPM parameters
min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500

# Set aggregation factor
#aggr <- 20

# Load climate layers
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

# Aggregate
#ppt_c_raster <- aggregate(ppt_c_raster, aggr)
#ppt_w_raster <- aggregate(ppt_w_raster, aggr)
#vpd_c_raster <- aggregate(vpd_c_raster, aggr)
#vpd_w_raster <- aggregate(vpd_w_raster, aggr)

# Crop BA raster
#fourState <- readOGR("D:/EvansLab/IPM/Data/BA", "4state")
#fourState <- gUnaryUnion(fourState)
#ba_raster <- mask(ba_raster, fourState)

# Prepare empty raster
lambda <- ppt_yr_raster
lambda <- setValues(lambda, NA)
#lambda <- aggregate(lambda, aggr)

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
    K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel_q, smodel=sbase_q, 
               rmodel=r_q, gSD=growSD_q,
               data=pred_data, rperturb = 0.01)
    lambda_val <- Re(eigen(K)$values[1])
    print(lambda_val)
    lambda[i,j] <- lambda_val
    }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(lambda)
}

# Export maps
unperturbed <- raster("./Output/BC/PIED.q_lambda.tif")
elasticity_recruit <- (lambda-unperturbed)/(unperturbed*0.01)
writeRaster(elasticity_growth, "./Code/Elasticities/elasticity_growth.tif", overwrite = T)
writeRaster(elasticity_survival, "./Code/Elasticities/elasticity_survival.tif", overwrite = T)
writeRaster(elasticity_recruit, "./Code/Elasticities/elasticity_recruit.tif", overwrite = T)

pdf("./Output/elasticities.pdf")
plot(elasticity_growth, main = "Growth"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(elasticity_survival, main = "Survival"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(elasticity_recruit, main = "Recruitment"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()
