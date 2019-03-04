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

#LTRE
## Collect dparam/dba for the ba-dependent parameters (these are the ba climate coefficients from Bayes models)
#dparam.denv<-c(gmodel_q@beta[4],gmodel_q@beta[5],gmodel_q@beta[6],gmodel_q@beta[7],gmodel_q@beta[8],gmodel_q@beta[9]) #growth
#dparam.denv<-c(sbase_q@beta[3],sbase_q@beta[7],sbase_q@beta[4],sbase_q@beta[8],sbase_q@beta[5],sbase_q@beta[9]) #survival
dparam.denv<-c(r_q$fit$par[2],r_q$fit$par[3],r_q$fit$par[4],r_q$fit$par[5],r_q$fit$par[6],r_q$fit$par[7]) #recruitment

## Collect dlambda/dparam for the climate-dependent parameters
## Sensitivities will vary with level of spei, so the output is a matrix with vital rates as rows and spei as columns

dlambda.dparam<-array(NA,dim=c(length(dparam.denv),nrow(ppt_yr_raster),ncol(ppt_yr_raster)))
perturb<-0.01
## here are the indices of the climate-dependent parameters
#dlambda.dparam.indices<-c(1,4,1,6,1,8) #growth
#dlambda.dparam.indices<-c(1,3,1,4,1,5) #survival
dlambda.dparam.indices<-c(1,2,1,4,1,6) #recruitment

cov<-c("ba","ba","ppt","ppt","t","t")
#1,5) #survival
#1,2) #recruitment
## Here, store the product of dlambda.dparam and dparam.ddens
LTRE.arr<-array(NA,dim=c(length(dparam.denv),nrow(ppt_yr_raster),ncol(ppt_yr_raster)))##6 rows for 6 env-dependent parameters
#growth_cont_ba <- ppt_yr_raster
#growth_cont_ba <- setValues(growth_cont_ba, NA)
#growth_cont_ppt <- growth_cont_ba
#growth_cont_ba <- growth_cont_ba
#surv_cont_ba <- ppt_yr_raster
#surv_cont_ba <- setValues(surv_cont_ba, NA)
#surv_cont_ppt <- surv_cont_ba
#surv_cont_ba <- surv_cont_ba
recr_cont_ba <- ppt_yr_raster
recr_cont_ba <- setValues(recr_cont_ba, NA)
recr_cont_ppt <- recr_cont_ba
recr_cont_ba <- recr_cont_ba
lambda<-recr_cont_ba

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
      #growth_cont_ba[i,j] <- NA
      #growth_cont_ppt[i,j] <- NA
      #growth_cont_t[i,j] <- NA
      #surv_cont_ba[i,j] <- NA
      #surv_cont_ppt[i,j] <- NA
      #surv_cont_t[i,j] <- NA
      recr_cont_ba[i,j] <- NA
      recr_cont_ppt[i,j] <- NA
      recr_cont_t[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    
    K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel_q, smodel=sbase_q, 
               rmodel=r_q, gSD=growSD_q,
               data=pred_data)
    lam<-Re(eigen(K)$values[1])
    lambda[i,j]<-lam
    
    pars<-seq(1:6)
    if (pred_data$BALIVE>190){
      pars<-pars[-c(1:2)]
      LTRE.arr[1:2,i,j]<-0
    }
    #if (pred_data$T_yr>190){
    #  pars<-pars[-5:6]
    #  LTRE.arr[5:6,i,j]<-0
    #}
    ## cycle over parameters; compute sensitivities and LTREs
    for(p in pars){
      #oldCoef <- as.numeric(gmodel_q@beta[dlambda.dparam.indices[p]])
      #gmodel_q@beta[dlambda.dparam.indices[p]] <- gmodel_q@beta[dlambda.dparam.indices[p]] + perturb
      #oldCoef <- as.numeric(sbase_q@beta[dlambda.dparam.indices[p]])
      #sbase_q@beta[dlambda.dparam.indices[p]] <- sbase_q@beta[dlambda.dparam.indices[p]] + perturb
      oldCoef <- as.numeric(r_q$fit$par[dlambda.dparam.indices[p]])
      r_q$fit$par[dlambda.dparam.indices[p]] <- r_q$fit$par[dlambda.dparam.indices[p]] + perturb
      K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel_q, smodel=sbase_q, 
                 rmodel=r_q, gSD=growSD_q,
                 data=pred_data)
      perturb.lam<-Re(eigen(K)$values[1])
      dlambda.dparam[p,i,j]<-(perturb.lam-lambda[i,j])/perturb
      LTRE.arr[p,i,j]<-dlambda.dparam[p,i,j]*dparam.denv[p]
      #gmodel_q@beta[dlambda.dparam.indices[p]] <- oldCoef
      #sbase_q@beta[dlambda.dparam.indices[p]] <- oldCoef
      r_q$fit$par[dlambda.dparam.indices[p]] <- oldCoef
    }
    # Calculate lambda
    #growth_cont_ba[i,j] <- sum(LTRE.arr[which(cov=="ba"),i,j])
    #growth_cont_ppt[i,j] <- sum(LTRE.arr[which(cov=="ppt"),i,j])
    #growth_cont_t[i,j] <- sum(LTRE.arr[which(cov=="t"),i,j])
    #print(growth_cont_ba[i,j])
    #surv_cont_ba[i,j] <- sum(LTRE.arr[which(cov=="ba"),i,j])
    #surv_cont_ppt[i,j] <- sum(LTRE.arr[which(cov=="ppt"),i,j])
    #surv_cont_t[i,j] <- sum(LTRE.arr[which(cov=="t"),i,j])
    #print(surv_cont_ba[i,j])
    recr_cont_ba[i,j] <- sum(LTRE.arr[which(cov=="ba"),i,j])
    recr_cont_ppt[i,j] <- sum(LTRE.arr[which(cov=="ppt"),i,j])
    recr_cont_t[i,j] <- sum(LTRE.arr[which(cov=="t"),i,j])
    print(recr_cont_ba[i,j])
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  #plot(growth_cont_ba)
  #plot(surv_cont_ba)
  plot(recr_cont_ba)
}

writeRaster(growth_cont_ba, "./Code/Elasticities/growth_cont_ba.tif", overwrite = T)
writeRaster(growth_cont_ppt, "./Code/Elasticities/growth_cont_ppt.tif", overwrite = T)
writeRaster(growth_cont_t, "./Code/Elasticities/growth_cont_t.tif", overwrite = T)

writeRaster(surv_cont_ba, "./Code/Elasticities/surv_cont_ba.tif", overwrite = T)
writeRaster(surv_cont_ppt, "./Code/Elasticities/surv_cont_ppt.tif", overwrite = T)
writeRaster(surv_cont_t, "./Code/Elasticities/surv_cont_t.tif", overwrite = T)

writeRaster(recr_cont_ba, "./Code/Elasticities/recr_cont_ba.tif", overwrite = T)
writeRaster(recr_cont_ppt, "./Code/Elasticities/recr_cont_ppt.tif", overwrite = T)
writeRaster(recr_cont_t, "./Code/Elasticities/recr_cont_t.tif", overwrite = T)

pdf("./Output/growth_cont.pdf")
plot(growth_cont_ba, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(growth_cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(growth_cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/surv_cont.pdf")
plot(surv_cont_ba, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(surv_cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(surv_cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/recr_cont.pdf")
plot(recr_cont_ba, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(recr_cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(recr_cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

ba_ltre<-sum(growth_cont_ba,surv_cont_ba,recr_cont_ba)
plot(ba_ltre, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
ppt_ltre<-sum(growth_cont_ppt,surv_cont_ppt,recr_cont_ppt)
plot(ppt_ltre, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
t_ltre<-sum(growth_cont_t,surv_cont_t,recr_cont_t)
plot(t_ltre, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)

#Calculate dlambda/dpred (the contributions of each vital rate should sum to these total changes)
perturb<-0.01
dlambda <- ppt_yr_raster
dlambda <- setValues(recr_cont, NA)
lambda <- ppt_yr_raster
lambda <- setValues(recr_cont, NA)

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
      dlambda[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    
    #K1<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel_q, smodel=sbase_q, 
    #           rmodel=r_q, gSD=growSD_q,
    #          data=pred_data)
    #lam<-Re(eigen(K1)$values[1])
    
    #lambda[i,j]<-lam
    ## perturb predictor; compute change in lambda
    
    pred_data$T_yr<-pred_data$T_yr+perturb
    K2<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel_q, smodel=sbase_q, 
                rmodel=r_q, gSD=growSD_q,
                data=pred_data)
    perturb.lam<-Re(eigen(K2)$values[1])
    
    dlambda[i,j]<-(perturb.lam-lambda[i,j])/perturb
    print(dlambda[i,j])
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(dlambda)
}

ba_ltre_total<-dlambda
ppt_ltre_total<-dlambda
t_ltre_total<-dlambda

pdf("./Output/ltre_totals.pdf")
plot(ba_ltre_total, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(ppt_ltre_total, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(t_ltre_total, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

writeRaster(ba_ltre_total, "./Code/Elasticities/ba_ltre_total.tif", overwrite = T)
writeRaster(ppt_ltre_total, "./Code/Elasticities/ppt_ltre_total.tif", overwrite = T)
writeRaster(t_ltre_total, "./Code/Elasticities/t_ltre_total.tif", overwrite = T)
