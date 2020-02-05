library(raster)
#library(rgdal)
#library(rgeos)
#library(dplyr)
library(glmmTMB)
library(mgcv)
library(tidyverse)
#library(brms)

# Load vital rate and IPM functions

PRISM.norm.path <-  "./ClimateData/PRISM/Normals/"

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
#load("./Code/IPM/RecruitRescaling8_gam.Rdata")

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
ba_raster <- raster("./BA/balive_RF.tif")

ppt_yr_raster <- resample(ppt_yr_raster, ba_raster)
t_yr_raster <- resample(t_yr_raster, ba_raster)
t_wd_raster <- resample(t_wd_raster, ba_raster)
t_c_raster <- resample(t_c_raster, ba_raster)
t_m_raster <- resample(t_m_raster, ba_raster)

#Load rasters of unperturbed lambda values
lambda_clim<-raster("./Output/tifs/PIED.climclamp_lambda.tif")
lambda_climcomp<-raster("./Output/tifs/PIED.climcomp_lambda.tif")
lambda_int<-raster("./Output/tifs/PIED.int_lambda.tif")

#LTRE
## Collect dparam/dba for the env-dependent parameters (these are the coefficients vital rate models)
# Growth
dparam.denv.g.clim<-c(gmodel.clim@beta[3],gmodel.clim@beta[6],gmodel.clim@beta[4],gmodel.clim@beta[7])
dparam.denv.g.climcomp<-c(gmodel.clim.comp@beta[3],gmodel.clim.comp@beta[7],gmodel.clim.comp@beta[4],gmodel.clim.comp@beta[8],gmodel.clim.comp@beta[5],gmodel.clim.comp@beta[9]) #growth
dparam.denv.g.int<-c(gmodel.int@beta[3],gmodel.int@beta[7],gmodel.int@beta[10],gmodel.int@beta[13],gmodel.int@beta[14],
                     gmodel.int@beta[4],gmodel.int@beta[8],gmodel.int@beta[11],gmodel.int@beta[13],gmodel.int@beta[15],
                     gmodel.int@beta[5],gmodel.int@beta[9],gmodel.int@beta[12],gmodel.int@beta[14],gmodel.int@beta[15]) 
# Survival
dparam.denv.s.clim<-c(smodel.clim@beta[3],smodel.clim@beta[6],smodel.clim@beta[4],smodel.clim@beta[7])
dparam.denv.s.climcomp<-c(smodel.clim.comp@beta[3],smodel.clim.comp@beta[7],smodel.clim.comp@beta[4],smodel.clim.comp@beta[8],smodel.clim.comp@beta[5],smodel.clim.comp@beta[9]) #survival
dparam.denv.s.int<-c(smodel.int@beta[3],smodel.int@beta[7],smodel.int@beta[10],smodel.int@beta[13],smodel.int@beta[14],
                     smodel.int@beta[4],smodel.int@beta[8],smodel.int@beta[11],smodel.int@beta[13],smodel.int@beta[15],
                     smodel.int@beta[5],smodel.int@beta[9],smodel.int@beta[12],smodel.int@beta[14],smodel.int@beta[15]) 
# Recruitment
dparam.denv.r.clim<-c(rmodel.clim$fit$par[2],rmodel.clim$fit$par[4],rmodel.clim$fit$par[3],rmodel.clim$fit$par[5])
dparam.denv.r.climcomp<-c(rmodel.clim.comp$fit$par[2],rmodel.clim.comp$fit$par[5],rmodel.clim.comp$fit$par[3],rmodel.clim.comp$fit$par[6],rmodel.clim.comp$fit$par[4],rmodel.clim.comp$fit$par[7]) #recruitment

## here are the indices of the climate-dependent parameters
# Growth
dlambda.dparam.indices.g.clim<-c(1,3,1,4)
dlambda.dparam.indices.g.climcomp<-c(1,3,1,4,1,5)
dlambda.dparam.indices.g.int<-c(1,3,2,4,5,1,4,2,3,5,1,5,2,3,4) 
#survival
dlambda.dparam.indices.s.clim<-c(1,3,1,4)
dlambda.dparam.indices.s.climcomp<-c(1,3,1,4,1,5) 
dlambda.dparam.indices.s.int<-c(1,3,2,4,5,1,4,2,3,5,1,5,2,3,4) 
# Recruitment
dlambda.dparam.indices.r.clim<-c(1,2,1,3)
dlambda.dparam.indices.r.climcomp<-c(1,2,1,3,1,4) 

# Load IPM and LTRE functions
source("./Code/IPM/BuildIPM.R")
source("./Code/IPM/BuildLTRE.R")

cov<-c("ppt","ppt","t","t")
grow_cont_clim<-ltre(dparam.denv=dparam.denv.g.clim,dlambda.dparam.indices=dlambda.dparam.indices.g.clim,
                       cov=cov,ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,
                       lambda=lambda_clim,vital.perturb="g",perturb=0.01,perturb.env=1,
                       gmodel.ltre=gmodel.clim,smodel.ltre=smodel.clim,
                       rmodel.ltre=rmodel.clim,gSD.ltre=growSD.clim,gba.clamp=F,gt.clamp=T,rba.clamp=F)

surv_cont_clim<-ltre(dparam.denv=dparam.denv.s.clim,dlambda.dparam.indices=dlambda.dparam.indices.s.clim,
                       cov=cov,ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,
                       lambda=lambda_clim,vital.perturb="s",perturb=0.01,perturb.env=1,
                       gmodel.ltre=gmodel.clim,smodel.ltre=smodel.clim,
                       rmodel.ltre=rmodel.clim,gSD.ltre=growSD.clim,gba.clamp=F,gt.clamp=T,rba.clamp=F)

recr_cont_clim<-ltre(dparam.denv=dparam.denv.r.clim,dlambda.dparam.indices=dlambda.dparam.indices.r.clim,
                       cov=cov,ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,
                       lambda=lambda_clim,vital.perturb="r",perturb=0.01,perturb.env=1,
                       gmodel.ltre=gmodel.clim,smodel.ltre=smodel.clim,
                       rmodel.ltre=rmodel.clim,gSD.ltre=growSD.clim,gba.clamp=F,gt.clamp=T,rba.clamp=F)

clim_total_ppt<-ltre_total(ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,lambda=lambda_clim,
                           gmodel.ltre=gmodel.clim,smodel.ltre=smodel.clim,rmodel.ltre=rmodel.clim,
                           gSD.ltre=growSD.clim,perturb.var="ppt",perturb.env=1)

clim_total_t<-ltre_total(ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,lambda=lambda_clim,
                         gmodel.ltre=gmodel.clim,smodel.ltre=smodel.clim,rmodel.ltre=rmodel.clim,
                         gSD.ltre=growSD.clim,perturb.var="t",perturb.env=1)

test_ppt<-sum(grow_cont_clim$cont_ppt,surv_cont_clim$cont_ppt,recr_cont_clim$cont_ppt)  
test_t<-sum(grow_cont_clim$cont_t,surv_cont_clim$cont_t,recr_cont_clim$cont_t)  

sens.ras<-ppt_yr_raster

writeRaster(grow_cont_clim$cont_ppt, "./Code/Elasticities/growth_cont_ppt_clim.tif", overwrite = T)
writeRaster(grow_cont_clim$cont_t, "./Code/Elasticities/growth_cont_t_clim.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_clim$sens[1,,]), "./Code/Elasticities/sens_grow_clim_ppt.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_clim$sens[2,,]), "./Code/Elasticities/sens_grow_clim_ppt2.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_clim$sens[3,,]), "./Code/Elasticities/sens_grow_clim_t.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_clim$sens[4,,]), "./Code/Elasticities/sens_grow_clim_t2.tif", overwrite = T)

writeRaster(surv_cont_clim$cont_ppt, "./Code/Elasticities/surv_cont_ppt_clim.tif", overwrite = T)
writeRaster(surv_cont_clim$cont_t, "./Code/Elasticities/surv_cont_t_clim.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_clim$sens[1,,]), "./Code/Elasticities/sens_surv_clim_ppt.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_clim$sens[2,,]), "./Code/Elasticities/sens_surv_clim_ppt2.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_clim$sens[3,,]), "./Code/Elasticities/sens_surv_clim_t.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_clim$sens[4,,]), "./Code/Elasticities/sens_surv_clim_t2.tif", overwrite = T)

writeRaster(recr_cont_clim$cont_ppt, "./Code/Elasticities/recr_cont_ppt_clim.tif", overwrite = T)
writeRaster(recr_cont_clim$cont_t, "./Code/Elasticities/recr_cont_t_clim.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_clim$sens[1,,]), "./Code/Elasticities/sens_recr_clim_ppt.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_clim$sens[2,,]), "./Code/Elasticities/sens_recr_clim_ppt2.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_clim$sens[3,,]), "./Code/Elasticities/sens_recr_clim_t.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_clim$sens[4,,]), "./Code/Elasticities/sens_recr_clim_t2.tif", overwrite = T)

pdf("./Output/growth_cont_clim.pdf")
plot(grow_cont_clim$cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(grow_cont_clim$cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/surv_cont_clim.pdf")
plot(surv_cont_clim$cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(surv_cont_clim$cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/recr_cont_clim.pdf")
plot(recr_cont_clim$cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(recr_cont_clim$cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

ppt_ltre<-sum(grow_cont_clim$cont_ppt,surv_cont_clim$cont_ppt,recr_cont_clim$cont_ppt)
plot(ppt_ltre, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
t_ltre<-sum(grow_cont_clim$cont_t,surv_cont_clim$cont_t,recr_cont_clim$cont_t)
plot(t_ltre, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)

writeRaster(clim_total_ppt, "./Code/Elasticities/ppt_ltre_clim.tif", overwrite = T)
writeRaster(clim_total_t, "./Code/Elasticities/t_ltre_clim.tif", overwrite = T)

# Climate + competition, no interactions
cov<-c("ba","ba","ppt","ppt","t","t")
grow_cont_climcomp<-ltre(dparam.denv=dparam.denv.g.climcomp,
                         dlambda.dparam.indices=dlambda.dparam.indices.g.climcomp,
                         cov=cov,ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,
                         lambda=lambda_climcomp,vital.perturb="g",perturb=0.01,perturb.env=1,
                         gmodel.ltre=gmodel.clim.comp,smodel.ltre=smodel.clim.comp,
                         rmodel.ltre=rmodel.clim.comp,gSD.ltre=growSD.clim.comp)

surv_cont_climcomp<-ltre(dparam.denv=dparam.denv.s.climcomp,
                     dlambda.dparam.indices=dlambda.dparam.indices.s.climcomp,
                     cov=cov,ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,
                     lambda=lambda_climcomp,vital.perturb="s",perturb=0.01,perturb.env=1, 
                     gmodel.ltre=gmodel.clim.comp,smodel.ltre=smodel.clim.comp,
                     rmodel.ltre=rmodel.clim.comp,gSD.ltre=growSD.clim.comp)

recr_cont_climcomp<-ltre(dparam.denv=dparam.denv.r.climcomp,
                         dlambda.dparam.indices=dlambda.dparam.indices.r.climcomp,
                         cov=cov,ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,
                         lambda=lambda_climcomp,vital.perturb="r",perturb=0.01,perturb.env=1,
                         gmodel.ltre=gmodel.clim.comp,smodel.ltre=smodel.clim.comp,
                         rmodel.ltre=rmodel.clim.comp,gSD.ltre=growSD.clim.comp)

climcomp_total_ppt<-ltre_total(ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,lambda=lambda_climcomp,
                           gmodel.ltre=gmodel.clim.comp,smodel.ltre=smodel.clim.comp,
                           rmodel.ltre=rmodel.clim.comp,
                           gSD.ltre=growSD.clim.comp,perturb.var="ppt",perturb.env=1)

climcomp_total_t<-ltre_total(ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,lambda=lambda_climcomp,
                         gmodel.ltre=gmodel.clim.comp,smodel.ltre=smodel.clim.comp,
                         rmodel.ltre=rmodel.clim.comp,
                         gSD.ltre=growSD.clim.comp,perturb.var="t",perturb.env=1)

climcomp_total_ba<-ltre_total(ppt=ppt_yr_raster,tmp=t_yr_raster,ba=ba_raster,lambda=lambda_climcomp,
                             gmodel.ltre=gmodel.clim.comp,smodel.ltre=smodel.clim.comp,
                             rmodel.ltre=rmodel.clim.comp,
                             gSD.ltre=growSD.clim.comp,perturb.var="ba",perturb.env=1)

test_ppt<-sum(grow_cont_climcomp$cont_ppt,surv_cont_climcomp$cont_ppt,recr_cont_climcomp$cont_ppt)  
test_t<-sum(grow_cont_climcomp$cont_t,surv_cont_climcomp$cont_t,recr_cont_climcomp$cont_t)  
test_ba<-sum(grow_cont_climcomp$cont_ba,surv_cont_climcomp$cont_ba,recr_cont_climcomp$cont_ba)  

writeRaster(grow_cont_climcomp$cont_ba, "./Code/Elasticities/growth_cont_ba_climcomp.tif", overwrite = T)
writeRaster(grow_cont_climcomp$cont_ppt, "./Code/Elasticities/growth_cont_ppt_climcomp.tif", overwrite = T)
writeRaster(grow_cont_climcomp$cont_t, "./Code/Elasticities/growth_cont_t_climcomp.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_climcomp$sens[1,,]), "./Code/Elasticities/sens_grow_climcomp_ba.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_climcomp$sens[2,,]), "./Code/Elasticities/sens_grow_climcomp_ba2.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_climcomp$sens[3,,]), "./Code/Elasticities/sens_grow_climcomp_ppt.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_climcomp$sens[4,,]), "./Code/Elasticities/sens_grow_climcomp_ppt2.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_climcomp$sens[5,,]), "./Code/Elasticities/sens_grow_climcomp_t.tif", overwrite = T)
writeRaster(setValues(sens.ras, grow_cont_climcomp$sens[6,,]), "./Code/Elasticities/sens_grow_climcomp_t2.tif", overwrite = T)

writeRaster(surv_cont_climcomp$cont_ba, "./Code/Elasticities/surv_cont_ba_climcomp.tif", overwrite = T)
writeRaster(surv_cont_climcomp$cont_ppt, "./Code/Elasticities/surv_cont_ppt_climcomp.tif", overwrite = T)
writeRaster(surv_cont_climcomp$cont_t, "./Code/Elasticities/surv_cont_t_climcomp.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_climcomp$sens[1,,]), "./Code/Elasticities/sens_surv_climcomp_ba.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_climcomp$sens[2,,]), "./Code/Elasticities/sens_surv_climcomp_ba2.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_climcomp$sens[3,,]), "./Code/Elasticities/sens_surv_climcomp_ppt.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_climcomp$sens[4,,]), "./Code/Elasticities/sens_surv_climcomp_ppt2.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_climcomp$sens[5,,]), "./Code/Elasticities/sens_surv_climcomp_t.tif", overwrite = T)
writeRaster(setValues(sens.ras, surv_cont_climcomp$sens[6,,]), "./Code/Elasticities/sens_surv_climcomp_t2.tif", overwrite = T)

writeRaster(recr_cont_climcomp$cont_ba, "./Code/Elasticities/recr_cont_ba_climcomp.tif", overwrite = T)
writeRaster(recr_cont_climcomp$cont_ppt, "./Code/Elasticities/recr_cont_ppt_climcomp.tif", overwrite = T)
writeRaster(recr_cont_climcomp$cont_t, "./Code/Elasticities/recr_cont_t_climcomp.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_climcomp$sens[1,,]), "./Code/Elasticities/sens_recr_climcomp_ba.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_climcomp$sens[2,,]), "./Code/Elasticities/sens_recr_climcomp_ba2.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_climcomp$sens[3,,]), "./Code/Elasticities/sens_recr_climcomp_ppt.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_climcomp$sens[4,,]), "./Code/Elasticities/sens_recr_climcomp_ppt2.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_climcomp$sens[5,,]), "./Code/Elasticities/sens_recr_climcomp_t.tif", overwrite = T)
writeRaster(setValues(sens.ras, recr_cont_climcomp$sens[6,,]), "./Code/Elasticities/sens_recr_climcomp_t2.tif", overwrite = T)

pdf("./Output/growth_cont_climcomp.pdf")
plot(grow_cont_climcomp$cont_ba, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(grow_cont_climcomp$cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(grow_cont_climcomp$cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/surv_cont_climcomp.pdf")
plot(surv_cont_climcomp$cont_ba, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(surv_cont_climcomp$cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(surv_cont_climcomp$cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

pdf("./Output/recr_cont_climcomp.pdf")
plot(recr_cont_climcomp$cont_ba, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(recr_cont_climcomp$cont_ppt, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(recr_cont_climcomp$cont_t, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

ba_ltre<-sum(growth_cont_ba,surv_cont_ba,recr_cont_ba)
plot(ba_ltre, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
ppt_ltre<-sum(growth_cont_ppt,surv_cont_ppt,recr_cont_ppt)
plot(ppt_ltre, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
t_ltre<-sum(growth_cont_t,surv_cont_t,recr_cont_t)
plot(t_ltre, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)

writeRaster(ppt_ltre, "./Code/Elasticities/ppt_ltre_climcomp.tif", overwrite = T)
writeRaster(t_ltre, "./Code/Elasticities/t_ltre_climcomp.tif", overwrite = T)

ba_ltre_total<-dlambda
ppt_ltre_total<-dlambda
t_ltre_total<-dlambda

pdf("./Output/ltre_totals_climcomp.pdf")
plot(ba_ltre_total, main = "Basal Area"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(ppt_ltre_total, main = "MAP"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(t_ltre_total, main = "MAT"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

writeRaster(climcomp_total_ba, "./Code/Elasticities/ba_ltre_total_clim.tif", overwrite = T)
writeRaster(climcomp_total_ppt, "./Code/Elasticities/ppt_ltre_total_clim.tif", overwrite = T)
writeRaster(climcomp_total_t, "./Code/Elasticities/t_ltre_total_clim.tif", overwrite = T)


### Elevation - lambda "LTRE"
source("./Code/IPM/BuildIPM.R")
mytheme2<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.text=element_text(size=12),legend.title=element_text(size=12),
                legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
                axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
                axis.line.x = element_line(color="black", size = 0.3),
                axis.line.y = element_line(color="black", size = 0.3))

FIA2 <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)

env<-cbind(FIA2$BALIVE,FIA2$PPT_yr_norm,FIA2$T_yr_norm)
elev_env<-manova(env ~ FIA2$elev)

FIA2$PLT_CN_factor<-as.factor(FIA2$plot)
k=5

elev_ba <-  bam(BALIVE ~ s(elev,k=k), data = FIA2)
elev_ppt <-  bam(PPT_yr_norm ~ s(elev,k=k), data = FIA2)
elev_t <-  bam(T_yr_norm ~ s(elev,k=k), data = FIA2)

elev_seq <- seq(min(FIA2$elev,na.rm=T),max(FIA2$elev,na.rm=T),by=100)

elev_env<-data.frame(Elevation = elev_seq,
                     BA = predict(elev_ba,newdata = data.frame(elev=elev_seq)),
                     MAP = predict(elev_ppt,newdata = data.frame(elev=elev_seq)),
                     MAT = predict(elev_t,newdata = data.frame(elev=elev_seq)))

save(FIA2,elev_ba,elev_ppt,elev_t,file="./Output/elev_models.rda")

min_elev_pied <- min(subset(FIA2,PApied==1)$elev)
max_elev_pied <- max(subset(FIA2,PApied==1)$elev)
save(min_elev_pied,max_elev_pied, file="./Output/elev_limits.rda")

min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500

lambda_elev<-numeric(0)
ssd_elev<-matrix(NA,n_dim,length(elev_seq))
  
for (i in 1:length(elev_seq)) {
    # Extract climate for cell
    pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"

    # Calculate lambda
    K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    lambda_val <- Re(eigen(K)$values[1])
    ssd<-Re(eigen(K)$vectors[,1])/sum(Re(eigen(K)$vectors[,1]))
    print(lambda_val)
    lambda_elev[i] <- lambda_val
    ssd_elev[,i] <- ssd
}

lambda_elev_ltre<-matrix(NA,length(elev_seq),3)
perturb<-1
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  for(j in 1:3){
    pred_data_new<-pred_data
    pred_data_new[1,j] <- pred_data_new[1,j]+perturb
  # Calculate lambda
  K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
             rmodel=rmodel.int.gam, gSD=growSD.int.gam,
             data=pred_data_new,
             s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  lambda_val <- Re(eigen(K)$values[1])
  lambda_elev_ltre[i,j] <- lambda_val
  }
  print(lambda_val)
}

d_lambda<-lambda_elev_ltre-lambda_elev

denv_delev<-matrix(NA,length(elev_seq),3)
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)" 
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  pred_data_new<-data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=(elev_seq[i]+1)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            T_yr=predict(elev_t, newdata = data.frame(elev=(elev_seq[i]+1)), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            BALIVE=predict(elev_ba, newdata = data.frame(elev=(elev_seq[i]+1)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
    for(j in 1:3){
    d <- pred_data_new[1,j]-pred_data[1,j]
    # Calculate lambda
    denv_delev[i,j] <- d
    }
  print(d)
}

elev_ltre<-d_lambda*denv_delev

test<-rowSums(elev_ltre)
test2<-(lambda_elev[2:length(elev_seq)]-lambda_elev[1:(length(elev_seq)-1)])/100
elev_seq2<-(elev_seq[1:(length(elev_seq)-1)]+elev_seq[2:length(elev_seq)])/2

plot(elev_seq2,test2)
lines(elev_seq,test)

lam_elev_data_c<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_c$Model="c"
lam_elev_data_cc<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_cc$Model="cc"
lam_elev_data_i<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_i$Model="i"

lam_elev_data=rbind(lam_elev_data_c,lam_elev_data_cc,lam_elev_data_i)

write.csv(lam_elev_data,"./Output/lam_elev_data.csv",row.names=F)

ltre_data<-data.frame(Elev=c(elev_seq,elev_seq,elev_seq,elev_seq),
                      Env_c=c(elev_ltre[,1],elev_ltre[,2],elev_ltre[,3],test),
                      dlam=c(d_lambda[,1],d_lambda[,2],d_lambda[,3],rowSums(d_lambda)),
                      denv=c(denv_delev[,1],denv_delev[,2],denv_delev[,3],rowSums(denv_delev)),
                      Predictor=c(rep("MAP",length(elev_seq)),rep("MAT",length(elev_seq)),rep("BA",length(elev_seq)),rep("Total",length(elev_seq))))
ltre_data$Env_cc<-c(elev_ltre[,1],elev_ltre[,2],elev_ltre[,3],test)
ltre_data$Env_i<-c(elev_ltre[,1],elev_ltre[,2],elev_ltre[,3],test)

write.csv(ltre_data,"./Output/ltre_data_env.csv",row.names=F)

# Contribution of vital rates
lambda_elev_ltre_vital<-matrix(NA,length(elev_seq),3)
perturb=0.01
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
    # Calculate lambda
    K_g<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, gperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    K_s<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                 rmodel=rmodel.int.gam, gSD=growSD.int.gam,
                 data=pred_data, sperturb=perturb,
                 s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    K_r<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                 rmodel=rmodel.int.gam, gSD=growSD.int.gam,
                 data=pred_data, rperturb=perturb,
                 s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    lambda_elev_ltre_vital[i,1] <- Re(eigen(K_g)$values[1])
    lambda_elev_ltre_vital[i,2] <- Re(eigen(K_s)$values[1])
    lambda_elev_ltre_vital[i,3] <- Re(eigen(K_r)$values[1])
    
  print(i)
}

d_lambda_vital<-(lambda_elev-lambda_elev_ltre_vital)/perturb

dvital_delev<-matrix(NA,length(elev_seq),3)
perturb=10
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)" 
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  # Perturb elevation and calculate new climate
  pred_data_new<-data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=(elev_seq[i]+perturb)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            T_yr=predict(elev_t, newdata = data.frame(elev=(elev_seq[i]+perturb)), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            BALIVE=predict(elev_ba, newdata = data.frame(elev=(elev_seq[i]+perturb)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  # Calculate change in vital rate
  g<-g.mean(model=gmodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data)
  g_perturb<-g.mean(model=gmodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data_new)
  s<-s.x(model=smodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data)
  s_perturb<-s.x(model=smodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data_new)
  r<-f.mean(model=rmodel.int.gam, data=pred_data)
  r_perturb<-f.mean(model=rmodel.int.gam, data=pred_data_new)

  dvital_delev[i,1] <- g-g_perturb
  dvital_delev[i,2] <- s-s_perturb
  dvital_delev[i,3] <- r-r_perturb
  
  print(i)
}

elev_ltre_vital<-d_lambda_vital*(dvital_delev)/10

test<-rowSums(elev_ltre_vital)
test2<-(lambda_elev[2:length(elev_seq)]-lambda_elev[1:(length(elev_seq)-1)])/100
elev_seq2<-(elev_seq[1:(length(elev_seq)-1)]+elev_seq[2:length(elev_seq)])/2

plot(elev_seq2,test2)
lines(elev_seq,test)

ltre_data_vital<-data.frame(Elev=c(elev_seq,elev_seq,elev_seq,elev_seq),
                      Env_c=c(elev_ltre_vital[,1],elev_ltre_vital[,2],elev_ltre_vital[,3],test),
                      dlam=c(d_lambda_vital[,1],d_lambda_vital[,2],d_lambda_vital[,3],rowSums(d_lambda_vital)),
                      denv=c((dvital_delev[,1])/10,(dvital_delev[,2])/10,(dvital_delev[,3])/10,rowSums(dvital_delev/10)),
                      Rate=c(rep("Growth",length(elev_seq)),rep("Survival",length(elev_seq)),rep("Recruitment",length(elev_seq)),rep("Total",length(elev_seq))))
ltre_data_vital$Env_cc<-c(elev_ltre_vital[,1],elev_ltre_vital[,2],elev_ltre_vital[,3],test)
ltre_data_vital$Env_i<-c(elev_ltre_vital[,1],elev_ltre_vital[,2],elev_ltre_vital[,3],test)

write.csv(ltre_data_vital,"./Output/ltre_data_vital.csv",row.names=F)

# Elasticity
elev_elast<-matrix(NA,length(elev_seq),3)
perturb=0.01
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  # Calculate lambda
  K_g<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, gperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_s<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, sperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_r<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, rperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  elev_elast[i,1] <- Re(eigen(K_g)$values[1])
  elev_elast[i,2] <- Re(eigen(K_s)$values[1])
  elev_elast[i,3] <- Re(eigen(K_r)$values[1])
  
  print(i)
}

#lam_elev_data<-read.csv("./Output/lam_elev_data.csv")

#lambda_elev<-lam_elev_data_i$Lambda #[which(lam_elev_data$Model=="cc")]
elast_vital<-(elev_elast-lambda_elev)/(lambda_elev*perturb)

test<-rowSums(elast_vital)-elast_vital[,1]

elast_data_vital<-data.frame(Elev=c(elev_seq,elev_seq,elev_seq),
                            Elast_c=c(elast_vital[,1],elast_vital[,2],elast_vital[,3]),
                            Rate=c(rep("Growth",length(elev_seq)),rep("Survival",length(elev_seq)),rep("Recruitment",length(elev_seq))))

elast_data_vital$Elast_cc<-c(elast_vital[,1],elast_vital[,2],elast_vital[,3])
elast_data_vital$Elast_i<-c(elast_vital[,1],elast_vital[,2],elast_vital[,3])

write.csv(elast_data_vital,"./Output/elast_vital.csv",row.names=F)


# Perturb vitals
perturb_seq=seq(0,-1,-0.01)
elev_perturb<-matrix(NA,length(perturb_seq),3)
pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=max_elev_pied), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                        T_yr=predict(elev_t, newdata = data.frame(elev=max_elev_pied), #,PLT_CN_factor=14546600020004
                                     type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                        BALIVE=predict(elev_ba, newdata = data.frame(elev=max_elev_pied), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"

for (i in 1:length(perturb_seq)) {

  # Calculate lambda
  K_g<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, gperturb=perturb_seq[i],
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_s<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, sperturb=perturb_seq[i],
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_r<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, rperturb=perturb_seq[i],
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  elev_perturb[i,1] <- Re(eigen(K_g)$values[1])
  elev_perturb[i,2] <- Re(eigen(K_s)$values[1])
  elev_perturb[i,3] <- Re(eigen(K_r)$values[1])
  
  print(i)
}


perturb_data_vital<-data.frame(Perturb=c(perturb_seq,perturb_seq,perturb_seq),
                             Elast_c=c(elev_perturb[,1],elev_perturb[,2],elev_perturb[,3]),
                             Rate=c(rep("Growth",length(perturb_seq)),rep("Survival",length(perturb_seq)),rep("Recruitment",length(perturb_seq))))

perturb_data_vital$Elast_cc<-c(elev_perturb[,1],elev_perturb[,2],elev_perturb[,3])
perturb_data_vital$Elast_i<-c(elev_perturb[,1],elev_perturb[,2],elev_perturb[,3])

perturb_data_vital2<-gather(perturb_data_vital,'Elast_c','Elast_cc','Elast_i',key="Model",value="Lambda")

write.csv(perturb_data_vital2,"./Output/perturb_vital.csv",row.names=F)

perturb_plot_growth<-ggplot(data=subset(perturb_data_vital2,Rate=="Growth"),aes(x=(-1*Perturb),y=Lambda,col=Model))+
  geom_abline(intercept=1,slope=0)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Elast_c","Elast_cc","Elast_i"), 
                      values=c("Elast_c"="#1b9e77","Elast_cc"="#d95f02","Elast_i"="#7570b3"),
                      labels=c("Clim","ClimComp","ClimCompInt"),
                      name="Model")+
  labs(x = "Percent perturbation", y ="Lambda")+
  mytheme2
ggsave(file="perturb_plot_growth.png", plot=perturb_plot_growth,width=4,height=3,units="in",dpi=600)

x_c_s<-(-1)*max(subset(perturb_data_vital2,Rate=="Survival" & Model=="Elast_c" & Lambda<1)$Perturb)
x_cc_s<-(-1)*max(subset(perturb_data_vital2,Rate=="Survival" & Model=="Elast_cc" & Lambda<1)$Perturb)
x_i_s<-(-1)*max(subset(perturb_data_vital2,Rate=="Survival" & Model=="Elast_i" & Lambda<1)$Perturb)

perturb_plot_survival<-ggplot(data=subset(perturb_data_vital2,Rate=="Survival"),aes(x=(-1*Perturb),y=Lambda,col=Model))+
  coord_cartesian(xlim = c(0,0.1))+
  geom_abline(intercept=1,slope=0)+
  geom_segment(x=x_c_s,xend=x_c_s,y=1,yend=0,col="#1b9e77",linetype="dashed")+
  geom_segment(x=x_cc_s,xend=x_cc_s,y=1,yend=0,col="#d95f02",linetype="dashed")+
  geom_segment(x=x_i_s,xend=x_i_s,y=1,yend=0,col="#7570b3",linetype="dashed")+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Elast_c","Elast_cc","Elast_i"), 
                      values=c("Elast_c"="#1b9e77","Elast_cc"="#d95f02","Elast_i"="#7570b3"),
                      labels=c("Clim","ClimComp","ClimCompInt"),
                      name="Model")+
  labs(x = "Percent perturbation", y ="Lambda")+
  mytheme2
ggsave(file="perturb_plot_survival.png", plot=perturb_plot_survival,width=4,height=3,units="in",dpi=600)

x_c_r<-(-1)*max(subset(perturb_data_vital2,Rate=="Recruitment" & Model=="Elast_c" & Lambda<1)$Perturb)
x_cc_r<-(-1)*max(subset(perturb_data_vital2,Rate=="Recruitment" & Model=="Elast_cc" & Lambda<1)$Perturb)
x_i_r<-(-1)*max(subset(perturb_data_vital2,Rate=="Recruitment" & Model=="Elast_i" & Lambda<1)$Perturb)

perturb_plot_recruit<-ggplot(data=subset(perturb_data_vital2,Rate=="Recruitment"),aes(x=(-1*Perturb),y=Lambda,col=Model))+
  geom_abline(intercept=1,slope=0)+
  geom_segment(x=x_c_r,xend=x_c_r,y=1,yend=0,col="#1b9e77",linetype="dashed")+
  geom_segment(x=x_cc_r,xend=x_cc_r,y=1,yend=0,col="#d95f02",linetype="dashed")+
  geom_segment(x=x_i_r,xend=x_i_r,y=1,yend=0,col="#7570b3",linetype="dashed")+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Elast_c","Elast_cc","Elast_i"), 
                      values=c("Elast_c"="#1b9e77","Elast_cc"="#d95f02","Elast_i"="#7570b3"),
                      labels=c("Clim","ClimComp","ClimCompInt"),
                      name="Model")+
  labs(x = "Percent perturbation", y ="Lambda")+
  mytheme2
ggsave(file="perturb_plot_recruit.png", plot=perturb_plot_recruit,width=4,height=3,units="in",dpi=600)
