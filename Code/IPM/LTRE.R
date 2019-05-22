library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(glmmTMB)

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
