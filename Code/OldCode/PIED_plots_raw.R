library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(glmmTMB)
library(coefplot)
library(ggplot2)
library(ggeffects)
library(cowplot)
library(effects)
library(DHARMa)
library(lme4)
library(grid)
library(MuMIn)

load("./Code/IPM/GrRescaling.Rdata")
load("./Code/IPM/SurvRescaling.Rdata")
load("./Code/IPM/RecruitRescaling.Rdata")
load("./Code/IPM/recrstats.rda")

# Data prep
grdata <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

# Only keep trees that didn't die
grdata <- subset(grdata, STATUSCD == 1) #18204

# Create increment columns
# note that growth increments need to be moved to the positive realm (by adding a constant)
# IF log transform is used
grdata$AGB_INCR <- grdata$DRYBIO_AG_DIFF / grdata$CENSUS_INTERVAL
grdata$DIA_INCR <- grdata$DIA_DIFF / grdata$CENSUS_INTERVAL
grdata$BA_INCR <- grdata$BA_DIFF / grdata$CENSUS_INTERVAL

grdata.scaled <- grdata %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                          -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                          -CENSUS_INTERVAL,
                                                          -AGB_INCR, -DIA_INCR, -BA_INCR))

survData <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

# Create increment columns
# not needed for survival/mort analysis
survData$AGB_INCR <- survData$DRYBIO_AG_DIFF / survData$CENSUS_INTERVAL
survData$DIA_INCR <- survData$DIA_DIFF / survData$CENSUS_INTERVAL
survData$BA_INCR <- survData$BA_DIFF / survData$CENSUS_INTERVAL
survData$log.size <- log(survData$PREVDIA)
survData$log.BALIVE <- log(survData$BALIVE)

# Recode status
survData$surv <- ifelse(survData$STATUSCD == 2, 0, 1)
survData$mort <- ifelse(survData$STATUSCD == 1, 0, 1)

# remove cases where BALIVE at time 1 = zero (should be impossible)
# survData <- subset(survData, log.BALIVE > 0) 
survData.2 <- subset(survData, BALIVE > 0) # goes from 20329 to 20161

# remove conditions where fire or harvest occurred
survData.3 <- survData[!(survData$DSTRBCD1 %in% c(30, 31, 32, 80)), ] # goes from 20329 to 19867

# standardize covariates
#ELS update: no AGENTCD, DSTRBCD1, DSTRBCD2, DSTRBCD3 in dataframe, so I removed them from the following code
survData.scaled <- survData %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                              -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                              -CENSUS_INTERVAL,
                                                              AGENTCD, DSTRBCD1, DSTRBCD2, DSTRBCD3,
                                                              -AGB_INCR, -DIA_INCR, -BA_INCR,
                                                              -surv, -mort))

survData2.scaled <- survData.2 %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                                 -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                                 -CENSUS_INTERVAL,
                                                                 AGENTCD, DSTRBCD1, DSTRBCD2, DSTRBCD3,
                                                                 -AGB_INCR, -DIA_INCR, -BA_INCR,
                                                                 -surv, -mort))

survData3.scaled <- survData.3 %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                                 -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                                 -CENSUS_INTERVAL,
                                                                 AGENTCD, DSTRBCD1, DSTRBCD2, DSTRBCD3,
                                                                 -AGB_INCR, -DIA_INCR, -BA_INCR,
                                                                 -surv, -mort))

rdata <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)
rdata <- subset(rdata, PIEDadults1 > 0)
rdata$BA.all <- rdata$BA.PIED + rdata$BA.notPIED
rdata.scaled <- rdata %>% mutate_at(scale, .vars = vars(-plot, -lat, -lon, -elev, -PApied,
                                                        -state, -county, -plotID, -CONDID, 
                                                        -measYear, -plotPrev, -PREV_MEASYEAR,
                                                        -CENSUS_INTERVAL, -recruits1, -recruits12,
                                                        -AGB_intra, -BA.PIED, -PIEDadults1,
                                                        -PIEDadults4, -PIEDadults8, -cumDIA.PIED))
ppt_yr_raster <- raster("./ClimateData/PRISM/Normals/PPT_year.tif")
t_yr_raster <- raster("./ClimateData/PRISM/Normals/T_year.tif")
ba_raster <- raster("./BA/balive_RF.tif")

#Growth
ncuts=30
chopsize_dia<-cut(grdata$PREVDIA,ncuts)
chopsize_ba<-cut(grdata$BALIVE,ncuts)
chopsize_PPT<-cut(grdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(grdata$T_yr_norm,ncuts)
grow_binned_dia<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_dia),mean,na.rm=T))
count_binned_dia<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_dia),length))
dia_binned<-as.vector(sapply(split(grdata$PREVDIA,chopsize_dia),mean,na.rm=T))
grow_binned_ba<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_ba),length))
ba_binned<-as.vector(sapply(split(grdata$BALIVE,chopsize_ba),mean,na.rm=T))
grow_binned_PPT<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_PPT),length))
PPT_binned<-as.vector(sapply(split(grdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
grow_binned_T<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(grdata$DIA_INCR,chopsize_T),length))
T_binned<-as.vector(sapply(split(grdata$T_yr_norm,chopsize_T),mean,na.rm=T))
g_binned<-as.data.frame(cbind(grow_binned_dia,grow_binned_ba,grow_binned_PPT,grow_binned_T,
                              dia_binned,ba_binned,PPT_binned,T_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))
names(g_binned)<-c("grow_dia","grow_ba","grow_PPT","grow_T","PREVDIA","BALIVE","PPT","T",
                   "count_dia","count_ba","count_PPT","count_T")

g_fun<-function(dia,ba,ppt,t,model,clampba=F,clampt=F){
  gdata=data.frame(PREVDIA=dia,BALIVE = ifelse(clampba==T & ba >190, 190, ba),
                   PPT_yr_norm = ppt, T_yr_norm = ifelse(clampt==T & t >10.6, 10.6, t))
  scaled.gdata = data.frame(scale(gdata,
                                  scale = gr.scaling$scale[match(names(gdata), 
                                                                 names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), 
                                                                   names(gr.scaling$center))]))  
  return(predict(model, newdata = scaled.gdata, re.form = NA))
}
means<-c(mean(grdata$PREVDIA,na.rm=T),mean(grdata$BALIVE,na.rm=T),
         mean(grdata$PPT_yr_norm,na.rm=T),mean(grdata$T_yr_norm,na.rm=T))
seq<-data.frame(dia=seq(0.5*min(grdata$PREVDIA),1.2*max(grdata$PREVDIA),length=50),
                ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(grdata$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))
grplot_data_clim<-cbind(data.frame(dia_pred=g_fun(seq$dia,means[2],means[3],means[4],gmodel.clim),
                                   ppt_pred=g_fun(means[1],means[2],seq$ppt,means[4],gmodel.clim),
                                   t_pred=g_fun(means[1],means[2],means[3],seq$t,gmodel.clim),
                                   dia_pred_c=g_fun(seq$dia,means[2],means[3],means[4],gmodel.clim,
                                                    clampba=T,clampt=T),
                                   ppt_pred_c=g_fun(means[1],means[2],seq$ppt,means[4],gmodel.clim,
                                                    clampba=T,clampt=T),
                                   t_pred_c=g_fun(means[1],means[2],means[3],seq$t,gmodel.clim,
                                                  clampba=T,clampt=T)),seq)
grplot_data_climcomp<-cbind(data.frame(dia_pred=g_fun(seq$dia,means[2],means[3],means[4],
                                                      gmodel.clim.comp),
                                       ba_pred=g_fun(means[1],seq$ba,means[3],means[4],gmodel.clim.comp),
                                       ppt_pred=g_fun(means[1],means[2],seq$ppt,means[4],gmodel.clim.comp),
                                       t_pred=g_fun(means[1],means[2],means[3],seq$t,gmodel.clim.comp),
                                       dia_pred_c=g_fun(seq$dia,means[2],means[3],means[4],gmodel.clim.comp,
                                                        clampba=T,clampt=F),
                                       ba_pred_c=g_fun(means[1],seq$ba,means[3],means[4],gmodel.clim.comp,
                                                       clampba=T,clampt=F),
                                       ppt_pred_c=g_fun(means[1],means[2],seq$ppt,means[4],gmodel.clim.comp,
                                                        clampba=T,clampt=F),
                                       t_pred_c=g_fun(means[1],means[2],means[3],seq$t,gmodel.clim.comp,
                                                      clampba=T,clampt=F)),seq)
grplot_data_int<-cbind(data.frame(dia_pred=g_fun(seq$dia,means[2],means[3],means[4],
                                                 gmodel.int),
                                  ba_pred=g_fun(means[1],seq$ba,means[3],means[4],gmodel.int),
                                  ppt_pred=g_fun(means[1],means[2],seq$ppt,means[4],gmodel.int),
                                  t_pred=g_fun(means[1],means[2],means[3],seq$t,gmodel.int),
                                  dia_pred_c=g_fun(seq$dia,means[2],means[3],means[4],gmodel.int,
                                                   clampba=T,clampt=F),
                                  ba_pred_c=g_fun(means[1],seq$ba,means[3],means[4],gmodel.int,
                                                  clampba=T,clampt=F),
                                  ppt_pred_c=g_fun(means[1],means[2],seq$ppt,means[4],gmodel.int,
                                                   clampba=T,clampt=F),
                                  t_pred_c=g_fun(means[1],means[2],means[3],seq$t,gmodel.int,
                                                 clampba=T,clampt=F)),seq)

#Survival
ncuts=50
chopsize_dia<-cut(survData$PREVDIA,ncuts)
chopsize_ba<-cut(survData$BALIVE,ncuts)
chopsize_PPT<-cut(survData$PPT_yr_norm,ncuts)
chopsize_PPTc<-cut(survData$PPT_c_norm,ncuts)
chopsize_PPTpf<-cut(survData$PPT_pf_norm,ncuts)
chopsize_PPTfs<-cut(survData$PPT_fs_norm,ncuts)
chopsize_PPTm<-cut(survData$PPT_m_norm,ncuts)
chopsize_T<-cut(survData$T_yr_norm,ncuts)
chopsize_Tc<-cut(survData$T_c_norm,ncuts)
chopsize_Tpf<-cut(survData$T_pf_norm,ncuts)
chopsize_Tfs<-cut(survData$T_fs_norm,ncuts)
chopsize_Tm<-cut(survData$T_m_norm,ncuts)
surv_binned_dia<-as.vector(sapply(split(survData$mort,chopsize_dia),mean,na.rm=T))
dia_binned<-as.vector(sapply(split(survData$PREVDIA,chopsize_dia),mean,na.rm=T))
count_binned_dia<-as.vector(sapply(split(survData$mort,chopsize_dia),length))
surv_binned_ba<-as.vector(sapply(split(survData$mort,chopsize_ba),mean,na.rm=T))
ba_binned<-as.vector(sapply(split(survData$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(survData$mort,chopsize_ba),length))
surv_binned_PPT<-as.vector(sapply(split(survData$mort,chopsize_PPT),mean,na.rm=T))
PPT_binned<-as.vector(sapply(split(survData$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(survData$mort,chopsize_PPT),length))
surv_binned_PPTc<-as.vector(sapply(split(survData$mort,chopsize_PPTc),mean,na.rm=T))
PPTc_binned<-as.vector(sapply(split(survData$PPT_c_norm,chopsize_PPTc),mean,na.rm=T))
surv_binned_PPTpf<-as.vector(sapply(split(survData$mort,chopsize_PPTpf),mean,na.rm=T))
PPTpf_binned<-as.vector(sapply(split(survData$PPT_pf_norm,chopsize_PPTpf),mean,na.rm=T))
surv_binned_PPTfs<-as.vector(sapply(split(survData$mort,chopsize_PPTfs),mean,na.rm=T))
PPTfs_binned<-as.vector(sapply(split(survData$PPT_fs_norm,chopsize_PPTfs),mean,na.rm=T))
surv_binned_PPTm<-as.vector(sapply(split(survData$mort,chopsize_PPTm),mean,na.rm=T))
PPTm_binned<-as.vector(sapply(split(survData$PPT_m_norm,chopsize_PPTm),mean,na.rm=T))
surv_binned_T<-as.vector(sapply(split(survData$mort,chopsize_T),mean,na.rm=T))
T_binned<-as.vector(sapply(split(survData$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(survData$mort,chopsize_T),length))
surv_binned_Tc<-as.vector(sapply(split(survData$mort,chopsize_Tc),mean,na.rm=T))
Tc_binned<-as.vector(sapply(split(survData$T_c_norm,chopsize_Tc),mean,na.rm=T))
surv_binned_Tpf<-as.vector(sapply(split(survData$mort,chopsize_Tpf),mean,na.rm=T))
Tpf_binned<-as.vector(sapply(split(survData$T_pf_norm,chopsize_Tpf),mean,na.rm=T))
surv_binned_Tfs<-as.vector(sapply(split(survData$mort,chopsize_Tfs),mean,na.rm=T))
Tfs_binned<-as.vector(sapply(split(survData$T_fs_norm,chopsize_Tfs),mean,na.rm=T))
surv_binned_Tm<-as.vector(sapply(split(survData$mort,chopsize_Tm),mean,na.rm=T))
Tm_binned<-as.vector(sapply(split(survData$T_m_norm,chopsize_Tm),mean,na.rm=T))

s_binned<-as.data.frame(cbind(surv_binned_dia,surv_binned_ba,surv_binned_PPT,surv_binned_PPTc,
                              surv_binned_PPTpf,surv_binned_PPTfs,surv_binned_PPTm,
                              surv_binned_T,surv_binned_Tc,surv_binned_Tpf,
                              surv_binned_Tfs,surv_binned_Tm,dia_binned,ba_binned,
                              PPT_binned,PPTc_binned,PPTpf_binned,PPTfs_binned,PPTm_binned,
                              T_binned,Tc_binned,Tpf_binned,Tfs_binned,Tm_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))
names(s_binned)<-c("mort_dia","mort_ba","mort_PPT","mort_PPTc","mort_PPTpf","mort_PPTfs",
                   "mort_PPTm","mort_T","mort_Tc","mort_Tpf","mort_Tfs","mort_Tm",
                   "PREVDIA","BALIVE","PPT","PPT_c","PPT_pf","PPT_fs","PPT_m",
                   "T","T_c","T_pf","T_fs","T_m","count_dia","count_ba","count_PPT","count_T")

s_fun<-function(dia,ba,ppt,t,ci,model,clampt=F){
  sdata=data.frame(PREVDIA=dia,BALIVE=ba,PPT_yr_norm = ppt, 
                   T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), 
                                                                   names(surv.scaling$scale))], 
                                  center = surv.scaling$center[match(names(sdata), 
                                                                     names(surv.scaling$center))])) 
  scaled.sdata = cbind(scaled.sdata,CENSUS_INTERVAL = ci)
  return(predict(model, newdata = scaled.sdata, type = "response", re.form = NA))
}

means<-c(mean(survData$PREVDIA,na.rm=T),mean(survData$BALIVE,na.rm=T),
         mean(survData$PPT_yr_norm,na.rm=T),mean(survData$T_yr_norm,na.rm=T),
         mean(survData$CENSUS_INTERVAL))
seq<-data.frame(dia=seq(0.5*min(survData$PREVDIA),1.2*max(survData$PREVDIA),length=50),
                ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(survData$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))

splot_data_clim<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                 smodel.clim),
                                  ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim),
                                  t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.clim),
                                  dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                   smodel.clim,clampt=T),
                                  ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                   smodel.clim,clampt=T),
                                  t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.clim,
                                                 clampt=T)),seq)
splot_data_climcomp<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.clim.comp),
                                      ba_pred=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                    smodel.clim.comp),
                                      ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                     smodel.clim.comp),
                                      t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                   smodel.clim.comp),
                                      dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                       smodel.clim.comp,clampt=T),
                                      ba_pred_c=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                      smodel.clim.comp,clampt=T),
                                      ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                       smodel.clim.comp,clampt=T),
                                      t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                     smodel.clim.comp,clampt=T)),seq)
splot_data_fire<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                 smodel.clim.comp.fire),
                                  ba_pred=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                smodel.clim.comp.fire),
                                  ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                 smodel.clim.comp.fire),
                                  t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                               smodel.clim.comp.fire),
                                  dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                   smodel.clim.comp.fire,clampt=T),
                                  ba_pred_c=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                  smodel.clim.comp.fire,clampt=T),
                                  ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                   smodel.clim.comp.fire,clampt=T),
                                  t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                 smodel.clim.comp.fire,clampt=T)),seq)
splot_data_int<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                smodel.int),
                                 ba_pred=s_fun(means[1],seq$ba,means[3],means[4],means[5],smodel.int),
                                 ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],smodel.int),
                                 t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.int),
                                 dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],smodel.int,
                                                  clampt=T),
                                 ba_pred_c=s_fun(means[1],seq$ba,means[3],means[4],means[5],smodel.int,
                                                 clampt=T),
                                 ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],smodel.int,
                                                  clampt=T),
                                 t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.int,
                                                clampt=T)),seq)

splot_data_clim.lin<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.clim.lin),
                                      ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim.lin),
                                      t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.clim.lin),
                                      dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                       smodel.clim.lin,clampt=T),
                                      ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                       smodel.clim.lin,clampt=T),
                                      t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.clim.lin,
                                                     clampt=T)),seq)
splot_data_climcomp.lin<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                         smodel.clim.comp.lin),
                                          ba_pred=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                        smodel.clim.comp.lin),
                                          ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                         smodel.clim.comp.lin),
                                          t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                       smodel.clim.comp.lin),
                                          dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                           smodel.clim.comp.lin,clampt=T),
                                          ba_pred_c=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                          smodel.clim.comp.lin,clampt=T),
                                          ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                           smodel.clim.comp.lin,clampt=T),
                                          t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                         smodel.clim.comp.lin,clampt=T)),seq)
splot_data_fire.lin<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.clim.comp.fire.lin),
                                      ba_pred=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                    smodel.clim.comp.fire.lin),
                                      ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                     smodel.clim.comp.fire.lin),
                                      t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                   smodel.clim.comp.fire.lin),
                                      dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                       smodel.clim.comp.fire.lin,clampt=T),
                                      ba_pred_c=s_fun(means[1],seq$ba,means[3],means[4],means[5],
                                                      smodel.clim.comp.fire.lin,clampt=T),
                                      ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],
                                                       smodel.clim.comp.fire.lin,clampt=T),
                                      t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],
                                                     smodel.clim.comp.fire.lin,clampt=T)),seq)
splot_data_int.lin<-cbind(data.frame(dia_pred=s_fun(seq$dia,means[2],means[3],means[4],means[5],
                                                    smodel.int.lin),
                                     ba_pred=s_fun(means[1],seq$ba,means[3],means[4],means[5],smodel.int.lin),
                                     ppt_pred=s_fun(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.lin),
                                     t_pred=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.int.lin),
                                     dia_pred_c=s_fun(seq$dia,means[2],means[3],means[4],means[5],smodel.int.lin,
                                                      clampt=T),
                                     ba_pred_c=s_fun(means[1],seq$ba,means[3],means[4],means[5],smodel.int.lin,
                                                     clampt=T),
                                     ppt_pred_c=s_fun(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.lin,
                                                      clampt=T),
                                     t_pred_c=s_fun(means[1],means[2],means[3],seq$t,means[5],smodel.int.lin,
                                                    clampt=T)),seq)

#Recruitment
ncuts=30
chopsize_ba<-cut(rdata$BALIVE,ncuts)
chopsize_PPT<-cut(rdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(rdata$T_yr_norm,ncuts)
recr_binned_ba<-as.vector(sapply(split(rdata$recruits1,chopsize_ba),mean,na.rm=T))
ba_binned<-as.vector(sapply(split(rdata$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(rdata$recruits1,chopsize_ba),length))
recr_binned_PPT<-as.vector(sapply(split(rdata$recruits1,chopsize_PPT),mean,na.rm=T))
PPT_binned<-as.vector(sapply(split(rdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(rdata$recruits1,chopsize_PPT),length))
recr_binned_T<-as.vector(sapply(split(rdata$recruits1,chopsize_T),mean,na.rm=T))
T_binned<-as.vector(sapply(split(rdata$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(rdata$recruits1,chopsize_T),length))
r_binned<-as.data.frame(cbind(recr_binned_ba,recr_binned_PPT,recr_binned_T,
                              ba_binned,PPT_binned,T_binned,
                              count_binned_ba,count_binned_PPT,count_binned_T))
names(r_binned)<-c("recr_ba","recr_PPT","recr_T","BALIVE","PPT","T","count_ba","count_PPT","count_T")

r_fun<-function(ba,ppt,t,pa,ci,model,clampba=F){
  rdata=data.frame(BALIVE = ifelse(clampba==T & ba >204, 204, 
                                   ifelse(clampba==T & ba <93, 93, ba)),
                   PPT_yr_norm = ppt, T_yr_norm = t)
  scaled.rdata = data.frame(scale(rdata, 
                                  scale = r.scaling$scale[match(names(rdata), 
                                                                names(r.scaling$scale))], 
                                  center = r.scaling$center[match(names(rdata), 
                                                                  names(r.scaling$center))])) 
  scaled.rdata$CENSUS_INTERVAL = ci
  scaled.rdata$PIEDadults1 = pa
  return(predict(model, newdata = scaled.rdata, type = "response"))
}

means<-c(mean(rdata$BALIVE,na.rm=T),
         mean(rdata$PPT_yr_norm,na.rm=T),mean(rdata$T_yr_norm,na.rm=T),
         mean(rdata$PIEDadults1),mean(rdata$CENSUS_INTERVAL))
seq<-data.frame(ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(rdata$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))

rplot_data_clim<-cbind(data.frame(ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                 rmodel.clim),
                                  t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim),
                                  ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                   rmodel.clim,clampba=T),
                                  t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim,
                                                 clampba=T)),seq)
rplot_data_climcomp<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                    rmodel.clim.comp),
                                      ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                     rmodel.clim.comp),
                                      t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],
                                                   rmodel.clim.comp),
                                      ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                      rmodel.clim.comp,clampba=T),
                                      ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                       rmodel.clim.comp,clampba=T),
                                      t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],
                                                     rmodel.clim.comp,clampba=T)),seq)
rplot_data_int<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                               rmodel.int),
                                 ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],rmodel.int),
                                 t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.int),
                                 ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                 rmodel.int,clampba=T),
                                 ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],rmodel.int,
                                                  clampba=T),
                                 t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.int,
                                                clampba=T)),seq)

rplot_data_clim.lin<-cbind(data.frame(ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                     rmodel.clim.lin),
                                      t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim.lin),
                                      ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                       rmodel.clim.lin,clampba=T),
                                      t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim.lin,
                                                     clampba=T)),seq)
rplot_data_climcomp.lin<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                        rmodel.clim.comp.lin),
                                          ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                         rmodel.clim.comp.lin),
                                          t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],
                                                       rmodel.clim.comp.lin),
                                          ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                          rmodel.clim.comp.lin,clampba=T),
                                          ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                           rmodel.clim.comp.lin,clampba=T),
                                          t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],
                                                         rmodel.clim.comp.lin,clampba=T)),seq)
rplot_data_int.lin<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                   rmodel.int.lin),
                                     ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],rmodel.int.lin),
                                     t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.int.lin),
                                     ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                     rmodel.int.lin,clampba=T),
                                     ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],rmodel.int.lin,
                                                      clampba=T),
                                     t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.int.lin,
                                                    clampba=T)),seq)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

load("./Output/lambda_effects.rda")
FIA <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)
# when running IPM without trees killed by fire, technically should filter out those trees
# prolly makes ~0 difference
FIA <- FIA[!(FIA$DSTRBCD1 %in% c(30, 31, 32, 80)), ]

FIA_lambda<-read.csv("./Output/FIA_lambda.csv")

lambda_c<-raster("./Output/tifs/PIED.clim_lambda.tif")
lambda_ccl<-raster("./Output/tifs/PIED.climclamp_lambda.tif")
lambda_cc<-raster("./Output/tifs/PIED.climcomp_lambda.tif")
lambda_ccf<-raster("./Output/tifs/PIED.climcompfire_lambda.tif")
lambda_i<-raster("./Output/tifs/PIED.int_lambda.tif")
lambda_min<-min(c(minValue(lambda_c),minValue(lambda_ccl),minValue(lambda_cc),minValue(lambda_ccf)
                  ,minValue(lambda_i)))
lambda_max<-max(c(maxValue(lambda_c),maxValue(lambda_ccl),maxValue(lambda_cc),maxValue(lambda_ccf),
                  maxValue(lambda_i)))

grow_c<-raster("./Output/tifs/PIED.clim_growth.tif")
grow_ccl<-raster("./Output/tifs/PIED.climclamp_growth.tif")
grow_cc<-raster("./Output/tifs/PIED.climcomp_growth.tif")
grow_ccf<-raster("./Output/tifs/PIED.climcompfire_growth.tif")
grow_i<-raster("./Output/tifs/PIED.int_growth.tif")
grow_min<-min(c(minValue(grow_c),minValue(grow_ccl),minValue(grow_cc),minValue(grow_ccf)
                ,minValue(grow_i)))
grow_max<-max(c(maxValue(grow_c),maxValue(grow_ccl),maxValue(grow_cc),maxValue(grow_ccf),
                maxValue(grow_i)))

surv_c<-raster("./Output/tifs/PIED.clim_survival.tif")
surv_ccl<-raster("./Output/tifs/PIED.climclamp_survival.tif")
surv_cc<-raster("./Output/tifs/PIED.climcomp_survival.tif")
surv_ccf<-raster("./Output/tifs/PIED.climcompfire_survival.tif")
surv_i<-raster("./Output/tifs/PIED.int_survival.tif")
surv_min<-min(c(minValue(surv_c),minValue(surv_ccl),minValue(surv_cc),minValue(surv_ccf)
                ,minValue(surv_i)))
surv_max<-max(c(maxValue(surv_c),maxValue(surv_ccl),maxValue(surv_cc),maxValue(surv_ccf),
                maxValue(surv_i)))

repr_c<-raster("./Output/tifs/PIED.clim_reproduction.tif")
repr_ccl<-raster("./Output/tifs/PIED.climclamp_reproduction.tif")
repr_cc<-raster("./Output/tifs/PIED.climcomp_reproduction.tif")
repr_ccf<-raster("./Output/tifs/PIED.climcompfire_reproduction.tif")
repr_i<-raster("./Output/tifs/PIED.int_reproduction.tif")
repr_min<-min(c(minValue(repr_c),minValue(repr_ccl),minValue(repr_cc),minValue(repr_ccf)
                ,minValue(repr_i)))
repr_max<-max(c(maxValue(repr_c),maxValue(repr_ccl),maxValue(repr_cc),maxValue(repr_ccf),
                maxValue(repr_i)))

extrap<-raster("./Output/tifs/extrap.tif")

FIA_pa <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)
FIA_pa.plot<-FIA_pa %>%
  group_by(plot) %>%
  summarise(lat=mean(lat),lon=mean(lon),PApied=mean(PApied))
FIA_pa.plot$PApied<-as.factor(FIA_pa.plot$PApied)

#Make points spatial and extract lambda for each point
Points <- SpatialPoints(coords = cbind(FIA_pa.plot$lon, FIA_pa.plot$lat), 
                        proj4string = CRS("+proj=longlat +datum=NAD83"))
FIA_pa.plot$lambda_c <- raster::extract(lambda_c, Points) 
FIA_pa.plot$lambda_ccl <- raster::extract(lambda_ccl, Points) 
FIA_pa.plot$lambda_cc <- raster::extract(lambda_cc, Points) 
FIA_pa.plot$lambda_ccf <- raster::extract(lambda_ccf, Points) 
FIA_pa.plot$lambda_i <- raster::extract(lambda_i, Points) 

l_means<-FIA_pa.plot %>%
  group_by(PApied) %>%
  summarise(lambda_c=median(lambda_c),lambda_ccl=median(lambda_ccl),lambda_cc=median(lambda_cc),
            lambda_ccf=median(lambda_ccf),lambda_i=median(lambda_i))

## Climate only
# Growth
grow_c_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_clim$dia),xmax=min(grdata$PREVDIA),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$dia),xmin=max(grdata$PREVDIA),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=-1.6),size=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_clim,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_clim,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  labs(x="Previous diameter", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_c_d.png", plot=grow_c_d,dpi=400)

grow_c_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=0),size=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_clim,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_clim,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_c_p.png", plot=grow_c_p,dpi=400)

grow_c_t<-ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim$t),xmax=min(grdata$T_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$t),xmin=max(grdata$T_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=-0.005),size=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_clim,aes(x=t,y=t_pred,linetype="No clamp"),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_clim,aes(x=t,y=t_pred_c,linetype="Clamp"),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment")+
  scale_linetype_manual("",breaks=c("No clamp","Clamp"),
                        values=c("No clamp"="solid","Clamp"="dotted"),
                        labels=c("No clamp","Clamp"))+
  theme(legend.key.width=unit(2,"cm"),legend.text=element_text(size=24),
        legend.key.size = unit(3, 'lines')) + 
  theme(legend.position = "none")

ggsave(file="all_PIED_manuscript_grow_c_t.png", plot=grow_c_t,dpi=400)

# Survival
surv_c_d <- ggplot(data=survData,aes(x=PREVDIA,y=mort))+
  geom_rect(aes(xmin=min(splot_data_clim$dia),xmax=min(survData$PREVDIA),ymin=-0.02,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$dia),xmin=max(survData$PREVDIA),ymin=-0.02,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=PREVDIA,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_clim,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_clim,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_c_d.png", plot=surv_c_d,dpi=400)

surv_c_p <- ggplot(data=survData,aes(x=PPT_yr_norm,y=mort))+
  geom_rect(aes(xmin=min(splot_data_clim$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=PPT_yr_norm,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_clim,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_clim,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_c_p.png", plot=surv_c_p,dpi=400)

surv_c_t <- ggplot(data=survData,aes(x=T_yr_norm,y=mort))+
  geom_rect(aes(xmin=min(splot_data_clim$t),xmax=min(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$t),xmin=max(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=T_yr_norm,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_clim,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_clim,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_c_t.png", plot=surv_c_t,dpi=400)

# Recruit
recr_c_p <- ggplot(data=rdata,aes(x=PPT_yr_norm,y=recruits1))+
  geom_rect(aes(xmin=min(rplot_data_clim$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_clim$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point()+
  geom_point(aes(y=-0.15),size=0.1)+
  geom_line(data=rplot_data_clim,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_clim,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_c_p.png", plot=recr_c_p,dpi=400)

recr_c_t <- ggplot(data=rdata,aes(x=T_yr_norm,y=recruits1))+
  geom_rect(aes(xmin=min(rplot_data_clim$t),xmax=min(rdata$T_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_clim$t),xmin=max(rdata$T_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point()+
  geom_point(aes(y=-0.15),size=0.1)+
  geom_line(data=rplot_data_clim,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_clim,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_c_t.png", plot=recr_c_p,dpi=400)

## Climate + competition

# Growth
grow_cc_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$dia),xmax=min(grdata$PREVDIA),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp$dia),xmin=max(grdata$PREVDIA),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_climcomp,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_climcomp,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  labs(x="Previous diameter", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_cc_d.png", plot=grow_cc_d,dpi=400)

grow_cc_b <- ggplot(data=grdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$ba),xmax=min(grdata$BALIVE),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_climcomp,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_climcomp,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  labs(x="Live basal area", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_cc_b.png", plot=grow_cc_b,dpi=400)

grow_cc_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_climcomp,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_climcomp,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_cc_p.png", plot=grow_cc_p,dpi=400)

grow_cc_t <- ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$t),xmax=min(grdata$T_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp$t),xmin=max(grdata$T_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_climcomp,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_climcomp,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_cc_t.png", plot=grow_cc_t,dpi=400)

# Survival
surv_cc_d <- ggplot(data=survData,aes(x=PREVDIA,y=mort))+
  geom_rect(aes(xmin=min(splot_data_climcomp$dia),xmax=min(survData$PREVDIA),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp$dia),xmin=max(survData$PREVDIA),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=PREVDIA,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_climcomp,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_climcomp,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  geom_line(data=splot_data_fire,aes(x=dia,y=dia_pred),col="#d95f02",size=1.25)+
  geom_line(data=splot_data_fire,aes(x=dia,y=dia_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_cc_d.png", plot=surv_cc_d,dpi=400)

surv_cc_b <- ggplot(data=survData,aes(x=BALIVE,y=mort))+
  geom_rect(aes(xmin=min(splot_data_climcomp$ba),xmax=min(survData$BALIVE),ymin=-0.02,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=BALIVE,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_climcomp,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_climcomp,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  geom_line(data=splot_data_fire,aes(x=ba,y=ba_pred),col="#d95f02",size=1.25)+
  geom_line(data=splot_data_fire,aes(x=ba,y=ba_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_cc_b.png", plot=surv_cc_b,dpi=400)

surv_cc_p <- ggplot(data=survData,aes(x=PPT_yr_norm,y=mort))+
  geom_rect(aes(xmin=min(splot_data_climcomp$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=PPT_yr_norm,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_climcomp,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_climcomp,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  geom_line(data=splot_data_fire,aes(x=ppt,y=ppt_pred),col="#d95f02",size=1.25)+
  geom_line(data=splot_data_fire,aes(x=ppt,y=ppt_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_cc_p.png", plot=surv_cc_p,dpi=400)

surv_cc_t <- ggplot(data=survData,aes(x=T_yr_norm,y=mort))+
  geom_rect(aes(xmin=min(splot_data_climcomp$t),xmax=min(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp$t),xmin=max(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=T_yr_norm,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_climcomp,aes(x=t,y=t_pred,linetype="No clamp",col="No fire"),size=1.25)+
  geom_line(data=splot_data_climcomp,aes(x=t,y=t_pred_c,linetype="Clamp",col="No fire"),size=1.25)+
  geom_line(data=splot_data_fire,aes(x=t,y=t_pred,linetype="No clamp",col="Fire"),size=1.25)+
  geom_line(data=splot_data_fire,aes(x=t,y=t_pred_c,linetype="Clamp",col="Fire"),size=1.25)+
  labs(x="30-year temperature norm", y="Mortality")+
  scale_linetype_manual("",breaks=c("No clamp","Clamp"),
                        values=c("No clamp"="solid","Clamp"="dotted"),
                        labels=c("No clamp","Clamp"))+
  scale_colour_manual("",breaks=c("No fire","Fire"),
                      values=c("No fire"="#1b9e77","Fire"="#d95f02"),
                      labels=c("No fire","Fire"))+
  theme(legend.key.width=unit(2,"cm"),legend.text=element_text(size=24),
        legend.key.size = unit(3, 'lines')) + theme(legend.position = "none")

ggsave(file="all_PIED_manuscript_surv_cc_t.png", plot=surv_cc_t,dpi=400)

# Recruit
recr_cc_b <- ggplot(data=rdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(rplot_data_climcomp$ba),xmax=min(rdata$BALIVE),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=recruits1))+
  geom_line(data=rplot_data_climcomp,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_climcomp,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  labs(x="Live basal area", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_cc_b.png", plot=recr_cc_b,dpi=400)

recr_cc_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climcomp$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climcomp$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=recruits1))+
  geom_line(data=rplot_data_climcomp,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_climcomp,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",
            size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_cc_p.png", plot=recr_cc_p,dpi=400)

recr_cc_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climcomp$t),xmax=min(rdata$T_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climcomp$t),xmin=max(rdata$T_yr_norm),
                ymin=min(rdata$recruits1),ymax=max(rdata$recruits1)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=recruits1))+
  geom_line(data=rplot_data_climcomp,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_climcomp,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_cc_t.png", plot=recr_cc_t,dpi=400)

## Climate + competition, interactions

# Growth
grow_i_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_int$dia),xmax=min(grdata$PREVDIA),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int$dia),xmin=max(grdata$PREVDIA),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_int,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_i_d.png", plot=grow_i_d,dpi=400)

grow_i_b <- ggplot(data=grdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(grplot_data_int$ba),xmax=min(grdata$BALIVE),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_int,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_i_b.png", plot=grow_i_b,dpi=400)

grow_i_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_int$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_int,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_i_p.png", plot=grow_i_p,dpi=400)

grow_i_t <- ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_int$t),xmax=min(grdata$T_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int$t),xmin=max(grdata$T_yr_norm),
                ymin=min(grdata$DIA_INCR,na.rm=T),ymax=max(grdata$DIA_INCR,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=DIA_INCR))+
  geom_line(data=grplot_data_int,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment")

ggsave(file="all_PIED_manuscript_grow_i_t.png", plot=grow_i_t,dpi=400)

# Survival
surv_i_d <- ggplot(data=survData,aes(x=PREVDIA,y=mort))+
  geom_rect(aes(xmin=min(splot_data_int$dia),xmax=min(survData$PREVDIA),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int$dia),xmin=max(survData$PREVDIA),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=PREVDIA,y=-0.05),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_int,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_i_d.png", plot=surv_i_d,dpi=400)

surv_i_b <- ggplot(data=survData,aes(x=BALIVE,y=mort))+
  geom_rect(aes(xmin=min(splot_data_int$ba),xmax=min(survData$BALIVE),ymin=-0.02,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=BALIVE,y=-0.08),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_int,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_i_b.png", plot=surv_i_b,dpi=400)

surv_i_p <- ggplot(data=survData,aes(x=PPT_yr_norm,y=mort))+
  geom_rect(aes(xmin=min(splot_data_int$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=PPT_yr_norm,y=-0.08),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_i_p.png", plot=surv_i_p,dpi=400)

surv_i_t <- ggplot(data=survData,aes(x=T_yr_norm,y=mort))+
  geom_rect(aes(xmin=min(splot_data_int$t),xmax=min(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int$t),xmin=max(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=survData,aes(x=T_yr_norm,y=-0.07),size=0.1)+
  geom_point()+
  geom_line(data=splot_data_int,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Mortality")

ggsave(file="all_PIED_manuscript_surv_i_t.png", plot=surv_i_t,dpi=400)

# Recruit
recr_i_b <- ggplot(data=rdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(rplot_data_int$ba),xmax=min(rdata$BALIVE),
                ymin=min(rdata$recruits1,na.rm=T),ymax=max(rdata$recruits1,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=recruits1))+
  geom_line(data=rplot_data_int,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_i_b.png", plot=recr_i_b,dpi=400)

recr_i_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_int$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=min(rdata$recruits1,na.rm=T),ymax=max(rdata$recruits1,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_int$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=min(rdata$recruits1,na.rm=T),ymax=max(rdata$recruits1,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=recruits1))+
  geom_line(data=rplot_data_int,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_i_p.png", plot=recr_i_p,dpi=400)

recr_i_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_int$t),xmax=min(rdata$T_yr_norm),
                ymin=min(rdata$recruits1,na.rm=T),ymax=max(rdata$recruits1,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_int$t),xmin=max(rdata$T_yr_norm),
                ymin=min(rdata$recruits1,na.rm=T),ymax=max(rdata$recruits1,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(aes(y=recruits1))+
  geom_line(data=rplot_data_int,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits")

ggsave(file="all_PIED_manuscript_recr_i_t.png", plot=recr_i_t,dpi=400)
