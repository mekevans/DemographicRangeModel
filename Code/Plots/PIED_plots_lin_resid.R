
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
survData<-survData[which(!is.na(match(survData$PLT_CN,survData3.scaled$PLT_CN))),]

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
grdata<-grdata[-which(grdata$PLT_CN==40383858010690|grdata$PLT_CN==40423861010690|
                        grdata$PLT_CN==40423872010690|grdata$PLT_CN==40424710010690|
                        grdata$PLT_CN==186092190020004|
                        grdata$PLT_CN==188784045020004|grdata$PLT_CN==188784634020004),]

g_fun<-function(dia,ba,ppt,t,plot,model,clampba=F,clampt=F){
  gdata=data.frame(PREVDIA=dia,BALIVE = ifelse(clampba==T & ba >190, 190, ba),
                   PPT_yr_norm = ppt, T_yr_norm = ifelse(clampt==T & t >10.6, 10.6, t))
  scaled.gdata = data.frame(scale(gdata,
                                  scale = gr.scaling$scale[match(names(gdata), 
                                                                 names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), 
                                                                   names(gr.scaling$center))]))  
  scaled.gdata$PLT_CN = plot
  return(predict(model, newdata = scaled.gdata))
}
g_fun_mean<-function(dia,ba,ppt,t,model,clampba=F,clampt=F){
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

#Calculate residuals
grdata$resid_cl<-grdata$DIA_INCR-g_fun(grdata$PREVDIA,grdata$BALIVE,grdata$PPT_yr_norm,grdata$T_yr_norm,
                                      plot=grdata$PLT_CN,gmodel.clim.lin)
grdata$resid_ccl<-grdata$DIA_INCR-g_fun(grdata$PREVDIA,grdata$BALIVE,grdata$PPT_yr_norm,grdata$T_yr_norm,
                                        plot=grdata$PLT_CN,gmodel.clim.comp.lin)
grdata$resid_il<-grdata$DIA_INCR-g_fun(grdata$PREVDIA,grdata$BALIVE,grdata$PPT_yr_norm,grdata$T_yr_norm,
                                       plot=grdata$PLT_CN,gmodel.int.lin)

ncuts=30
chopsize_dia<-cut(grdata$PREVDIA,ncuts)
chopsize_ba<-cut(grdata$BALIVE,ncuts)
chopsize_PPT<-cut(grdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(grdata$T_yr_norm,ncuts)

count_binned_dia<-as.vector(sapply(split(grdata$resid_cl,chopsize_dia),length))
dia_binned<-as.vector(sapply(split(grdata$PREVDIA,chopsize_dia),mean,na.rm=T))
grow_binned_dia_cl<-as.vector(sapply(split(grdata$resid_cl,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.clim.lin))
grow_binned_dia_ccl<-as.vector(sapply(split(grdata$resid_ccl,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.clim.comp.lin))
grow_binned_dia_il<-as.vector(sapply(split(grdata$resid_il,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.int.lin))

count_binned_ba<-as.vector(sapply(split(grdata$resid_cl,chopsize_ba),length))
ba_binned<-as.vector(sapply(split(grdata$BALIVE,chopsize_ba),mean,na.rm=T))
grow_binned_ba_cl<-as.vector(sapply(split(grdata$resid_cl,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.clim.lin))
grow_binned_ba_ccl<-as.vector(sapply(split(grdata$resid_ccl,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.clim.comp.lin))
grow_binned_ba_il<-as.vector(sapply(split(grdata$resid_il,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.int.lin))

count_binned_PPT<-as.vector(sapply(split(grdata$resid_cl,chopsize_PPT),length))
PPT_binned<-as.vector(sapply(split(grdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
grow_binned_PPT_cl<-as.vector(sapply(split(grdata$resid_cl,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.clim.lin))
grow_binned_PPT_ccl<-as.vector(sapply(split(grdata$resid_ccl,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.clim.comp.lin))
grow_binned_PPT_il<-as.vector(sapply(split(grdata$resid_il,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.int.lin))

count_binned_T<-as.vector(sapply(split(grdata$resid_cl,chopsize_T),length))
T_binned<-as.vector(sapply(split(grdata$T_yr_norm,chopsize_T),mean,na.rm=T))
grow_binned_T_cl<-as.vector(sapply(split(grdata$resid_cl,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.clim.lin))
grow_binned_T_ccl<-as.vector(sapply(split(grdata$resid_ccl,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.clim.comp.lin))
grow_binned_T_il<-as.vector(sapply(split(grdata$resid_il,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.int.lin))

g_binned<-as.data.frame(cbind(grow_binned_dia_cl,grow_binned_ba_cl,grow_binned_PPT_cl,grow_binned_T_cl,
                              grow_binned_dia_ccl,grow_binned_ba_ccl,grow_binned_PPT_ccl,grow_binned_T_ccl,
                              grow_binned_dia_il,grow_binned_ba_il,grow_binned_PPT_il,grow_binned_T_il,
                              dia_binned,ba_binned,PPT_binned,T_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))
names(g_binned)<-c("grow_dia_cl","grow_ba_cl","grow_PPT_cl","grow_T_cl",
                   "grow_dia_ccl","grow_ba_ccl","grow_PPT_ccl","grow_T_ccl",
                   "grow_dia_il","grow_ba_il","grow_PPT_il","grow_T_il",
                   "PREVDIA","BALIVE","PPT","T",
                   "count_dia","count_ba","count_PPT","count_T")

grplot_data_clim.lin<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.lin),
                                   ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.lin),
                                   t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.lin),
                                   dia_pred_c=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.lin,
                                                    clampba=T,clampt=T),
                                   ppt_pred_c=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.lin,
                                                    clampba=T,clampt=T),
                                   t_pred_c=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.lin,
                                                  clampba=T,clampt=T)),seq)
grplot_data_climcomp.lin<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],
                                                      gmodel.clim.comp.lin),
                                       ba_pred=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.clim.comp.lin),
                                       ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.comp.lin),
                                       t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.comp.lin),
                                       dia_pred_c=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.comp.lin,
                                                        clampba=T,clampt=F),
                                       ba_pred_c=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.clim.comp.lin,
                                                       clampba=T,clampt=F),
                                       ppt_pred_c=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.comp.lin,
                                                        clampba=T,clampt=F),
                                       t_pred_c=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.comp.lin,
                                                      clampba=T,clampt=F)),seq)
grplot_data_int.lin<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],
                                                 gmodel.int.lin),
                                  ba_pred=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.int.lin),
                                  ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.int.lin),
                                  t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.int.lin),
                                  dia_pred_c=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.int.lin,
                                                   clampba=T,clampt=F),
                                  ba_pred_c=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.int.lin,
                                                  clampba=T,clampt=F),
                                  ppt_pred_c=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.int.lin,
                                                   clampba=T,clampt=F),
                                  t_pred_c=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.int.lin,
                                                 clampba=T,clampt=F)),seq)

#Survival
s_fun<-function(dia,ba,ppt,t,ci,plot,model,clampt=F){
  sdata=data.frame(PREVDIA=dia,BALIVE=ba,PPT_yr_norm = ppt, 
                   T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), 
                                                                   names(surv.scaling$scale))], 
                                  center = surv.scaling$center[match(names(sdata), 
                                                                     names(surv.scaling$center))])) 
  scaled.sdata = cbind(scaled.sdata,CENSUS_INTERVAL = ci)
  scaled.sdata$PLT_CN=plot
  return(predict(model, newdata = scaled.sdata, type = "response"))
}
s_fun_mean<-function(dia,ba,ppt,t,ci,model,clampt=F){
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

#Calculate residuals
survData$resid_c_lin<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                          survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                          survData$PLT_CN,smodel.clim.lin)
survData$resid_cc_lin<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                           survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                           survData$PLT_CN,smodel.clim.comp.lin)
survData$resid_ccf_lin<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                            survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                            survData$PLT_CN,smodel.clim.comp.fire.lin)
survData$resid_i_lin<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                          survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                          survData$PLT_CN,smodel.int.lin)

ncuts=50
chopsize_dia<-cut(survData$PREVDIA,ncuts)
chopsize_ba<-cut(survData$BALIVE,ncuts)
chopsize_PPT<-cut(survData$PPT_yr_norm,ncuts)
chopsize_T<-cut(survData$T_yr_norm,ncuts)

dia_binned<-as.vector(sapply(split(survData$PREVDIA,chopsize_dia),mean,na.rm=T))
count_binned_dia<-as.vector(sapply(split(survData$mort,chopsize_dia),length))
surv_binned_dia_cl<-as.vector(sapply(split(survData$resid_c_lin,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.lin))
surv_binned_dia_ccl<-as.vector(sapply(split(survData$resid_cc_lin,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.comp.lin))
surv_binned_dia_ccfl<-as.vector(sapply(split(survData$resid_ccf_lin,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.comp.fire.lin))
surv_binned_dia_il<-as.vector(sapply(split(survData$resid_i_lin,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.int.lin))

ba_binned<-as.vector(sapply(split(survData$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(survData$mort,chopsize_ba),length))
surv_binned_ba_cl<-as.vector(sapply(split(survData$resid_c_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.lin))
surv_binned_ba_ccl<-as.vector(sapply(split(survData$resid_cc_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.comp.lin))
surv_binned_ba_ccfl<-as.vector(sapply(split(survData$resid_ccf_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.comp.fire.lin))
surv_binned_ba_il<-as.vector(sapply(split(survData$resid_i_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.int.lin))

PPT_binned<-as.vector(sapply(split(survData$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(survData$mort,chopsize_PPT),length))
surv_binned_PPT_cl<-as.vector(sapply(split(survData$resid_c_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.lin))
surv_binned_PPT_ccl<-as.vector(sapply(split(survData$resid_cc_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.comp.lin))
surv_binned_PPT_ccfl<-as.vector(sapply(split(survData$resid_ccf_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.comp.fire.lin))
surv_binned_PPT_il<-as.vector(sapply(split(survData$resid_i_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.int.lin))

T_binned<-as.vector(sapply(split(survData$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(survData$mort,chopsize_T),length))
surv_binned_T_cl<-as.vector(sapply(split(survData$resid_c_lin,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.clim.lin))
surv_binned_T_ccl<-as.vector(sapply(split(survData$resid_cc_lin,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.clim.comp.lin))
surv_binned_T_ccfl<-as.vector(sapply(split(survData$resid_ccf_lin,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.clim.comp.fire.lin))
surv_binned_T_il<-as.vector(sapply(split(survData$resid_i_lin,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.int.lin))

s_binned<-as.data.frame(cbind(surv_binned_dia_cl,surv_binned_ba_cl,surv_binned_PPT_cl,surv_binned_T_cl,
                              surv_binned_dia_ccl,surv_binned_ba_ccl,surv_binned_PPT_ccl,surv_binned_T_ccl,
                              surv_binned_dia_ccfl,surv_binned_ba_ccfl,surv_binned_PPT_ccfl,surv_binned_T_ccfl,
                              surv_binned_dia_il,surv_binned_ba_il,surv_binned_PPT_il,surv_binned_T_il,
                              dia_binned,ba_binned,PPT_binned,T_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))

names(s_binned)<-c("mort_dia_cl","mort_ba_cl","mort_PPT_cl","mort_T_cl",
                   "mort_dia_ccl","mort_ba_ccl","mort_PPT_ccl","mort_T_ccl",
                   "mort_dia_ccfl","mort_ba_ccfl","mort_PPT_ccfl","mort_T_ccfl",
                   "mort_dia_il","mort_ba_il","mort_PPT_il","mort_T_il",
                   "PREVDIA","BALIVE","PPT","T",
                   "count_dia","count_ba","count_PPT","count_T")

splot_data_clim.lin<-cbind(data.frame(dia_pred=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.clim.lin),
                                      ppt_pred=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim.lin),
                                      t_pred=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.clim.lin),
                                      dia_pred_c=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                       smodel.clim.lin,clampt=T),
                                      ppt_pred_c=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                       smodel.clim.lin,clampt=T),
                                      t_pred_c=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.clim.lin,
                                                     clampt=T)),seq)
splot_data_climcomp.lin<-cbind(data.frame(dia_pred=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                         smodel.clim.comp.lin),
                                          ba_pred=s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],
                                                        smodel.clim.comp.lin),
                                          ppt_pred=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                         smodel.clim.comp.lin),
                                          t_pred=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],
                                                       smodel.clim.comp.lin),
                                          dia_pred_c=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                           smodel.clim.comp.lin,clampt=T),
                                          ba_pred_c=s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],
                                                          smodel.clim.comp.lin,clampt=T),
                                          ppt_pred_c=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                           smodel.clim.comp.lin,clampt=T),
                                          t_pred_c=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],
                                                         smodel.clim.comp.lin,clampt=T)),seq)
splot_data_fire.lin<-cbind(data.frame(dia_pred=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.clim.comp.fire.lin),
                                      ba_pred=s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],
                                                    smodel.clim.comp.fire.lin),
                                      ppt_pred=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                     smodel.clim.comp.fire.lin),
                                      t_pred=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],
                                                   smodel.clim.comp.fire.lin),
                                      dia_pred_c=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                       smodel.clim.comp.fire.lin,clampt=T),
                                      ba_pred_c=s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],
                                                      smodel.clim.comp.fire.lin,clampt=T),
                                      ppt_pred_c=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                       smodel.clim.comp.fire.lin,clampt=T),
                                      t_pred_c=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],
                                                     smodel.clim.comp.fire.lin,clampt=T)),seq)
splot_data_int.lin<-cbind(data.frame(dia_pred=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                    smodel.int.lin),
                                     ba_pred=s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],smodel.int.lin),
                                     ppt_pred=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.lin),
                                     t_pred=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.int.lin),
                                     dia_pred_c=s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],smodel.int.lin,
                                                      clampt=T),
                                     ba_pred_c=s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],smodel.int.lin,
                                                     clampt=T),
                                     ppt_pred_c=s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.lin,
                                                      clampt=T),
                                     t_pred_c=s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.int.lin,
                                                    clampt=T)),seq)

#Recruitment
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

#Calculate residuals
rdata$resid_c_lin<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                         rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.clim.lin)
rdata$resid_cc_lin<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                          rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.clim.comp.lin)
rdata$resid_i_lin<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                         rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.int.lin)

ncuts=30
chopsize_ba<-cut(rdata$BALIVE,ncuts)
chopsize_PPT<-cut(rdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(rdata$T_yr_norm,ncuts)

ba_binned<-as.vector(sapply(split(rdata$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(rdata$recruits1,chopsize_ba),length))
recr_binned_ba_cl<-as.vector(sapply(split(rdata$resid_c_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.clim.lin))
recr_binned_ba_ccl<-as.vector(sapply(split(rdata$resid_cc_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.clim.comp.lin))
recr_binned_ba_il<-as.vector(sapply(split(rdata$resid_i_lin,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.int.lin))

PPT_binned<-as.vector(sapply(split(rdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(rdata$recruits1,chopsize_PPT),length))
recr_binned_PPT_cl<-as.vector(sapply(split(rdata$resid_c_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.clim.lin))
recr_binned_PPT_ccl<-as.vector(sapply(split(rdata$resid_cc_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.clim.comp.lin))
recr_binned_PPT_il<-as.vector(sapply(split(rdata$resid_i_lin,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.int.lin))

T_binned<-as.vector(sapply(split(rdata$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(rdata$recruits1,chopsize_T),length))
recr_binned_T_cl<-as.vector(sapply(split(rdata$resid_c_lin,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.clim.lin))
recr_binned_T_ccl<-as.vector(sapply(split(rdata$resid_cc_lin,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.clim.comp.lin))
recr_binned_T_il<-as.vector(sapply(split(rdata$resid_i_lin,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.int.lin))

r_binned<-as.data.frame(cbind(recr_binned_ba_cl,recr_binned_PPT_cl,recr_binned_T_cl,
                              recr_binned_ba_ccl,recr_binned_PPT_ccl,recr_binned_T_ccl,
                              recr_binned_ba_il,recr_binned_PPT_il,recr_binned_T_il,
                              ba_binned,PPT_binned,T_binned,
                              count_binned_ba,count_binned_PPT,count_binned_T))
names(r_binned)<-c("recr_ba_cl","recr_PPT_cl","recr_T_cl",
                   "recr_ba_ccl","recr_PPT_ccl","recr_T_ccl",
                   "recr_ba_il","recr_PPT_il","recr_T_il",
                   "BALIVE","PPT","T","count_ba","count_PPT","count_T")

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


mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               legend.text=element_text(size=11),legend.title=element_text(size=12),
               legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
               axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
               axis.line.x = element_line(color="black", size = 0.3),
               axis.line.y = element_line(color="black", size = 0.3))

## Climate only
# Growth
grow_c_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_clim.lin$dia),xmax=min(grdata$PREVDIA),ymin=-1.5,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim.lin$dia),xmin=max(grdata$PREVDIA),ymin=-1.5,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-1.6),size=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_cl,size=count_dia))+
  geom_line(data=grplot_data_clim.lin,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_clim.lin,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  labs(x="Previous diameter", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_c_d.png", plot=grow_c_d,dpi=400)

grow_c_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim.lin$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=0.015,ymax=max(grplot_data_clim.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim.lin$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=0.005,ymax=max(grplot_data_clim.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=0),size=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_cl,size=count_PPT))+
  geom_line(data=grplot_data_clim.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_clim.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_c_p.png", plot=grow_c_p,dpi=400)

grow_c_t<-ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim.lin$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(g_binned$grow_T_cl,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim.lin$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(g_binned$grow_T_cl,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-0.005),size=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_cl,size=count_T))+
  geom_line(data=grplot_data_clim.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_clim.lin,aes(x=t,y=t_pred_c),linetype="dotted",col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_c_t.png", plot=grow_c_t,dpi=400)

# Survival
surv_c_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_clim.lin$dia),xmax=min(survData$PREVDIA),ymin=-0.02,ymax=0.8),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim.lin$dia),xmin=max(survData$PREVDIA),ymin=-0.02,ymax=0.8),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=PREVDIA,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_cl,size=count_dia))+
  geom_line(data=splot_data_clim.lin,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_clim.lin,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top") + mytheme

ggsave(file="lin_PIED_manuscript_surv_c_d.png", plot=surv_c_d,dpi=400)

surv_c_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_clim.lin$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=0,ymax=0.25),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim.lin$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=0,ymax=0.25),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=PPT_yr_norm,y=-0.01),size=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_cl,size=count_PPT))+
  geom_line(data=splot_data_clim.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_clim.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_c_p.png", plot=surv_c_p,dpi=400)

surv_c_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_clim.lin$t),xmax=min(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim.lin$t),xmin=max(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=T_yr_norm,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_cl,size=count_T))+
  geom_line(data=splot_data_clim.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_clim.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_c_t.png", plot=surv_c_t,dpi=400)

# Recruit
recr_c_p <- ggplot(data=rdata.scaled,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_clim.lin$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_clim.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_clim.lin$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_clim.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_cl,size=count_PPT))+
  #geom_point(aes(y=-0.15),size=0.1)+
  #geom_line(data=rplot_data_clim,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_clim,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  geom_line(data=rplot_data_clim.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_clim.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_c_p.png", plot=recr_c_p,dpi=400)

recr_c_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_clim.lin$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=max(r_binned$recr_T_cl,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_clim.lin$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=max(r_binned$recr_T_cl,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_cl,size=count_T))+
  #geom_point(aes(y=-0.1),size=0.1)+
  #geom_line(data=rplot_data_clim,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_clim,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  geom_line(data=rplot_data_clim.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_clim.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_c_t.png", plot=recr_c_t,dpi=400)

## Climate + competition

# Growth
grow_cc_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_climcomp.lin$dia),xmax=min(grdata$PREVDIA),ymin=-1.5,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp.lin$dia),xmin=max(grdata$PREVDIA),ymin=-1.5,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-1.6),size=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_ccl,size=count_dia))+
  geom_line(data=grplot_data_climcomp.lin,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_climcomp.lin,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  labs(x="Previous diameter", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_cc_d.png", plot=grow_cc_d,dpi=400)

grow_cc_b <- ggplot(data=grdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(grplot_data_climcomp.lin$ba),xmax=min(grdata$BALIVE),
                ymin=min(grplot_data_climcomp.lin$ba_pred),ymax=0.1),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-0.03),size=0.1)+
  geom_point(data=g_binned,aes(x=BALIVE,y=grow_ba_ccl,size=count_ba))+
  geom_line(data=grplot_data_climcomp.lin,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_climcomp.lin,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  labs(x="Live basal area", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_cc_b.png", plot=grow_cc_b,dpi=400)

grow_cc_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climcomp.lin$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=0,ymax=max(grplot_data_climcomp.lin$ppt_pred)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp.lin$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=0,ymax=max(grplot_data_climcomp.lin$ppt_pred)),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=0),size=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_ccl,size=count_PPT))+
  geom_line(data=grplot_data_climcomp.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_climcomp.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_cc_p.png", plot=grow_cc_p,dpi=400)

grow_cc_t <- ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climcomp.lin$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_climcomp.lin$t_pred,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp.lin$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_climcomp.lin$t_pred,na.rm=T)),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-0.005),size=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_ccl,size=count_T))+
  geom_line(data=grplot_data_climcomp.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_climcomp.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_cc_t.png", plot=grow_cc_t,dpi=400)

# Survival
surv_cc_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_climcomp.lin$dia),xmax=min(survData$PREVDIA),
                ymin=-0.02,ymax=max(splot_data_climcomp.lin$dia_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp.lin$dia),xmin=max(survData$PREVDIA),
                ymin=-0.02,ymax=max(splot_data_climcomp.lin$dia_pred)),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=PREVDIA,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_ccl,size=count_dia))+
  geom_line(data=splot_data_climcomp.lin,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_climcomp.lin,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=dia,y=dia_pred),col="#d95f02",size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=dia,y=dia_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_cc_d.png", plot=surv_cc_d,dpi=400)

surv_cc_b <- ggplot(data=survData,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(splot_data_climcomp.lin$ba),xmax=min(survData$BALIVE),ymin=-0.02,ymax=0.4),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=BALIVE,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=BALIVE,y=mort_ba_ccl,size=count_ba))+
  geom_line(data=splot_data_climcomp.lin,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_climcomp.lin,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=ba,y=ba_pred),col="#d95f02",size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=ba,y=ba_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_cc_b.png", plot=surv_cc_b,dpi=400)

surv_cc_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_climcomp.lin$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=-0.02,ymax=0.4),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp.lin$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=-0.02,ymax=0.4),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=PPT_yr_norm,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_ccl,size=count_PPT))+
  geom_line(data=splot_data_climcomp.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_climcomp.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=ppt,y=ppt_pred),col="#d95f02",size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=ppt,y=ppt_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_cc_p.png", plot=surv_cc_p,dpi=400)

surv_cc_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_climcomp.lin$t),xmax=min(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp.lin$t),xmin=max(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=T_yr_norm,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_ccl,size=count_T))+
  geom_line(data=splot_data_climcomp.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_climcomp.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=t,y=t_pred),col="#d95f02",size=1.25)+
  #geom_line(data=splot_data_fire.lin,aes(x=t,y=t_pred_c),col="#d95f02",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_cc_t.png", plot=surv_cc_t,dpi=400)

# Recruit
recr_cc_b <- ggplot(data=rdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(rplot_data_climcomp.lin$ba),xmax=min(rdata$BALIVE),
                ymin=-0.05,ymax=2),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=BALIVE,y=recr_ba_ccl,size=count_ba))+
  #geom_point(aes(y=-0.2),size=0.1)+
  geom_line(data=rplot_data_climcomp.lin,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_climcomp.lin,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  labs(x="Live basal area", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_cc_b.png", plot=recr_cc_b,dpi=400)

recr_cc_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climcomp.lin$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_climcomp.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climcomp.lin$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_climcomp.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_ccl,size=count_PPT))+
  #geom_point(aes(y=-0.15),size=0.1)+
  geom_line(data=rplot_data_climcomp.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_climcomp.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",
  #          size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_cc_p.png", plot=recr_cc_p,dpi=400)

recr_cc_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climcomp.lin$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=0.8),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climcomp.lin$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=0.8),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_ccl,size=count_T))+
  #geom_point(aes(y=-0.1),size=0.1)+
  geom_line(data=rplot_data_climcomp.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_climcomp.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_cc_t.png", plot=recr_cc_t,dpi=400)

## Climate + competition, interactions

# Growth
grow_i_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_int.lin$dia),xmax=min(grdata$PREVDIA),ymin=-1.53,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int.lin$dia),xmin=max(grdata$PREVDIA),ymin=-1.53,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-1.6),size=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_il,size=count_dia))+
  geom_line(data=grplot_data_int.lin,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_int.lin,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_i_d.png", plot=grow_i_d,dpi=400)

grow_i_b <- ggplot(data=grdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(grplot_data_int.lin$ba),xmax=min(grdata$BALIVE),
                ymin=-0.03,ymax=0.1),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-0.035),size=0.1)+
  geom_point(data=g_binned,aes(x=BALIVE,y=grow_ba_il,size=count_ba))+
  geom_line(data=grplot_data_int.lin,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_int.lin,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_i_b.png", plot=grow_i_b,dpi=400)

grow_i_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_int.lin$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=0.015,ymax=max(grplot_data_int.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int.lin$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=0.015,ymax=max(grplot_data_int.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=0.01),size=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_il,size=count_PPT))+
  geom_line(data=grplot_data_int.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_int.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_i_p.png", plot=grow_i_p,dpi=400)

grow_i_t <- ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_int.lin$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_int.lin$t_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int.lin$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_int.lin$t_pred)),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(aes(y=-0.005),size=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_il,size=count_T))+
  geom_line(data=grplot_data_int.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=grplot_data_int.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_grow_i_t.png", plot=grow_i_t,dpi=400)

# Survival
surv_i_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_int.lin$dia),xmax=min(survData$PREVDIA),
                ymin=-0.02,ymax=max(splot_data_int.lin$dia_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int.lin$dia),xmin=max(survData$PREVDIA),
                ymin=-0.02,ymax=max(splot_data_int.lin$dia_pred)),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=PREVDIA,y=-0.05),size=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_il,size=count_dia))+
  geom_line(data=splot_data_int.lin,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_int.lin,aes(x=dia,y=dia_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Previous diameter", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_i_d.png", plot=surv_i_d,dpi=400)

surv_i_b <- ggplot(data=survData,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(splot_data_int.lin$ba),xmax=min(survData$BALIVE),ymin=-0.02,ymax=0.4),
            fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=BALIVE,y=-0.08),size=0.1)+
  geom_point(data=s_binned,aes(x=BALIVE,y=mort_ba_il,size=count_ba))+
  geom_line(data=splot_data_int.lin,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_int.lin,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_i_b.png", plot=surv_i_b,dpi=400)

surv_i_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_int.lin$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=-0.02,ymax=0.4),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int.lin$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=-0.02,ymax=0.4),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=PPT_yr_norm,y=-0.08),size=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_il,size=count_PPT))+
  #geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  geom_line(data=splot_data_int.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_int.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_i_p.png", plot=surv_i_p,dpi=400)

surv_i_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_int.lin$t),xmax=min(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int.lin$t),xmin=max(survData$T_yr_norm),
                ymin=-0.02,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  #geom_point(data=survData,aes(x=T_yr_norm,y=-0.07),size=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_il,size=count_T))+
  geom_line(data=splot_data_int.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=splot_data_int.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Mortality")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_surv_i_t.png", plot=surv_i_t,dpi=400)

# Recruit
recr_i_b <- ggplot(data=rdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(rplot_data_int.lin$ba),xmax=min(rdata$BALIVE),
                ymin=-0.4,ymax=3),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=BALIVE,y=recr_ba_il,size=count_ba))+
  #geom_point(aes(y=-1),size=0.1)+
  geom_line(data=rplot_data_int.lin,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_int.lin,aes(x=ba,y=ba_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="Live basal area", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_i_b.png", plot=recr_i_b,dpi=400)

recr_i_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_int.lin$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_int.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_int.lin$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_int.lin$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_il,size=count_PPT))+
  #geom_point(aes(y=-0.15),size=0.1)+
  geom_line(data=rplot_data_int.lin,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_int.lin,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_i_p.png", plot=recr_i_p,dpi=400)

recr_i_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_int.lin$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_int.lin$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_il,size=count_T))+
  #geom_point(aes(y=-0.1),size=0.1)+
  geom_line(data=rplot_data_int.lin,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  #geom_line(data=rplot_data_int.lin,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top")+mytheme

ggsave(file="lin_PIED_manuscript_recr_i_t.png", plot=recr_i_t,dpi=400)


## Partial residual plots

# Climate only
est<-Effect("PREVDIA", partial.residuals=T, gmodel.clim)
grow_c_d_resid<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, gmodel.clim)
grow_c_p_resid<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, gmodel.clim)
grow_c_t_resid<-plot(est)

est<-Effect("PREVDIA", partial.residuals=T, smodel.clim.lin)
surv_c_d_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, smodel.clim.lin)
surv_c_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, smodel.clim.lin)
surv_c_t_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, rmodel.clim.lin)
recr_c_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, rmodel.clim.lin)
recr_c_t_resid_lin<-plot(est)

# Climate + comp, no interactions
est<-Effect("PREVDIA", partial.residuals=T, gmodel.clim.comp)
grow_cc_d_resid<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, gmodel.clim.comp)
grow_cc_b_resid<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, gmodel.clim.comp)
grow_cc_p_resid<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, gmodel.clim.comp)
grow_cc_t_resid<-plot(est)

est<-Effect("PREVDIA", partial.residuals=T, smodel.clim.comp.lin)
surv_cc_d_resid_lin<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, smodel.clim.comp.lin)
surv_cc_b_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, smodel.clim.comp.lin)
surv_cc_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, smodel.clim.comp.lin)
surv_cc_t_resid_lin<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, rmodel.clim.comp.lin)
recr_cc_b_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, rmodel.clim.comp.lin)
recr_cc_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, rmodel.clim.comp.lin)
recr_cc_t_resid_lin<-plot(est)

# Climat + comp, no interactions, fire
est<-Effect("PREVDIA", partial.residuals=T, smodel.clim.comp.fire.lin)
surv_ccf_d_resid_lin<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, smodel.clim.comp.fire.lin)
surv_ccf_d_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, smodel.clim.comp.fire.lin)
surv_ccf_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, smodel.clim.comp.fire.lin)
surv_ccf_t_resid_lin<-plot(est)

# Climate + comp, interactions
est<-Effect("PREVDIA", partial.residuals=T, gmodel.int)
grow_i_d_resid<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, gmodel.int)
grow_i_b_resid<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, gmodel.int)
grow_i_p_resid<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, gmodel.int)
grow_i_t_resid<-plot(est)

est<-Effect("PREVDIA", partial.residuals=T, smodel.int.lin)
surv_i_d_resid_lin<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, smodel.int.lin)
surv_i_b_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, smodel.int.lin)
surv_i_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, smodel.int.lin)
surv_i_t_resid_lin<-plot(est)

est<-Effect("BALIVE", partial.residuals=T, rmodel.int.lin)
recr_i_b_resid_lin<-plot(est)

est<-Effect("PPT_yr_norm", partial.residuals=T, rmodel.int.lin)
recr_i_p_resid_lin<-plot(est)

est<-Effect("T_yr_norm", partial.residuals=T, rmodel.int.lin)
recr_i_t_resid_lin<-plot(est)
