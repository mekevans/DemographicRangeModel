
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

load("./Code/IPM/GrRescaling_gam.Rdata")
load("./Code/IPM/SurvRescaling_gam.Rdata")
#load("./Code/IPM/SurvRescalingFire_gam.Rdata")
load("./Code/IPM/RecruitRescaling_gam.Rdata")
load("./Code/IPM/recrstats.rda")

ppt_yr_raster <- raster("./ClimateData/PRISM/Normals/PPT_year.tif")
t_yr_raster <- raster("./ClimateData/PRISM/Normals/T_year.tif")
ba_raster <- raster("./BA/balive_RF.tif")

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
grdata<-grdata[which(!is.na(match(grdata$PLT_CN,grdata.scaled$PLT_CN))),]
grdata=grdata[-which(grdata$DIA_INCR<(-1)),]

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

#Growth
grdata$PLT_CN_factor<-as.factor(grdata$PLT_CN)
grdata<-grdata[-which(grdata$PLT_CN==40383858010690|grdata$PLT_CN==40423861010690|
                        grdata$PLT_CN==40423872010690|grdata$PLT_CN==40424710010690|
                        grdata$PLT_CN==186092190020004|
                        grdata$PLT_CN==188784045020004|grdata$PLT_CN==188784634020004),]

g_fun<-function(dia,ba,ppt,t,plot,model,clampba=F,clampt=F){
  gdata=data.frame(PREVDIA=dia,BALIVE = ba, #ifelse(clampba==T & ba >190, 190, ba),
                   #PPT_yr_norm = ppt, T_yr_norm = t #ifelse(clampt==T & t >10.6, 10.6, t),
                   PPT_yr_norm = ppt, T_yr_norm = t #ifelse(clampt==T & t >10.6, 10.6, t),
  )
  scaled.gdata = data.frame(scale(gdata,
                                  scale = gr.scaling$scale[match(names(gdata), 
                                                                 names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), 
                                                                   names(gr.scaling$center))]))
  scaled.gdata$PLT_CN_factor = plot
  return(predict(model, newdata = scaled.gdata, re.form = NA))
}

g_fun_mean<-function(dia,ba,ppt,t,model,clampba=F,clampt=F){
  gdata=data.frame(PREVDIA=dia,BALIVE = ba, #ifelse(clampba==T & ba >190, 190, ba),
                   #PPT_yr_norm = ppt, T_yr_norm = t #ifelse(clampt==T & t >10.6, 10.6, t),
                   PPT_yr_norm = ppt, T_yr_norm = t #ifelse(clampt==T & t >10.6, 10.6, t),
  )
  scaled.gdata = data.frame(scale(gdata,
                                  scale = gr.scaling$scale[match(names(gdata), 
                                                                 names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), 
                                                                   names(gr.scaling$center))]))
  scaled.gdata$PLT_CN_factor = 14546600020004
  return(predict(model, newdata = scaled.gdata, re.form = NA, exclude = "s(PLT_CN_factor)"))
}

means<-c(mean(grdata$PREVDIA,na.rm=T),mean(grdata$BALIVE,na.rm=T),
         mean(grdata$PPT_yr_norm,na.rm=T),mean(grdata$T_yr_norm,na.rm=T))
twentyfive<-c(quantile(grdata$PREVDIA,na.rm=T)[2],quantile(grdata$BALIVE,na.rm=T)[2],
         quantile(grdata$PPT_yr_norm,na.rm=T)[2],quantile(grdata$T_yr_norm,na.rm=T)[2])
seventyfive<-c(quantile(grdata$PREVDIA,na.rm=T)[4],quantile(grdata$BALIVE,na.rm=T)[4],
               quantile(grdata$PPT_yr_norm,na.rm=T)[4],quantile(grdata$T_yr_norm,na.rm=T)[4])

seq<-data.frame(dia=seq(0.5*min(grdata$PREVDIA),1.2*max(grdata$PREVDIA),length=50),
                ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(grdata$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))

#Calculate residuals
grdata$resid_c<-grdata$DIA_INCR-g_fun(dia=grdata$PREVDIA,ba=grdata$BALIVE,ppt=grdata$PPT_yr_norm,t=grdata$T_yr_norm,
                                      plot=grdata$PLT_CN_factor,model=gmodel.clim.gam)
grdata$resid_ci<-grdata$DIA_INCR-g_fun(grdata$PREVDIA,grdata$BALIVE,grdata$PPT_yr_norm,grdata$T_yr_norm,
                                       plot=grdata$PLT_CN_factor,gmodel.clim.int.gam)
grdata$resid_cc<-grdata$DIA_INCR-g_fun(grdata$PREVDIA,grdata$BALIVE,grdata$PPT_yr_norm,grdata$T_yr_norm,
                                       plot=grdata$PLT_CN_factor,gmodel.clim.comp.gam)
grdata$resid_i<-grdata$DIA_INCR-g_fun(grdata$PREVDIA,grdata$BALIVE,grdata$PPT_yr_norm,grdata$T_yr_norm,
                                      plot=grdata$PLT_CN_factor,gmodel.int.gam)

ncuts=30
chopsize_dia<-cut(grdata$PREVDIA,ncuts)
chopsize_ba<-cut(grdata$BALIVE,ncuts)
chopsize_PPT<-cut(grdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(grdata$T_yr_norm,ncuts)

count_binned_dia<-as.vector(sapply(split(grdata$resid_c,chopsize_dia),length))
dia_binned<-as.vector(sapply(split(grdata$PREVDIA,chopsize_dia),mean,na.rm=T))
grow_binned_dia_c<-as.vector(sapply(split(grdata$resid_c,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.clim.gam))
grow_binned_dia_ci<-as.vector(sapply(split(grdata$resid_ci,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.clim.int.gam))
grow_binned_dia_cc<-as.vector(sapply(split(grdata$resid_cc,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.clim.comp.gam))
grow_binned_dia_i<-as.vector(sapply(split(grdata$resid_i,chopsize_dia),mean,na.rm=T))+
  as.vector(g_fun_mean(dia_binned,means[2],means[3],means[4],gmodel.int.gam))

count_binned_ba<-as.vector(sapply(split(grdata$resid_c,chopsize_ba),length))
ba_binned<-as.vector(sapply(split(grdata$BALIVE,chopsize_ba),mean,na.rm=T))
grow_binned_ba_c<-as.vector(sapply(split(grdata$resid_c,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.clim.gam))
grow_binned_ba_ci<-as.vector(sapply(split(grdata$resid_ci,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.clim.int.gam))
grow_binned_ba_cc<-as.vector(sapply(split(grdata$resid_cc,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.clim.comp.gam))
grow_binned_ba_i<-as.vector(sapply(split(grdata$resid_i,chopsize_ba),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],ba_binned,means[3],means[4],gmodel.int.gam))

count_binned_PPT<-as.vector(sapply(split(grdata$resid_c,chopsize_PPT),length))
PPT_binned<-as.vector(sapply(split(grdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
grow_binned_PPT_c<-as.vector(sapply(split(grdata$resid_c,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.clim.gam))
grow_binned_PPT_ci<-as.vector(sapply(split(grdata$resid_ci,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.clim.int.gam))
grow_binned_PPT_cc<-as.vector(sapply(split(grdata$resid_cc,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.clim.comp.gam))
grow_binned_PPT_i<-as.vector(sapply(split(grdata$resid_i,chopsize_PPT),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],PPT_binned,means[4],gmodel.int.gam))

count_binned_T<-as.vector(sapply(split(grdata$resid_c,chopsize_T),length))
T_binned<-as.vector(sapply(split(grdata$T_yr_norm,chopsize_T),mean,na.rm=T))
grow_binned_T_c<-as.vector(sapply(split(grdata$resid_c,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.clim.gam))
grow_binned_T_ci<-as.vector(sapply(split(grdata$resid_ci,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.clim.int.gam))
grow_binned_T_cc<-as.vector(sapply(split(grdata$resid_cc,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.clim.comp.gam))
grow_binned_T_i<-as.vector(sapply(split(grdata$resid_i,chopsize_T),mean,na.rm=T))+
  as.vector(g_fun_mean(means[1],means[2],means[3],T_binned,gmodel.int.gam))

g_binned<-as.data.frame(cbind(grow_binned_dia_c,grow_binned_ba_c,grow_binned_PPT_c,grow_binned_T_c,
                              grow_binned_dia_ci,grow_binned_ba_ci,grow_binned_PPT_ci,grow_binned_T_ci,
                              grow_binned_dia_cc,grow_binned_ba_cc,grow_binned_PPT_cc,grow_binned_T_cc,
                              grow_binned_dia_i,grow_binned_ba_i,grow_binned_PPT_i,grow_binned_T_i,
                              dia_binned,ba_binned,PPT_binned,T_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))
names(g_binned)<-c("grow_dia_c","grow_ba_c","grow_PPT_c","grow_T_c",
                   "grow_dia_ci","grow_ba_ci","grow_PPT_ci","grow_T_ci",
                   "grow_dia_cc","grow_ba_cc","grow_PPT_cc","grow_T_cc",
                   "grow_dia_i","grow_ba_i","grow_PPT_i","grow_T_i",
                   "PREVDIA","BALIVE","PPT","T",
                   "count_dia","count_ba","count_PPT","count_T")

grplot_data_clim<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.gam),
                                   ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.gam),
                                   t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.gam)),seq)
grplot_data_climint<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.int.gam),
                                   ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.int.gam),
                                   t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.int.gam),
                                   dia_pred_c=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.int.gam,
                                                    clampba=T,clampt=T),
                                   ppt_pred_c=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.int.gam,
                                                    clampba=T,clampt=T),
                                   t_pred_c=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.int.gam,
                                                  clampba=T,clampt=T)),seq)
grplot_data_climcomp<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],
                                                      gmodel.clim.comp.gam),
                                       ba_pred=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.clim.comp.gam),
                                       ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.comp.gam),
                                       t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.comp.gam),
                                       dia_pred_c=g_fun_mean(seq$dia,means[2],means[3],means[4],gmodel.clim.comp.gam,
                                                        clampba=T,clampt=F),
                                       ba_pred_c=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.clim.comp.gam,
                                                       clampba=T,clampt=F),
                                       ppt_pred_c=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.clim.comp.gam,
                                                        clampba=T,clampt=F),
                                       t_pred_c=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.clim.comp.gam,
                                                      clampba=T,clampt=F)),seq)
grplot_data_int<-cbind(data.frame(dia_pred=g_fun_mean(seq$dia,means[2],means[3],means[4],
                                                 gmodel.int.gam),
                                  ba_pred=g_fun_mean(means[1],seq$ba,means[3],means[4],gmodel.int.gam),
                                  ppt_pred=g_fun_mean(means[1],means[2],seq$ppt,means[4],gmodel.int.gam),
                                  t_pred=g_fun_mean(means[1],means[2],means[3],seq$t,gmodel.int.gam),
                                  ppt_pred_25=g_fun_mean(means[1],means[2],seq$ppt,twentyfive[4],gmodel.int.gam),
                                  t_pred_25=g_fun_mean(means[1],means[2],twentyfive[3],seq$t,gmodel.int.gam),
                                  ppt_pred_75=g_fun_mean(means[1],means[2],seq$ppt,seventyfive[4],gmodel.int.gam),
                                  t_pred_75=g_fun_mean(means[1],means[2],seventyfive[3],seq$t,gmodel.int.gam)),seq)

#Survival
#survData$PPT_yr_norm<-survData$PPT_yr
#survData$T_yr_norm<-survData$T_yr
survData$PLT_CN_factor<-as.factor(survData$PLT_CN)

s_fun<-function(dia,ba,ppt,t,ci,plot,model,clampt=F){
  sdata=data.frame(PREVDIA=dia,BALIVE=ba,
                   PPT_yr_norm = ppt,T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  #PPT_yr_norm = ppt,T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), 
                                                                   names(surv.scaling$scale))], 
                                  center = surv.scaling$center[match(names(sdata), 
                                                                     names(surv.scaling$center))])) 
scaled.sdata$PPT_yr_norm=as.matrix(scaled.sdata$PPT_yr_norm)
scaled.sdata$T_yr_norm=as.matrix(scaled.sdata$T_yr_norm)
scaled.sdata = cbind(scaled.sdata,CENSUS_INTERVAL = ci)
  scaled.sdata$PLT_CN_factor = plot
  return(predict(model, newdata = scaled.sdata, type = "response", re.form = NA))
}

s_fun_mean<-function(dia,ba,ppt,t,ci,model,clampt=F){
  sdata=data.frame(PREVDIA=dia,BALIVE=ba,
                   PPT_yr_norm = ppt, T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  #PPT_yr_norm = ppt, T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), 
                                                                   names(surv.scaling$scale))], 
                                  center = surv.scaling$center[match(names(sdata), 
                                                                     names(surv.scaling$center))])) 
  scaled.sdata$PPT_yr_norm=as.matrix(scaled.sdata$PPT_yr_norm)
  scaled.sdata$T_yr_norm=as.matrix(scaled.sdata$T_yr_norm)
  scaled.sdata = cbind(scaled.sdata,CENSUS_INTERVAL = ci)
  scaled.sdata$PLT_CN_factor = 14546600020004
  return(predict(model, newdata = scaled.sdata, type = "response", re.form = NA, exclude = "s(PLT_CN_factor)"))
}

means<-c(mean(survData$PREVDIA,na.rm=T),mean(survData$BALIVE,na.rm=T),
         mean(survData$PPT_yr_norm,na.rm=T),mean(survData$T_yr_norm,na.rm=T),
         mean(survData$CENSUS_INTERVAL))
twentyfive<-c(quantile(survData$PREVDIA,na.rm=T)[2],quantile(survData$BALIVE,na.rm=T)[2],
              quantile(survData$PPT_yr_norm,na.rm=T)[2],quantile(survData$T_yr_norm,na.rm=T)[2])
seventyfive<-c(quantile(survData$PREVDIA,na.rm=T)[4],quantile(survData$BALIVE,na.rm=T)[4],
               quantile(survData$PPT_yr_norm,na.rm=T)[4],quantile(survData$T_yr_norm,na.rm=T)[4])

seq<-data.frame(dia=seq(0.5*min(survData$PREVDIA),1.2*max(survData$PREVDIA),length=50),
                ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(survData$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))

#Calculate residuals
survData$resid_c<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                      survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                      survData$PLT_CN_factor,smodel.clim.gam)
survData$resid_ci<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                       survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                       survData$PLT_CN_factor,smodel.clim.int.gam)
survData$resid_cc<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                       survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                       survData$PLT_CN_factor,smodel.clim.comp.gam)
survData$resid_i<-survData$mort-s_fun(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                      survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                      survData$PLT_CN_factor,smodel.int.gam)

ncuts=30
chopsize_dia<-cut(survData$PREVDIA,ncuts)
chopsize_ba<-cut(survData$BALIVE,ncuts)
chopsize_PPT<-cut(survData$PPT_yr_norm,ncuts)
chopsize_T<-cut(survData$T_yr_norm,ncuts)

dia_binned<-as.vector(sapply(split(survData$PREVDIA,chopsize_dia),mean,na.rm=T))
count_binned_dia<-as.vector(sapply(split(survData$mort,chopsize_dia),length))
surv_binned_dia_c<-1-(as.vector(sapply(split(survData$resid_c,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.gam)))
surv_binned_dia_ci<-1-(as.vector(sapply(split(survData$resid_ci,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.int.gam)))
surv_binned_dia_cc<-1-(as.vector(sapply(split(survData$resid_cc,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.comp.gam)))
surv_binned_dia_i<-1-(as.vector(sapply(split(survData$resid_i,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean(dia_binned,means[2],means[3],means[4],means[5],smodel.int.gam)))

ba_binned<-as.vector(sapply(split(survData$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(survData$mort,chopsize_ba),length))
surv_binned_ba_c<-1-(as.vector(sapply(split(survData$resid_c,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.gam)))
surv_binned_ba_ci<-1-(as.vector(sapply(split(survData$resid_ci,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.int.gam)))
surv_binned_ba_cc<-1-(as.vector(sapply(split(survData$resid_cc,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.comp.gam)))
surv_binned_ba_i<-1-(as.vector(sapply(split(survData$resid_i,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],ba_binned,means[3],means[4],means[5],smodel.int.gam)))

PPT_binned<-as.vector(sapply(split(survData$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(survData$mort,chopsize_PPT),length))
surv_binned_PPT_c<-1-(as.vector(sapply(split(survData$resid_c,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.gam)))
surv_binned_PPT_ci<-1-(as.vector(sapply(split(survData$resid_ci,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.int.gam)))
surv_binned_PPT_cc<-1-(as.vector(sapply(split(survData$resid_cc,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.comp.gam)))
surv_binned_PPT_i<-1-(as.vector(sapply(split(survData$resid_i,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],PPT_binned,means[4],means[5],smodel.int.gam)))

T_binned<-as.vector(sapply(split(survData$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(survData$mort,chopsize_T),length))
surv_binned_T_c<-1-(as.vector(sapply(split(survData$resid_c,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.clim.gam)))
surv_binned_T_ci<-1-(as.vector(sapply(split(survData$resid_ci,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.clim.int.gam)))
surv_binned_T_cc<-1-(as.vector(sapply(split(survData$resid_cc,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.clim.comp.gam)))
surv_binned_T_i<-1-(as.vector(sapply(split(survData$resid_i,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean(means[1],means[2],means[3],T_binned,means[5],smodel.int.gam)))

s_binned<-as.data.frame(cbind(surv_binned_dia_c,surv_binned_ba_c,surv_binned_PPT_c,surv_binned_T_c,
                              surv_binned_dia_ci,surv_binned_ba_ci,surv_binned_PPT_ci,surv_binned_T_ci,
                              surv_binned_dia_cc,surv_binned_ba_cc,surv_binned_PPT_cc,surv_binned_T_cc,
                              surv_binned_dia_i,surv_binned_ba_i,surv_binned_PPT_i,surv_binned_T_i,
                              dia_binned,ba_binned,PPT_binned,T_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))

names(s_binned)<-c("mort_dia_c","mort_ba_c","mort_PPT_c","mort_T_c",
                   "mort_dia_ci","mort_ba_ci","mort_PPT_ci","mort_T_ci",
                   "mort_dia_cc","mort_ba_cc","mort_PPT_cc","mort_T_cc",
                   "mort_dia_i","mort_ba_i","mort_PPT_i","mort_T_i",
                   "PREVDIA","BALIVE","PPT","T",
                   "count_dia","count_ba","count_PPT","count_T")

splot_data_clim<-cbind(data.frame(dia_pred=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                 smodel.clim.gam),
                                  ppt_pred=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim.gam),
                                  t_pred=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.clim.gam),
                                  dia_pred_c=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                   smodel.clim.gam,clampt=T),
                                  ppt_pred_c=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                   smodel.clim.gam,clampt=T),
                                  t_pred_c=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.clim.gam,
                                                 clampt=T)),seq)
splot_data_climint<-cbind(data.frame(dia_pred=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                 smodel.clim.int.gam),
                                  ppt_pred=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim.int.gam),
                                  t_pred=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.clim.int.gam),
                                  dia_pred_c=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                   smodel.clim.int.gam,clampt=T),
                                  ppt_pred_c=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                   smodel.clim.int.gam,clampt=T),
                                  t_pred_c=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.clim.int.gam,
                                                 clampt=T)),seq)
splot_data_climcomp<-cbind(data.frame(dia_pred=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.clim.comp.gam),
                                      ba_pred=1-s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],
                                                    smodel.clim.comp.gam),
                                      ppt_pred=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                     smodel.clim.comp.gam),
                                      t_pred=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],
                                                   smodel.clim.comp.gam),
                                      dia_pred_c=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                       smodel.clim.comp.gam,clampt=T),
                                      ba_pred_c=1-s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],
                                                      smodel.clim.comp.gam,clampt=T),
                                      ppt_pred_c=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],
                                                       smodel.clim.comp.gam,clampt=T),
                                      t_pred_c=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],
                                                     smodel.clim.comp.gam,clampt=T)),seq)
splot_data_int<-cbind(data.frame(dia_pred=1-s_fun_mean(seq$dia,means[2],means[3],means[4],means[5],
                                                smodel.int.gam),
                                 ba_pred=1-s_fun_mean(means[1],seq$ba,means[3],means[4],means[5],smodel.int.gam),
                                 ppt_pred=1-s_fun_mean(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.gam),
                                 t_pred=1-s_fun_mean(means[1],means[2],means[3],seq$t,means[5],smodel.int.gam),
                                 ppt_pred_25=1-s_fun_mean(means[1],means[2],seq$ppt,twentyfive[4],means[5],smodel.int.gam),
                                 t_pred_25=1-s_fun_mean(means[1],means[2],twentyfive[3],seq$t,means[5],smodel.int.gam),
                                 ppt_pred_75=1-s_fun_mean(means[1],means[2],seq$ppt,seventyfive[4],means[5],smodel.int.gam),
                                 t_pred_75=1-s_fun_mean(means[1],means[2],seventyfive[3],seq$t,means[5],smodel.int.gam)),seq)

#Survival with fire
survData$PLT_CN_factor<-as.factor(survData$PLT_CN)

s_fun_fire<-function(dia,ba,ppt,t,ci,plot,model,clampt=F){
  sdata=data.frame(PREVDIA=dia,BALIVE=ba,PPT_yr_norm = ppt, 
                   T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling.fire$scale[match(names(sdata), 
                                                                        names(surv.scaling.fire$scale))], 
                                  center = surv.scaling.fire$center[match(names(sdata), 
                                                                          names(surv.scaling.fire$center))])) 
  scaled.sdata = cbind(scaled.sdata,CENSUS_INTERVAL = ci)
  scaled.sdata$PLT_CN_factor = plot
  return(predict(model, newdata = scaled.sdata, type = "response", re.form = NA))
}

s_fun_mean_fire<-function(dia,ba,ppt,t,ci,model,clampt=F){
  sdata=data.frame(PREVDIA=dia,BALIVE=ba,PPT_yr_norm = ppt, 
                   T_yr_norm = ifelse(clampt==T & t >12.5, 12.5, t))
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling.fire$scale[match(names(sdata), 
                                                                   names(surv.scaling.fire$scale))], 
                                  center = surv.scaling.fire$center[match(names(sdata), 
                                                                     names(surv.scaling.fire$center))])) 
  scaled.sdata = cbind(scaled.sdata,CENSUS_INTERVAL = ci)
  scaled.sdata$PLT_CN_factor = 14546600020004
  return(predict(model, newdata = scaled.sdata, type = "response", re.form = NA, exclude = "s(PLT_CN_factor)"))
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
survData$resid_c_fire<-survData$mort-s_fun_fire(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                      survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                      survData$PLT_CN_factor,smodel.clim.fire.gam)
survData$resid_ci_fire<-survData$mort-s_fun_fire(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                       survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                       survData$PLT_CN_factor,smodel.clim.int.fire.gam)
survData$resid_cc_fire<-survData$mort-s_fun_fire(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                       survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                       survData$PLT_CN_factor,smodel.clim.comp.fire.gam)
survData$resid_i_fire<-survData$mort-s_fun_fire(survData$PREVDIA,survData$BALIVE,survData$PPT_yr_norm,
                                      survData$T_yr_norm,survData$CENSUS_INTERVAL,
                                      survData$PLT_CN_factor,smodel.int.fire.gam)

ncuts=30
chopsize_dia<-cut(survData$PREVDIA,ncuts)
chopsize_ba<-cut(survData$BALIVE,ncuts)
chopsize_PPT<-cut(survData$PPT_yr_norm,ncuts)
chopsize_T<-cut(survData$T_yr_norm,ncuts)

dia_binned<-as.vector(sapply(split(survData$PREVDIA,chopsize_dia),mean,na.rm=T))
count_binned_dia<-as.vector(sapply(split(survData$mort,chopsize_dia),length))
surv_binned_dia_c<-1-(as.vector(sapply(split(survData$resid_c_fire,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.fire.gam)))
surv_binned_dia_ci<-1-(as.vector(sapply(split(survData$resid_ci_fire,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.int.fire.gam)))
surv_binned_dia_cc<-1-(as.vector(sapply(split(survData$resid_cc_fire,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(dia_binned,means[2],means[3],means[4],means[5],smodel.clim.comp.fire.gam)))
surv_binned_dia_i<-1-(as.vector(sapply(split(survData$resid_i_fire,chopsize_dia),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(dia_binned,means[2],means[3],means[4],means[5],smodel.int.fire.gam)))

ba_binned<-as.vector(sapply(split(survData$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(survData$mort,chopsize_ba),length))
surv_binned_ba_c<-1-(as.vector(sapply(split(survData$resid_c_fire,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.fire.gam)))
surv_binned_ba_ci<-1-(as.vector(sapply(split(survData$resid_ci_fire,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.int.fire.gam)))
surv_binned_ba_cc<-1-(as.vector(sapply(split(survData$resid_cc_fire,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],ba_binned,means[3],means[4],means[5],smodel.clim.comp.fire.gam)))
surv_binned_ba_i<-1-(as.vector(sapply(split(survData$resid_i_fire,chopsize_ba),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],ba_binned,means[3],means[4],means[5],smodel.int.fire.gam)))

PPT_binned<-as.vector(sapply(split(survData$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(survData$mort,chopsize_PPT),length))
surv_binned_PPT_c<-1-(as.vector(sapply(split(survData$resid_c_fire,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.fire.gam)))
surv_binned_PPT_ci<-1-(as.vector(sapply(split(survData$resid_ci_fire,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.int.fire.gam)))
surv_binned_PPT_cc<-1-(as.vector(sapply(split(survData$resid_cc_fire,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],PPT_binned,means[4],means[5],smodel.clim.comp.fire.gam)))
surv_binned_PPT_i<-1-(as.vector(sapply(split(survData$resid_i_fire,chopsize_PPT),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],PPT_binned,means[4],means[5],smodel.int.fire.gam)))

T_binned<-as.vector(sapply(split(survData$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(survData$mort,chopsize_T),length))
surv_binned_T_c<-1-(as.vector(sapply(split(survData$resid_c_fire,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],means[3],T_binned,means[5],smodel.clim.fire.gam)))
surv_binned_T_ci<-1-(as.vector(sapply(split(survData$resid_ci_fire,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],means[3],T_binned,means[5],smodel.clim.int.fire.gam)))
surv_binned_T_cc<-1-(as.vector(sapply(split(survData$resid_cc_fire,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],means[3],T_binned,means[5],smodel.clim.comp.fire.gam)))
surv_binned_T_i<-1-(as.vector(sapply(split(survData$resid_i_fire,chopsize_T),mean,na.rm=T))+
  as.vector(s_fun_mean_fire(means[1],means[2],means[3],T_binned,means[5],smodel.int.fire.gam)))

s_binned_fire<-as.data.frame(cbind(surv_binned_dia_c,surv_binned_ba_c,surv_binned_PPT_c,surv_binned_T_c,
                              surv_binned_dia_ci,surv_binned_ba_ci,surv_binned_PPT_ci,surv_binned_T_ci,
                              surv_binned_dia_cc,surv_binned_ba_cc,surv_binned_PPT_cc,surv_binned_T_cc,
                              surv_binned_dia_i,surv_binned_ba_i,surv_binned_PPT_i,surv_binned_T_i,
                              dia_binned,ba_binned,PPT_binned,T_binned,
                              count_binned_dia,count_binned_ba,count_binned_PPT,count_binned_T))

names(s_binned_fire)<-c("mort_dia_c","mort_ba_c","mort_PPT_c","mort_T_c",
                   "mort_dia_ci","mort_ba_ci","mort_PPT_ci","mort_T_ci",
                   "mort_dia_cc","mort_ba_cc","mort_PPT_cc","mort_T_cc",
                   "mort_dia_i","mort_ba_i","mort_PPT_i","mort_T_i",
                   "PREVDIA","BALIVE","PPT","T",
                   "count_dia","count_ba","count_PPT","count_T")

splot_data_clim_fire<-cbind(data.frame(dia_pred=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                      smodel.clim.fire.gam),
                                  ppt_pred=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim.fire.gam),
                                  t_pred=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],smodel.clim.fire.gam),
                                  dia_pred_c=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                        smodel.clim.fire.gam,clampt=T),
                                  ppt_pred_c=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],
                                                        smodel.clim.fire.gam,clampt=T),
                                  t_pred_c=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],smodel.clim.fire.gam,
                                                      clampt=T)),seq)
splot_data_climint_fire<-cbind(data.frame(dia_pred=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                         smodel.clim.int.fire.gam),
                                     ppt_pred=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],smodel.clim.int.fire.gam),
                                     t_pred=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],smodel.clim.int.fire.gam),
                                     dia_pred_c=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                           smodel.clim.int.fire.gam,clampt=T),
                                     ppt_pred_c=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],
                                                           smodel.clim.int.fire.gam,clampt=T),
                                     t_pred_c=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],smodel.clim.int.fire.gam,
                                                         clampt=T)),seq)
splot_data_climcomp_fire<-cbind(data.frame(dia_pred=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                          smodel.clim.comp.fire.gam),
                                      ba_pred=1-s_fun_mean_fire(means[1],seq$ba,means[3],means[4],means[5],
                                                         smodel.clim.comp.fire.gam),
                                      ppt_pred=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],
                                                          smodel.clim.comp.fire.gam),
                                      t_pred=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],
                                                        smodel.clim.comp.fire.gam),
                                      dia_pred_c=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                            smodel.clim.comp.fire.gam,clampt=T),
                                      ba_pred_c=1-s_fun_mean_fire(means[1],seq$ba,means[3],means[4],means[5],
                                                           smodel.clim.comp.fire.gam,clampt=T),
                                      ppt_pred_c=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],
                                                            smodel.clim.comp.fire.gam,clampt=T),
                                      t_pred_c=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],
                                                          smodel.clim.comp.fire.gam,clampt=T)),seq)
splot_data_int_fire<-cbind(data.frame(dia_pred=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],
                                                     smodel.int.fire.gam),
                                 ba_pred=1-s_fun_mean_fire(means[1],seq$ba,means[3],means[4],means[5],smodel.int.fire.gam),
                                 ppt_pred=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.fire.gam),
                                 t_pred=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],smodel.int.fire.gam),
                                 dia_pred_c=1-s_fun_mean_fire(seq$dia,means[2],means[3],means[4],means[5],smodel.int.fire.gam,
                                                       clampt=T),
                                 ba_pred_c=1-s_fun_mean_fire(means[1],seq$ba,means[3],means[4],means[5],smodel.int.fire.gam,
                                                      clampt=T),
                                 ppt_pred_c=1-s_fun_mean_fire(means[1],means[2],seq$ppt,means[4],means[5],smodel.int.fire.gam,
                                                       clampt=T),
                                 t_pred_c=1-s_fun_mean_fire(means[1],means[2],means[3],seq$t,means[5],smodel.int.fire.gam,
                                                     clampt=T)),seq)

#Recruitment
r_fun<-function(ba,ppt,t,pa,ci,model,clampba=F){
  rdata=data.frame(BALIVE = ifelse(clampba==T & ba >204, 204, 
                                   ifelse(clampba==T & ba <93, 93, ba)),
                   PPT_yr_norm = ppt, T_yr_norm = t)
  #PPT_yr_norm = ppt, T_yr_norm = t)
  scaled.rdata = data.frame(scale(rdata, 
                                  scale = r.scaling$scale[match(names(rdata), 
                                                                names(r.scaling$scale))], 
                                  center = r.scaling$center[match(names(rdata), 
                                                                  names(r.scaling$center))])) 
  scaled.rdata$CENSUS_INTERVAL = ci
  scaled.rdata$PIEDadults1 = pa
  scaled.rdata$off = log(scaled.rdata$CENSUS_INTERVAL) + log(scaled.rdata$PIEDadults1)
  return((predict(model, newdata = scaled.rdata, type = "response")))
}

means<-c(mean(rdata$BALIVE,na.rm=T),
         mean(rdata$PPT_yr_norm,na.rm=T),mean(rdata$T_yr_norm,na.rm=T),
         mean(rdata$PIEDadults1),mean(rdata$CENSUS_INTERVAL))
twentyfive<-c(quantile(rdata$BALIVE,na.rm=T)[2],
              quantile(rdata$PPT_yr_norm,na.rm=T)[2],quantile(rdata$T_yr_norm,na.rm=T)[2])
seventyfive<-c(quantile(rdata$BALIVE,na.rm=T)[4],
               quantile(rdata$PPT_yr_norm,na.rm=T)[4],quantile(rdata$T_yr_norm,na.rm=T)[4])

seq<-data.frame(ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(rdata$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))

#Calculate residuals
rdata$resid_c<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                     rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.clim.gam)
rdata$resid_ci<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                     rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.clim.int.gam)
rdata$resid_cc<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                      rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.clim.comp.gam)
rdata$resid_i<-rdata$recruits1-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                     rdata$PIEDadults1,rdata$CENSUS_INTERVAL,rmodel.int.gam)

ncuts=30
chopsize_ba<-cut(rdata$BALIVE,ncuts)
chopsize_PPT<-cut(rdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(rdata$T_yr_norm,ncuts)

ba_binned<-as.vector(sapply(split(rdata$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(rdata$recruits1,chopsize_ba),length))
recr_binned_ba_c<-as.vector(sapply(split(rdata$resid_c,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.clim.gam))
recr_binned_ba_ci<-as.vector(sapply(split(rdata$resid_ci,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.clim.int.gam))
recr_binned_ba_cc<-as.vector(sapply(split(rdata$resid_cc,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.clim.comp.gam))
recr_binned_ba_i<-as.vector(sapply(split(rdata$resid_i,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],means[5],rmodel.int.gam))

PPT_binned<-as.vector(sapply(split(rdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(rdata$recruits1,chopsize_PPT),length))
recr_binned_PPT_c<-as.vector(sapply(split(rdata$resid_c,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.clim.gam))
recr_binned_PPT_ci<-as.vector(sapply(split(rdata$resid_ci,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.clim.int.gam))
recr_binned_PPT_cc<-as.vector(sapply(split(rdata$resid_cc,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.clim.comp.gam))
recr_binned_PPT_i<-as.vector(sapply(split(rdata$resid_i,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],means[5],rmodel.int.gam))

T_binned<-as.vector(sapply(split(rdata$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(rdata$recruits1,chopsize_T),length))
recr_binned_T_c<-as.vector(sapply(split(rdata$resid_c,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.clim.gam))
recr_binned_T_ci<-as.vector(sapply(split(rdata$resid_ci,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.clim.int.gam))
recr_binned_T_cc<-as.vector(sapply(split(rdata$resid_cc,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.clim.comp.gam))
recr_binned_T_i<-as.vector(sapply(split(rdata$resid_i,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],means[5],rmodel.int.gam))

r_binned<-as.data.frame(cbind(recr_binned_ba_c,recr_binned_PPT_c,recr_binned_T_c,
                              recr_binned_ba_ci,recr_binned_PPT_ci,recr_binned_T_ci,
                              recr_binned_ba_cc,recr_binned_PPT_cc,recr_binned_T_cc,
                              recr_binned_ba_i,recr_binned_PPT_i,recr_binned_T_i,
                              ba_binned,PPT_binned,T_binned,
                              count_binned_ba,count_binned_PPT,count_binned_T))
names(r_binned)<-c("recr_ba_c","recr_PPT_c","recr_T_c",
                   "recr_ba_ci","recr_PPT_ci","recr_T_ci",
                   "recr_ba_cc","recr_PPT_cc","recr_T_cc",
                   "recr_ba_i","recr_PPT_i","recr_T_i",
                   "BALIVE","PPT","T","count_ba","count_PPT","count_T")

rplot_data_clim<-cbind(data.frame(ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                 rmodel.clim.gam),
                                  t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim.gam),
                                  ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                   rmodel.clim.gam,clampba=T),
                                  t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim.gam,
                                                 clampba=T)),seq)
rplot_data_climint<-cbind(data.frame(ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                 rmodel.clim.int.gam),
                                  t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim.int.gam),
                                  ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                   rmodel.clim.int.gam,clampba=T),
                                  t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.clim.int.gam,
                                                 clampba=T)),seq)
rplot_data_climcomp<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                    rmodel.clim.comp.gam),
                                      ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                     rmodel.clim.comp.gam),
                                      t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],
                                                   rmodel.clim.comp.gam),
                                      ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                                      rmodel.clim.comp.gam,clampba=T),
                                      ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],means[5],
                                                       rmodel.clim.comp.gam,clampba=T),
                                      t_pred_c=r_fun(means[1],means[2],seq$t,means[4],means[5],
                                                     rmodel.clim.comp.gam,clampba=T)),seq)
rplot_data_int<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],means[5],
                                               rmodel.int.gam),
                                 ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],means[5],rmodel.int.gam),
                                 t_pred=r_fun(means[1],means[2],seq$t,means[4],means[5],rmodel.int.gam),
                                 ppt_pred_25=r_fun(means[1],seq$ppt,twentyfive[3],means[4],means[5],rmodel.int.gam),
                                 t_pred_25=r_fun(means[1],twentyfive[2],seq$t,means[4],means[5],rmodel.int.gam),
                                 ppt_pred_75=r_fun(means[1],seq$ppt,seventyfive[3],means[4],means[5],rmodel.int.gam),
                                 t_pred_75=r_fun(means[1],seventyfive[2],seq$t,means[4],means[5],rmodel.int.gam)),seq)


save(grdata,survData,rdata,g_binned,s_binned,r_binned,
     grplot_data_clim,grplot_data_climint,grplot_data_climcomp,grplot_data_int, 
     splot_data_clim,splot_data_climint,splot_data_climcomp,splot_data_int, 
     rplot_data_clim,rplot_data_climint,rplot_data_climcomp,rplot_data_int, 
     file="./Output/vital_effects_gam.rda")
