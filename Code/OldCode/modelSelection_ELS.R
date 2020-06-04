#### Filter out fire and harvest for all models
#### Don't log transform
####### Emily's PIED survival model exploration:
library(MuMIn)
library(coefplot)
library(ggplot2)
library(ggeffects)
library(cowplot)
library(dplyr)
library(effects)
library(DHARMa)
library(lme4)

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

##Start with this basic model:
sbase<-glmer(mort ~ PREVDIA + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19593.2
sbase_balive<-glmer(mort ~ PREVDIA + BALIVE + 
               (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
             family = binomial(link = cloglog), data = survData3.scaled,
             control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19592.9
sbase_ppt<-glmer(mort ~ PREVDIA + PPT_yr_norm + 
               (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
             family = binomial(link = cloglog), data = survData3.scaled,
             control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
sbase_t<-glmer(mort ~ PREVDIA + T_yr_norm + 
               (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
             family = binomial(link = cloglog), data = survData3.scaled,
             control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
sbase_clim<-glmer(mort ~ PREVDIA + PPT_yr_norm + T_yr_norm + 
               (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
             family = binomial(link = cloglog), data = survData3.scaled,
             control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
sbase_balive_clim<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + 
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

mod.comp1<-model.sel(sbase,sbase_balive,sbase_ppt,sbase_t,sbase_clim,sbase_balive_clim)

#Add interactions

#Add quadratics

#Add interactions and quadratics

## Compare diameter vs basal area as size predictor
sbase_ba<-glmer(mort ~ BAt1 + BALIVE + PPT_yr_norm + T_yr_norm + 
                  (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                family = binomial(link = cloglog), data = survData.scaled,
                control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 20047.3
mod.comp_size<-model.sel(sbase_dia,sbase_ba) #Diameter better
res = simulateResiduals(sbase_dia)
plot(res)
res = simulateResiduals(sbase_ba)
plot(res)
#resid plots look similar

## Compare size and log(size)
sbase_logd<-glmer(mort ~ log.size + BALIVE + PPT_yr_norm + T_yr_norm + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19909.0
sbase_logba<-glmer(mort ~ log.size + log.BALIVE + PPT_yr_norm + T_yr_norm + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData2.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000))) #Gives the invalid grouping factor error
mod.comp_log<-model.sel(sbase_dia,sbase_log) #log better
res = simulateResiduals(sbase_log)
plot(res) #looks similar

## Compare all data vs ba>0 vs no fire/harvest
sbase_ba0<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + 
                              (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                            family = binomial(link = cloglog), data = survData2.scaled,
                            control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19756.9
sbase_noFH<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + 
                              (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                            family = binomial(link = cloglog), data = survData3.scaled,
                            control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19592.9
mod.comp_data<-model.sel(sbase_dia,sbase_ba0,sbase_noFH) #no fire/harvest has lower AIC, but different data, so comparison doesn't make much sense
res = simulateResiduals(sbase_noFH)
plot(res)

## Compare climate norm to census intervals
sbase_Pcen<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr + T_yr_norm + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19593.7
sbase_Tcen<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr + 
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19588.3
sbase_cen<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr + T_yr + 
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19587.7
mod.comp_cen<-model.sel(sbase_dia,sbase_Pcen,sbase_Tcen,sbase_cen) #cen and Tcen were best

## Compare annual and seasonal climate
sbase_c<-glmer(mort ~ PREVDIA + BALIVE + PPT_c + T_c + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19568.1
sbase_wd<-glmer(mort ~ PREVDIA + BALIVE + PPT_wd + T_wd + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19580.7
sbase_pf<-glmer(mort ~ PREVDIA + BALIVE + PPT_pf + T_pf + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19610.8
sbase_fs<-glmer(mort ~ PREVDIA + BALIVE + PPT_fs + T_fs + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19562.3
sbase_m<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_m + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19585.0
mod.comp_clim<-model.sel(sbase_cen,sbase_c,sbase_wd,sbase_pf,sbase_fs,sbase_m) #fs best, seemed to be driven by temp

sbase_c_ppt<-glmer(mort ~ PREVDIA + BALIVE + PPT_c + T_fs + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19560.0
sbase_wd_ppt<-glmer(mort ~ PREVDIA + BALIVE + PPT_wd + T_fs + 
                  (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                family = binomial(link = cloglog), data = survData3.scaled,
                control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19565.0
sbase_pf_ppt<-glmer(mort ~ PREVDIA + BALIVE + PPT_pf + T_fs + 
                  (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                family = binomial(link = cloglog), data = survData3.scaled,
                control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19564.7
sbase_m_ppt<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19538.5
mod.comp_clim2<-model.sel(sbase_fs,sbase_c_ppt,sbase_wd_ppt,sbase_pf_ppt,sbase_m_ppt) #m was best

## Add anomalies
sbase_Panom<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                     PPT_dr_anom +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19538.7
sbase_Tanom<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                     T_dr_anom +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19535.0
sbase_anom<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                     PPT_dr_anom + T_dr_anom +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19533.8
mod.comp_anom<-model.sel(sbase_m_ppt,sbase_Panom,sbase_Tanom,sbase_anom) #both anom and Tanom were best

## Which season drought anomaly? - precip
#season_dr
sbase_anom_Pc<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                    PPT_c_dr_anom + T_dr_anom +
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19533.5
sbase_anom_Pm<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_m_dr_anom + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19529.9
sbase_anom_Pfs<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_fs_dr_anom + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19523.6
sbase_anom_Ppf<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_pf_dr_anom + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19519.7
mod.comp_sanom1<-model.sel(sbase_anom,sbase_anom_Pc,sbase_anom_Pm,sbase_anom_Pfs,sbase_anom_Ppf) #Ppf was best
#season
sbase_anom_Pc2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_c_anom + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19537.0
sbase_anom_Pm2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_m_anom + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19536.9
sbase_anom_Pfs2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                        PPT_fs_anom + T_dr_anom +
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19537.0
sbase_anom_Ppf2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                        PPT_pf_anom + T_dr_anom +
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19533.9
mod.comp_sanom2<-model.sel(sbase_anom_Ppf,sbase_anom_Pc2,sbase_anom_Pm2,sbase_anom_Pfs2,sbase_anom_Ppf2) #Ppf still best
#no anom
sbase_anom_Pc3<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_c_dr + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19532.7
sbase_anom_Pm3<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                       PPT_m_dr + T_dr_anom +
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19530.7
sbase_anom_Pfs3<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                        PPT_fs_dr + T_dr_anom +
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19525.3
sbase_anom_Ppf3<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                        PPT_pf_dr + T_dr_anom +
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19531.7
mod.comp_sanom3<-model.sel(sbase_anom_Ppf,sbase_anom_Pc3,sbase_anom_Pm3,sbase_anom_Pfs3,sbase_anom_Ppf3) #Ppf still best

## Which season drought anomaly? - temp
#season
sbase_anom_Tc2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                        PPT_pf_dr_anom + T_c_anom +
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19516.7
sbase_anom_Tm2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                        PPT_pf_dr_anom + T_m_anom +
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19528.4
sbase_anom_Tfs2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                         PPT_pf_dr_anom + T_fs_anom +
                         (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                       family = binomial(link = cloglog), data = survData3.scaled,
                       control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19525.6
sbase_anom_Tpf2<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + 
                         PPT_pf_dr_anom + T_pf_anom +
                         (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                       family = binomial(link = cloglog), data = survData3.scaled,
                       control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19526.3
mod.comp_sanom4<-model.sel(sbase_anom_Ppf,sbase_anom_Tc2,sbase_anom_Tm2,sbase_anom_Tfs2,sbase_anom_Tpf2) #Tc2 was best

## Add quadratics
sbase_q1a<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(PREVDIA^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19466.6
sbase_q1b<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(BALIVE^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19516.2
sbase_q1c<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(PPT_yr_norm^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19500.2
sbase_q1d<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_dr_anom + 
                   I(T_yr_norm^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19489.2
sbase_q1e<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(PPT_pf_dr_anom^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19518.2
sbase_q1f<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(T_c_anom^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19518.7
mod.comp_quad1<-model.sel(sbase_anom_Tc2,sbase_q1a,sbase_q1b,sbase_q1c,sbase_q1d,sbase_q1e,sbase_q1f) #q1a best, all quadratic models except f and marginally e better than no quadratic
sbase_q_ef<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                 I(PREVDIA^2) +I(BALIVE^2) +I(PPT_yr_norm^2) +I(T_yr_norm^2) + 
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19427.0
sbase_q_f<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(PREVDIA^2) +I(BALIVE^2) +I(PPT_yr_norm^2) +I(T_yr_norm^2) + I(PPT_pf_dr_anom^2) + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19428.7
sbase_q<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                 I(PREVDIA^2) +I(BALIVE^2) +I(PPT_yr_norm^2) +I(T_yr_norm^2) + I(PPT_pf_dr_anom^2) + I(T_c_anom^2) +
                 (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
               family = binomial(link = cloglog), data = survData3.scaled,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19429.6
mod.comp_quad2<-model.sel(sbase_q1a,sbase_q_ef,sbase_q_f,sbase_q) #_ef and _f were best

## Add interactions
sbase_int<-glmer(mort ~ (PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom)^2 + 
                   I(PREVDIA^2) +I(BALIVE^2) +I(PPT_m_norm^2) +I(T_fs_norm^2) + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19393.1
mod.comp_int<-model.sel(sbase_q_ef,sbase_int) #Interactions better
# Only significant interactions:
sbase_int_sig<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                   I(PREVDIA^2) +I(BALIVE^2) +I(PPT_m_norm^2) +I(T_fs_norm^2) + I(PPT_pf_dr_anom^2) + 
                     PREVDIA:BALIVE + PREVDIA:PPT_m + PREVDIA:T_c_anom + T_fs:PPT_pf_dr_anom + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19379.5
# + next "most" significant
sbase_int_mar<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                       I(PREVDIA^2) +I(BALIVE^2) +I(PPT_m_norm^2) +I(T_fs_norm^2) + I(PPT_pf_dr_anom^2) + 
                       PREVDIA:BALIVE + PREVDIA:PPT_m + PREVDIA:T_c_anom + BALIVE:T_c_anom + T_fs:PPT_pf_dr_anom + 
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19377.5
mod.comp_int2<-model.sel(sbase_int,sbase_int_sig,sbase_int_mar) #mar was best

########Best model so far###########
sbase_int_mar<-glmer(mort ~ PREVDIA + BALIVE + PPT_m + T_fs + PPT_pf_dr_anom + T_c_anom + 
                       I(PREVDIA^2) +I(BALIVE^2) +I(PPT_m_norm^2) +I(T_fs_norm^2) + I(PPT_pf_dr_anom^2) + 
                       PREVDIA:BALIVE + PREVDIA:PPT_m + PREVDIA:T_c_anom + BALIVE:T_c_anom + T_fs:PPT_pf_dr_anom + 
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#AIC = 19712.3
res = simulateResiduals(sbase_int_sig)
plot(res)