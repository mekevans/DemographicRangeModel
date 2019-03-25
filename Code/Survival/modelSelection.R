library(coefplot)
library(ggplot2)
library(ggeffects)
library(cowplot)
library(dplyr)
library(effects)
library(DHARMa)

# Read in data
survData <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

# Create increment columns
# not needed for survival/mort analysis
survData$AGB_INCR <- survData$DRYBIO_AG_DIFF / survData$CENSUS_INTERVAL
survData$DIA_INCR <- survData$DIA_DIFF / survData$CENSUS_INTERVAL
survData$BA_INCR <- survData$BA_DIFF / survData$CENSUS_INTERVAL

# distribution of size...update to evaluate PREVDIA
hist(survData$PREVDIA) # size is lognormal-ish, with a wonky bit at the small end of the scale (due to min size threshold?)
hist(survData$PREVDIA, breaks = c(seq(0, 1220, by = 10)), xlim = c(0, 500))
hist(survData$BAt1, breaks = c(seq(0, 1220, by = 10)), xlim = c(0, 1200))
hist(log(survData$BAt1)) # not quite lognormal, but worth trying in the models below (heavy in the left tail)

# distribution of other predictors
hist(survData$PPT_yr) # not too bad...Poisson-ish but with a large mean count
hist(log(survData$PPT_yr)) # more normal-looking
hist(survData$T_yr) # looks normal
hist(survData$PPT_yr_norm) # not too bad...Poisson-ish but with a large mean count
hist(log(survData$PPT_yr_norm)) # more normal-looking
hist(survData$T_yr_norm) # looks normal
hist(survData$BALIVE) # not too bad...Poisson-ish but with a large mean count
hist(log(survData$BALIVE)) # log transform has a heavy left tail

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


library(lme4)
#library(lmerTest)
library(MuMIn) # use MuMin to choose between models (AICc)

#Compare different size predictors
smodel.0a <- glmer(mort ~ BAt1 + I(BAt1^2) + BALIVE + PPT_yr + T_yr + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.0b <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + T_yr + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.0c <- glmer(mort ~ PREV_DRYBIO_AG + I(PREV_DRYBIO_AG^2) + BALIVE + PPT_yr + T_yr + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

plot(allEffects(smodel.0a))
plot(allEffects(smodel.0b))
plot(allEffects(smodel.0c))

res = simulateResiduals(smodel.0a)
plot(res)
#p=0.0837
res = simulateResiduals(smodel.0b)
plot(res)
#p=0.05726
res = simulateResiduals(smodel.0c)
plot(res)
#p=0.04531

plotResiduals(survData.scaled$PREVDIA, res$scaledResiduals, quantreg = T, main = "PREVDIA") #resid plot doesn't look good
plotResiduals(survData.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE") #resid plots a bit off here as well
plotResiduals(survData.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm") # better, but not exactly on the dashed lines
plotResiduals(survData.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm") # looks good

#Demonstrate the effect of quadratics
smodel.lin <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                    family = binomial(link = cloglog), data = survData3.scaled,
                    control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.q <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + 
                    I(PREVDIA^2) + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
mod.comp0<-model.sel(smodel.lin,smodel.q)
# strong preference for quadratics (delta AIC = 93.38)

# compare VPD and temperature
smodel.1a <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + VPD_yr_norm + 
                    I(PREVDIA^2) + I(BALIVE^2) + I(PPT_yr_norm^2) + I(VPD_yr_norm^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData3.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.1b <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + 
                     I(PREVDIA^2) + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
mod.comp1<-model.sel(smodel.1a,smodel.1b)
# strong prefence for temperature (temperature also preferred in growth model)

# compare annual vs. 3 vs. 4 seasons...likes 4-season best
smodel.2a <- glmer(mort ~ PREVDIA + BALIVE + PPT_c + PPT_wd + PPT_m + T_c + T_wd + T_m + 
                     I(PREVDIA^2) + I(BALIVE^2) + I(PPT_c^2) + I(PPT_wd^2) + I(PPT_m^2) + 
                     I(T_c^2) + I(T_wd^2) + (T_m^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.2b <- glmer(mort ~ PREVDIA + BALIVE + PPT_c + PPT_pf + PPT_fs + PPT_m + T_c + T_pf + T_fs + T_m + 
                     I(PREVDIA^2) + I(BALIVE^2) + I(PPT_c^2) + I(PPT_pf^2) + I(PPT_fs^2) + I(PPT_m^2) + 
                     I(T_c^2) + I(T_pf^2) + I(T_fs^2) + (T_m^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

mod.comp2 <- model.sel(smodel.1b, smodel.2a, smodel.2b)

# compare normals vs. census interval vs. anomalies...
smodel.3a <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_yr_norm^2) + I(T_yr_norm^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3b <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3c <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_anom + T_yr_anom + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_yr_anom^2) + I(T_yr_anom^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3d <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_anom + PPT_pf_anom + PPT_fs_anom + PPT_m_anom + 
                     T_c_anom + T_pf_anom + T_fs_anom + T_m_anom + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_c_anom^2) + I(PPT_pf_anom^2) + I(PPT_fs_anom^2) + I(PPT_m_anom^2) + 
                     I(T_c_anom^2) + I(T_pf_anom^2) + I(T_fs_anom^2) + (T_m_anom^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3e <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_anom + T_yr_anom + PPT_yr_norm + T_yr_norm + 
                     I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_yr_anom^2) + I(T_yr_anom^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3f <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_anom + PPT_pf_anom + PPT_fs_anom + PPT_m_anom + 
                     T_c_anom + T_pf_anom + T_fs_anom + T_m_anom + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_c_anom^2) + I(PPT_pf_anom^2) + I(PPT_fs_anom^2) + I(PPT_m_anom^2) + 
                     I(T_c_anom^2) + I(T_pf_anom^2) + I(T_fs_anom^2) + (T_m_anom^2) + 
                     PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + 
                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3g <- glmer(mort ~ PREVDIA + BALIVE + PPTex_yr_anom + Tex_yr_anom + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPTex_yr_anom^2) + I(Tex_yr_anom^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3h <- glmer(mort ~ PREVDIA + BALIVE + PPTex_c_anom + PPTex_pf_anom + PPTex_fs_anom + PPTex_m_anom + 
                     Tex_c_anom + Tex_pf_anom + Tex_fs_anom + Tex_m_anom + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPTex_c_anom^2) + I(PPTex_pf_anom^2) + I(PPTex_fs_anom^2) + I(PPTex_m_anom^2) + 
                     I(Tex_c_anom^2) + I(Tex_pf_anom^2) + I(Tex_fs_anom^2) + (Tex_m_anom^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3i <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_yr_norm^2) + I(T_yr_norm^2) + PPTex_yr_anom + Tex_yr_anom + 
                     I(PPTex_yr_anom^2) + I(Tex_yr_anom^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.3j <- glmer(mort ~ PREVDIA + BALIVE + PPTex_c_anom + PPTex_pf_anom + PPTex_fs_anom + PPTex_m_anom + 
                     Tex_c_anom + Tex_pf_anom + Tex_fs_anom + Tex_m_anom + I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPTex_c_anom^2) + I(PPTex_pf_anom^2) + I(PPTex_fs_anom^2) + I(PPTex_m_anom^2) + 
                     I(Tex_c_anom^2) + I(Tex_pf_anom^2) + I(Tex_fs_anom^2) + (Tex_m_anom^2) + 
                     PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + 
                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

mod.comp3 <- model.sel(smodel.2b, smodel.3a, smodel.3b, smodel.3c, smodel.3d, smodel.3e, smodel.3f, 
                       smodel.3g, smodel.3h, smodel.3i, smodel.3j)
# smodel.3j and smodel.3f were best among these, but gave warning message (In commonArgs(par, fn, control, environment()) :
# maxfun < 10 * length(par)^2 is not recommended.), probably because model is so complex, so proceeding with 
# the next best, smodel.3b (delta AIC = 9.56)

# add drought anomalies
# There were convergence issues using 4-season climate data plus drought, 
# so instead I am comparing annual plus drought to 4-season without drought,
# using annual norms and anomalies, which was the best annual model from mod.comp3

#smodel.4a <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
#                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + PPT_drought + Tmean_drought + 
#                     I(PREVDIA^2) + I(BALIVE^2) + 
#                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
#                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) + 
#                     I(PPT_drought^2) + I(Tmean_drought^2) + 
#                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
#                   family = binomial(link = cloglog), data = survData3.scaled,
#                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#smodel.4b <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
#                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + 
#                     PPT_pf_dr + PPT_c_dr + PPT_fs_dr + PPT_m_dr + Tmean_drought + 
#                     I(PREVDIA^2) + I(BALIVE^2) + 
#                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
#                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) + 
#                     I(PPT_c_dr^2) + I(PPT_fs_dr^2) + I(PPT_m_dr^2) + I(Tmean_drought^2) + 
#                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
#                   family = binomial(link = cloglog), data = survData3.scaled,
#                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#smodel.4c <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_anom + PPT_pf_anom + PPT_fs_anom + PPT_m_anom + 
#                     T_c_anom + T_pf_anom + T_fs_anom + T_m_anom + PPT_dr_anom + T_dr_anom + 
#                     I(PREVDIA^2) + I(BALIVE^2) + I(T_c_anom^2) + I(T_pf_anom^2) + I(T_fs_anom^2) + 
#                     (T_m_anom^2) + PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
#                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + I(PREVDIA^2) + I(BALIVE^2) + 
#                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
#                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) +
#                     I(PPT_dr_anom^2) + I(T_dr_anom^2) + 
#                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
#                   family = binomial(link = cloglog), data = survData3.scaled,
#                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
#smodel.4d <- glmer(mort ~ PREVDIA + BALIVE + PPT_c_anom + PPT_pf_anom + PPT_fs_anom + PPT_m_anom + 
#                     T_c_anom + T_pf_anom + T_fs_anom + T_m_anom + 
#                     PPT_pf_dr_anom + PPT_c_dr_anom + PPT_fs_dr_anom + PPT_m_dr_anom + T_dr_anom + 
#                     I(PREVDIA^2) + I(BALIVE^2) + I(T_c_anom^2) + I(T_pf_anom^2) + I(T_fs_anom^2) + 
#                     (T_m_anom^2) + PPT_c_norm + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
#                     T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + I(PREVDIA^2) + I(BALIVE^2) + 
#                     I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
#                     I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) +
#                     I(PPT_c_dr_anom^2) + I(PPT_fs_dr_anom^2) + I(PPT_m_dr_anom^2) + I(T_dr_anom^2) + 
#                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
#                   family = binomial(link = cloglog), data = survData3.scaled,
#                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

smodel.4a <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + PPT_yr_anom + T_yr_anom + 
                     I(PREVDIA^2) + I(BALIVE^2) + PPT_drought + Tmean_drought + 
                     I(PPT_yr_norm^2) + I(T_yr_norm^2) + I(PPT_yr_anom^2) + I(T_yr_anom^2) + 
                     I(PPT_drought^2) + I(Tmean_drought^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.4b <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + PPT_yr_anom + T_yr_anom + 
                     I(PREVDIA^2) + I(BALIVE^2) + 
                     PPT_pf_dr + PPT_c_dr + PPT_fs_dr + PPT_m_dr + Tmean_drought + 
                     I(PPT_yr_norm^2) + I(T_yr_norm^2) + I(PPT_yr_anom^2) + I(T_yr_anom^2) + 
                     I(PPT_pf_dr^2) + I(PPT_c_dr^2) + I(PPT_fs_dr^2) + I(PPT_m_dr^2) + I(Tmean_drought^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.4c <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + PPT_yr_anom + T_yr_anom + 
                     I(PREVDIA^2) + I(BALIVE^2) + PPT_dr_anom + T_dr_anom + 
                     I(PPT_yr_norm^2) + I(T_yr_norm^2) + I(PPT_yr_anom^2) + I(T_yr_anom^2) + 
                     I(PPT_dr_anom^2) + I(T_dr_anom^2) + (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.4d <- glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + PPT_yr_anom + T_yr_anom + 
                     I(PREVDIA^2) + I(BALIVE^2) + 
                     PPT_pf_dr_anom + PPT_c_dr_anom + PPT_fs_dr_anom + PPT_m_dr_anom + T_dr_anom + 
                     I(PPT_yr_norm^2) + I(T_yr_norm^2) + I(PPT_yr_anom^2) + I(T_yr_anom^2) + 
                     I(PPT_pf_dr_anom^2) + I(PPT_c_dr_anom^2) + I(PPT_fs_dr_anom^2) + I(PPT_m_dr_anom^2) + I(T_dr_anom^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

mod.comp4<-model.sel(smodel.3b,smodel.4a,smodel.4b,smodel.4c,smodel.4d)
# smodel.3b is still best

# add 2-way interactions, excluding quadratics
# There were convergence issues using 4-season climate data plus drought, 
# so instead I am comparing annual plus interactions to 4-season without interactions,
# using annual norms and anomalies, which was the best annual model from mod.comp3

smodel.5a <- glmer(mort ~ (PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + PPT_yr_anom + T_yr_anom)^2 + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel.5b <- glmer(mort ~ (PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + PPT_yr_anom + T_yr_anom)^2 + 
                      I(PREVDIA^2) + I(BALIVE^2) + 
                     I(PPT_yr_norm^2) + I(PPT_yr_anom^2) + I(T_yr_norm^2) + I(T_yr_anom^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

mod.comp5 <- model.sel(smodel.3b, smodel.5a, smodel.5b)
# smodel.3b is still best


#### NEW STUFF 10/25/18
### dealing with std'ized covariates

# Models to export:
smodel.clim<-glmer(mort ~ PREVDIA + PPT_yr_norm + T_yr_norm + 
                           I(PREVDIA^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) +
                           (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                         family = binomial(link = cloglog), data = survData3.scaled,
                         control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

smodel.clim.comp<-glmer(mort ~ PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm + 
                          I(PREVDIA^2) + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
                         (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData3.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

smodel.int<-glmer(mort ~ (PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm)^2 + 
                     I(PREVDIA^2) + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData3.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

smodel.best<-smodel.3b

# specify the predictors in the "best" model (or candidate best)
#surv.predictors <- c("PREVDIA", "T_yr_norm", "PPT_yr_norm")
surv.predictors <- c("PREVDIA", "T_yr_norm", "PPT_yr_norm", "BALIVE", "PPT_c_norm", "PPT_pf_norm", 
                     "PPT_fs_norm", "PPT_m_norm", "T_c_norm", "T_pf_norm", "T_fs_norm", "T_m_norm")
# eventually rewrite this so that it can handle alternative "best" models

get_scale = function(data, predictors) {
  sc = list("scale" = NULL, "center"  = NULL)
  for (i in predictors) {
    sc$scale[i] = attributes(data[, i])$"scaled:scale"
    sc$center[i] = attributes(data[, i])$"scaled:center"
  }
  return(sc)
}

surv.scaling = get_scale(survData3.scaled, surv.predictors)

# remove scaling information from the dataset so that the model doesnt expect scaled data in predict()
for (i in surv.predictors) {
  attributes(survData.scaled[, i]) = NULL
}

# export model for coefficients and scaling information -------------------
#save(smodel4.q, surv.scaling, file = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Code/IPM/SurvRescaling.Rdata")
#save(smodel4, surv.scaling, file = "./Code/IPM/SurvRescalingNoFire.Rdata")
#save(smodel3, surv.scaling, file = "./Code/IPM/SurvRescalingBA.Rdata")
save(smodel.clim,smodel.clim.comp,smodel.int,smodel.best, surv.scaling, file = "./Code/IPM/SurvRescalingNoFire.Rdata")
