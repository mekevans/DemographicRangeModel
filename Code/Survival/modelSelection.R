library(coefplot)
library(ggplot2)
library(ggeffects)
library(cowplot)
library(dplyr)
library(effects)

# Read in data
survData <- read.csv("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

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

# Recode status
survData$surv <- ifelse(survData$STATUSCD == 2, 0, 1)
survData$mort <- ifelse(survData$STATUSCD == 1, 0, 1)

# transformed predictors
survData$log.PPT_yr <- log(survData$PPT_yr)
survData$log.PPT_yr_norm <- log(survData$PPT_yr_norm)

survData$log.size <- log(survData$PREVDIA)
survData$log.BALIVE <- log(survData$BALIVE)

# survData <- subset(survData, log.BALIVE > 0) 
survData.2 <- subset(survData, BALIVE > 0) # goes from 20329 to 20161

# standardize covariates
survData.scaled <- survData %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                          -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                          -CENSUS_INTERVAL, 
                                                          -AGB_INCR, -DIA_INCR, -BA_INCR,
                                                          -surv, -mort))

survData2.scaled <- survData.2 %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                              -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                              -CENSUS_INTERVAL, 
                                                              -AGB_INCR, -DIA_INCR, -BA_INCR,
                                                              -surv, -mort))

library(lme4)
# model with PREVDIA instead of BAt1
# AIC = 17318.9
smodel1 <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + 
                   PPT_yr_norm + PPT_dr_anom + T_yr_norm + T_dr_anom + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel1))

res = simulateResiduals(smodel1)

plot(res)

plotResiduals(survData.scaled$PREVDIA, res$scaledResiduals, quantreg = T, main = "PREVDIA") #resid plot doesn't look good
plotResiduals(survData.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE") #resid plots a bit off here as well
plotResiduals(survData.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm") # better, but not exactly on the dashed lines
plotResiduals(survData.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm") # looks good

# not updated: remove T_dr_anom and delta AIC = 5 (17387.6)
smodel1.a <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + 
                   PPT_yr_norm + PPT_dr_anom + T_yr_norm + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# with quadratic terms, AIC = 17268.0
smodel1.q <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + 
                     T_dr_anom + I(T_dr_anom^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# remove PPT anomaly; AIC = 17265.8
smodel2.q <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     T_yr_norm + I(T_yr_norm^2) + 
                     T_dr_anom + I(T_dr_anom^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# remove T anomaly; AIC = 17265.1
smodel3.q <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     T_yr_norm + I(T_yr_norm^2) +  
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# remove BALIVE; AIC = 17282.9
smodel4.q <- glmer(mort ~ PREVDIA + I(PREVDIA^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     T_yr_norm + I(T_yr_norm^2) +  
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# interactions, AIC = 17335.7
smodel6 <- glmer(mort ~ (PREVDIA + BALIVE + PPT_yr_norm + T_yr_norm)^2 + 
                   I(PREVDIA^2) + I(BALIVE^2) + I(T_yr_norm^2) + I(PPT_yr_norm^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# this is after filtering the BALIVE = 0 data points
# needs to be updated - rerun with PREVDIA as size
# AIC = 17096.5
smodel12 <- glmer(mort ~ (log.size + log.BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                    I(log.size^2) + I(log.BALIVE^2) + I(T_yr_norm^2) +
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData2.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel12))

res = simulateResiduals(smodel12)

plot(res)

plotResiduals(survData2.scaled$log.size, res$scaledResiduals, quantreg = T, main = "log.size")
plotResiduals(survData2.scaled$log.BALIVE, res$scaledResiduals, quantreg = T, main = "log.BALIVE")
plotResiduals(survData2.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm")
plotResiduals(survData2.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "v")


# AIC = 17169.1 
smodel1.alt <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + 
                     T_dr_anom + I(T_dr_anom^2) + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData2.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel1.alt))
res = simulateResiduals(smodel1.alt)
plot(res)
plotResiduals(survData2.scaled$PREVDIA, res$scaledResiduals, quantreg = T, main = "PREVDIA") # not good
plotResiduals(survData2.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE") # a bit wonky (lines not horizontal)
plotResiduals(survData2.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm")
plotResiduals(survData2.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm")

# AIC = 17167.4
smodel2 <- glmer(mort ~ PREVDIA + I(PREVDIA^2) + 
                       PPT_yr_norm + I(PPT_yr_norm^2) +
                       T_yr_norm + I(T_yr_norm^2) + 
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData2.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel2))

# AIC = 17164.6
smodel3 <- glmer(mort ~ (PREVDIA + PPT_yr_norm + T_yr_norm)^2 + 
                   I(PREVDIA^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData2.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
res = simulateResiduals(smodel3) # deviation of residuals is significant (doesn't <look> worse than others)
plot(res) # deviation not significant
plotResiduals(survData2.scaled$PREVDIA, res$scaledResiduals, quantreg = T, main = "PREVDIA") # does NOT look good
plotResiduals(survData2.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE") # almost OK
plotResiduals(survData2.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm") # pretty good
plotResiduals(survData2.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm")


# use log-transformed size and BALIVE
# AIC = 17098.8
smodel1.log <- glmer(mort ~ log.size + I(log.size^2) + log.BALIVE + I(log.BALIVE^2) +
                       PPT_yr_norm + I(PPT_yr_norm^2) +
                       PPT_dr_anom + I(PPT_dr_anom^2) + 
                       T_yr_norm + I(T_yr_norm^2) + 
                       T_dr_anom + I(T_dr_anom^2) + 
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData2.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel1.log))
res = simulateResiduals(smodel1.log) # deviation of residuals is significant (doesn't <look> worse than others)
plot(res)
plotResiduals(survData2.scaled$log.size, res$scaledResiduals, quantreg = T, main = "log.size") # looks better to me (than not log transformed)
plotResiduals(survData2.scaled$log.BALIVE, res$scaledResiduals, quantreg = T, main = "log.BALIVE") # not good
plotResiduals(survData2.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm")
plotResiduals(survData2.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "v")

### try model with PREVDIA as predictor of mortality (tree size)
# AIC = 17276.7
smodel4 <- glmer(mort ~ (PREVDIA + PPT_yr_norm + T_yr_norm)^2 + 
                   I(PREVDIA^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
res = simulateResiduals(smodel4) # deviation of residuals is significant (doesn't <look> worse than others)
plot(res) # deviation not significant, but p=0.07221 (nearly significant)
plotResiduals(survData2.scaled$PREVDIA, res$scaledResiduals, quantreg = T, main = "PREVDIA") # not great, but maybe better that BAt1?
plotResiduals(survData2.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE") # almost OK
plotResiduals(survData2.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm") # pretty good
plotResiduals(survData2.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm")



####################################################################################
# without log transform (two filters)
# AIC = 17139.4
smodel12.3nolog <- glmer(mort ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                      I(BAt1^2) + I(BALIVE^2) + I(T_yr_norm^2) +
                      (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                    family = binomial(link = cloglog), data = survData3.scaled,
                    control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel12.3nolog))

res = simulateResiduals(smodel12.3nolog)

plot(res)

plotResiduals(survData3.scaled$BAt1, res$scaledResiduals, quantreg = T, main = "size")
plotResiduals(survData3.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE")
plotResiduals(survData3.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm")
plotResiduals(survData3.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "v")

# remove precip anomaly (which is NS); AIC = 17139.6
smodel13.3nolog <- glmer(mort ~ (BAt1 + BALIVE + PPT_yr_norm + T_yr_norm)^2 + 
                           I(BAt1^2) + I(BALIVE^2) + I(T_yr_norm^2) +
                           (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                         family = binomial(link = cloglog), data = survData3.scaled,
                         control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# remove BALIVE (which is NS); AIC = 17140.6 
smodel14.3nolog <- glmer(mort ~ (BAt1 + PPT_yr_norm + T_yr_norm)^2 + 
                           I(BAt1^2) + I(T_yr_norm^2) +
                           (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                         family = binomial(link = cloglog), data = survData3.scaled,
                         control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# add quadratic on PPT_yr_norm; AIC = 17127.3
smodel15.3nolog <- glmer(mort ~ (BAt1 + PPT_yr_norm + T_yr_norm)^2 + 
                           I(BAt1^2) + I(T_yr_norm^2) + I(PPT_yr_norm^2) + 
                           (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                         family = binomial(link = cloglog), data = survData3.scaled,
                         control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel15.3nolog))
res = simulateResiduals(smodel12.3nolog)

plot(res)

plotResiduals(survData3.scaled$BAt1, res$scaledResiduals, quantreg = T, main = "size")
plotResiduals(survData3.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE")
plotResiduals(survData3.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm")
plotResiduals(survData3.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "v")

# add in PPT quadratic with (ns) BALIVE - just for visualization of effects
# AIC = 17126.6
smodel16.3nolog <- glmer(mort ~ (BAt1 + BALIVE + PPT_yr_norm + T_yr_norm)^2 + 
                           I(BAt1^2) + I(BALIVE^2) + I(T_yr_norm^2) + I(PPT_yr_norm^2) + 
                           (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                         family = binomial(link = cloglog), data = survData3.scaled,
                         control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel16.3nolog))


# competing models that have only quadratic terms, no interactions
# AIC = 17132.9
smodel1.3alt <- glmer(mort ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                       PPT_yr_norm + I(PPT_yr_norm^2) +
                       PPT_dr_anom + I(PPT_dr_anom^2) + 
                       T_yr_norm + I(T_yr_norm^2) + 
                       T_dr_anom + I(T_dr_anom^2) + 
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel1.3alt))

smodel2.3alt <- glmer(mort ~ BAt1 + I(BAt1^2) +
                        PPT_yr_norm + I(PPT_yr_norm^2) +
                        T_yr_norm + I(T_yr_norm^2) + 
                        (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                      family = binomial(link = cloglog), data = survData3.scaled,
                      control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel2.3alt))

# use log-transformed size and BALIVE
# AIC = 17068.8
smodel1.3log <- glmer(mort ~ log.size + I(log.size^2) + log.BALIVE + I(log.BALIVE^2) +
                       PPT_yr_norm + I(PPT_yr_norm^2) +
                       PPT_dr_anom + I(PPT_dr_anom^2) + 
                       T_yr_norm + I(T_yr_norm^2) + 
                       T_dr_anom + I(T_dr_anom^2) + 
                       (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                     family = binomial(link = cloglog), data = survData3.scaled,
                     control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel1.3log))

res = simulateResiduals(smodel1.3log)

plot(res)

plotResiduals(survData3.scaled$BAt1, res$scaledResiduals, quantreg = T, main = "BAt1")
plotResiduals(survData3.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE")
plotResiduals(survData3.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm")
plotResiduals(survData3.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "v")








# AIC = 17844.7
smodel.q2 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BAt1:I(BALIVE^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17843.6
smodel.q3 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BAt1:T_yr_norm +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17847.9
smodel.q4 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BAt1:I(T_yr_norm^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17847.2
smodel.q5 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BAt1:PPT_dr_anom +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17847.7
smodel.q6 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BAt1:I(PPT_dr_anom^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17847.4
smodel.q7 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BALIVE:PPT_dr_anom +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17847.9
smodel.q8 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + I(BALIVE^2) +
                     PPT_yr_norm + I(PPT_yr_norm^2) +
                     PPT_dr_anom + I(PPT_dr_anom^2) + 
                     T_yr_norm + I(T_yr_norm^2) + BALIVE:I(PPT_dr_anom^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))


# narrow down which season drought anomaly is important
# cool season? AIC = 17897
smodel1.b <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                   PPT_yr_norm + PPT_c_dr_anom + T_yr_norm + T_dr_anom + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# monsoon? AIC = 17897
smodel1.c <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_m_dr_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# foresummer? AIC = 17892
smodel1.d <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_fs_dr_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# previous fall? AIC = 17885
smodel1.e <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_pf_dr_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# cool season? AIC = 17901
smodel1.f <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_c_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# monsoon? AIC = 17899
smodel1.g <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_m_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# foresummer? AIC = 17901
smodel1.h <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_fs_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# BEST SO FAR previous fall AIC = 17873
smodel1.i <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_pf_anom + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# cool season? AIC = 17897
smodel1.j <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_c_dr + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# monsoon? AIC = 17891
smodel1.k <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_m_dr + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# foresummer? AIC = 17898
smodel1.l <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_fs_dr + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# previous fall? AIC = 17892
smodel1.m <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                     PPT_yr_norm + PPT_pf_dr + T_yr_norm + T_dr_anom + 
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17884 no improvement
smodel2 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                   PPT_c_norm + PPT_wd_norm + PPT_m_norm + PPT_dr_anom + T_yr_norm + T_dr_anom + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# doesn't like the following model, failed to converge, "nearly unidentifiable"
smodel3 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                   PPT_c + PPT_wd + PPT_m + PPT_dr_anom + T_yr_norm + T_dr_anom + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# AIC = 17895
smodel4 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                   PPT_drought + T_yr_norm + T_dr_anom + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# doesn't like this model, failed to converge, "nearly unidentifiable"
smodel5 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE + 
                   PPT_drought + T_yr_norm + 
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

# AIC = 17918
smodel20 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE +
                PPT_c + PPT_wd + PPT_m + VPD_c + VPD_wd + VPD_m +
                  (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                family = binomial(link = cloglog), data = survData.scaled,
                control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
smodel21 <- glmer(surv ~ BAt1 + I(BAt1^2) + BALIVE +
                    PPT_c + PPT_wd + PPT_m + T_c + T_wd + T_m +
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))


# AIC = 17863 (got rid of T_dr_anom)
smodel7 <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                   I(BAt1^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))


# add quadratic terms into (nearly) best model with 2-way interactions
# this won't run...failed to converge and model is nearly unidentifiable
smodel7.b <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                   I(BAt1^2) + I(BALIVE^2) + I(PPT_yr_norm^2) + I(PPT_dr_anom^2) + I(T_yr_norm^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# not worse, and effect of BALIVE^2 is significant (AIC = 17844.2)
smodel7.c <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                     I(BAt1^2) + I(BALIVE^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# 10 AIC units better! 17834.2
# AIC = 17800.3 after removing largest BAt1 values
# model elicits warnings after removing smallest BALIVE values
smodel7.d <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                     I(BAt1^2) + I(BALIVE^2) + I(T_yr_norm^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel7.d), cex = 0.5)
# AIC = 17753.6 (after removing BAt1 outliers)
# AIC = 17549.9 (after also removing 168 cases where log.BALIVE < 0)
smodel11 <- glmer(surv ~ (log.size + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                     I(log.size^2) + I(BALIVE^2) + I(T_yr_norm^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17539.9
smodel12 <- glmer(surv ~ (log.size + log.BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                    I(log.size^2) + I(log.BALIVE^2) + I(T_yr_norm^2) +
                    (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                  family = binomial(link = cloglog), data = survData.scaled,
                  control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
plot(allEffects(smodel12), cex = 0.5)
#pdf
plot(effect("PPT_yr_norm", smodel12))

#####################################################
## THIS IS THE BEST MODEL AS OF 7/26/2018 (smodel12)

# slightly worse, AIC = 17835.5; PPT_dr_anom^2 is not significant
smodel7.e <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                     I(BAt1^2) + I(BALIVE^2) + I(T_yr_norm^2) + I(PPT_dr_anom^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# two warnings, not sure this model can be trusted...
# AIC is 3 units lower, and it says PPT_yr_norm^2 is significant
# maybe worth investigating with a GAM?
smodel7.f <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_dr_anom + T_yr_norm)^2 + 
                     I(BAt1^2) + I(BALIVE^2) + I(T_yr_norm^2) + I(PPT_yr_norm^2) +
                     (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                   family = binomial(link = cloglog), data = survData.scaled,
                   control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))



# doesn't like this model, failed to converge, "nearly unidentifiable"
smodel8 <- glmer(surv ~ (BAt1 + BALIVE + PPT_yr_norm + PPT_pf_anom + T_yr_norm)^2 + 
                   I(BAt1^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))


# support for 3 seasons (PPT)? AIC = 17854
smodel9 <- glmer(surv ~ (BAt1 + BALIVE + 
                          PPT_c_norm + PPT_wd_norm + PPT_m_norm + PPT_dr_anom + 
                  T_yr_norm + T_dr_anom)^2 + 
                  I(BAt1^2) +
                  (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                family = binomial(link = cloglog), data = survData.scaled,
                control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))
# AIC = 17854.7
smodel10 <- glmer(surv ~ (BAt1 + BALIVE + 
                           PPT_c_norm + PPT_wd_norm + PPT_m_norm + PPT_dr_anom + 
                           T_yr_norm)^2 + 
                   I(BAt1^2) +
                   (1|PLT_CN) + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = survData.scaled,
                 control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=10000)))

plot(allEffects(smodel1))


#### NEW STUFF 10/25/18
### dealing with std'ized covariates

# specify the predictors in the "best" model (or candidate best)
surv.predictors <- c("PREVDIA", "T_yr_norm", "PPT_yr_norm") # "I(PREVDIA)", "I(T_yr_norm)", "I(PPT_yr_norm)"
# eventually rewrite this so that it can handle alternative "best" models

get_scale = function(data, predictors) {
  sc = list("scale" = NULL, "center"  = NULL)
  for (i in predictors) {
    sc$scale[i] = attributes(data[, i])$"scaled:scale"
    sc$center[i] = attributes(data[, i])$"scaled:center"
  }
  return(sc)
}

surv.scaling = get_scale(survData.scaled, surv.predictors)

# remove scaling information from the dataset so that the model doesnt expect scaled data in predict()
for (i in surv.predictors) {
  attributes(survData.scaled[, i]) = NULL
}

# export model for coefficients and scaling information -------------------
save(smodel4.q, surv.scaling, file = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Code/IPM/SurvRescaling.Rdata")



### OLD STUFF - Michiel


# Model: size + elevation + BA + climate
smodel <- update(smodel, . ~ . - PPT_c)
smodel <- update(smodel, . ~ . - PPT_w)
smodel <- update(smodel, . ~ . - VPD_c)
smodel <- update(smodel, . ~ . + I(PREV_DRYBIO_AG^2))
smodel <- update(smodel, . ~ . + I(PREV_DRYBIO_AG^2) + I(PREV_DRYBIO_AG^3))
# smodel <- update(smodel, . ~ . + I(VPD_w^2) + I(VPD_w^3))
smodel <- update(smodel, . ~ . - VPD_w)
smodel <- update(smodel, . ~ . + log(VPD_w))
summary(smodel)
save(smodel, file = "D:/EvansLab/Final/Models/BC/surv.rda")
pdf("D:/EvansLab/Final/Output/BC/SurvivalModel.pdf")
coefplot(smodel, intercept = F)
preds <- ggpredict(smodel, terms = c("PREV_DRYBIO_AG"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. biomass")
preds <- ggpredict(smodel, terms = c("elev"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. elevation")
preds <- ggpredict(smodel, terms = c("baLive"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. plot BA")
preds <- ggpredict(smodel, terms = c("VPD_w"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. VPD warm")
dev.off()
# Publication figure
preds <- ggpredict(smodel, terms = c("PREV_DRYBIO_AG"))
A <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Biomass")
preds <- ggpredict(smodel, terms = c("elev"))
B <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Elevation")
preds <- ggpredict(smodel, terms = c("baLive"))
C <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Plot basal area")
preds <- ggpredict(smodel, terms = c("VPD_w"))
D <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Warm-season VPD")
all <- plot_grid(A, B, C, D, labels = c("A", "B", "C", "D"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/FigS6.png", all, base_aspect_ratio = 1.5)

# Model: size + elevation + BA
smodel <- glm(dead ~ PREV_DRYBIO_AG + I(PREV_DRYBIO_AG^2) + I(PREV_DRYBIO_AG^3) + elev + baLive, family = "binomial", data = survData)
summary(smodel)
save(smodel, file = "D:/EvansLab/Final/Models/B/surv.rda")
pdf("D:/EvansLab/Final/Output/B/SurvivalModel.pdf")
coefplot(smodel, intercept = F)
preds <- ggpredict(smodel, terms = c("PREV_DRYBIO_AG"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. biomass")
preds <- ggpredict(smodel, terms = c("elev"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. elevation")
preds <- ggpredict(smodel, terms = c("baLive"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. plot BA")
dev.off()
# Publication figure
preds <- ggpredict(smodel, terms = c("PREV_DRYBIO_AG"))
A <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Biomass")
preds <- ggpredict(smodel, terms = c("elev"))
B <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Elevation")
preds <- ggpredict(smodel, terms = c("baLive"))
C <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Plot basal area")
all <- plot_grid(A, B, C, labels = c("A", "B", "C"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/FigS4.png", all, base_aspect_ratio = 1.5)

# Model: size + elevation + climate
smodel <- glm(dead ~ PREV_DRYBIO_AG + elev + PPT_c + PPT_w + VPD_c + VPD_w, family = "binomial", data = survData)
smodel <- update(smodel, . ~ . - PPT_c)
smodel <- update(smodel, . ~ . - PPT_w)
smodel <- update(smodel, . ~ . - VPD_c)
smodel <- update(smodel, . ~ . + I(PREV_DRYBIO_AG^2) + I(PREV_DRYBIO_AG^3))
summary(smodel)
save(smodel, file = "D:/EvansLab/Final/Models/C/surv.rda")
pdf("D:/EvansLab/Final/Output/C/SurvivalModel.pdf")
coefplot(smodel, intercept = F)
preds <- ggpredict(smodel, terms = c("PREV_DRYBIO_AG"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. biomass")
preds <- ggpredict(smodel, terms = c("elev"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. elevation")
preds <- ggpredict(smodel, terms = c("VPD_w"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Mortality vs. VPD warm")
dev.off()
# Publication figure
preds <- ggpredict(smodel, terms = c("PREV_DRYBIO_AG"))
A <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Biomass")
preds <- ggpredict(smodel, terms = c("elev"))
B <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Elevation")
preds <- ggpredict(smodel, terms = c("VPD_w"))
C <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Mortality") + xlab("Warm-season VPD")
all <- plot_grid(A, B, C, labels = c("A", "B", "C"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/FigS5.png", all, base_aspect_ratio = 1.5)

