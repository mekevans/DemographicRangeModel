library(ggeffects)
library(ggplot2)
library(coefplot)
library(effects)

# Read and process data
rdata <- read.csv("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)

#rdata <- read.csv("C:/Users/mekevans/Documents/Cdrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PIED_IPM/MEKEvans/Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)

# looking at the data
recruit.table <- table(rdata[, c("PIEDadults1", "recruits1")])
colSums(recruit.table)
table(rdata[, c("PIEDadults4", "recruits1")])

plot.info <- c("plot", "lat", "lon", "elev", "state", "county", "plotID", "PIEDadults1", "recruits1")
colonize.plots <- rdata[rdata$PIEDadults1 == 0 & rdata$recruits1 > 0, plot.info]

# Subset data: only plots that contained PIED trees at time 1 census
# AGB: note that it was previously incorrectly time 2
rdata <- subset(rdata, PIEDadults1 > 0)

# look at distribution of response
hist(rdata$recruits1)
hist(rdata$recruits1, breaks = c(seq(0, 8, by = 1)), xlim = c(0, 8))
hist(rdata$PIEDadults1, breaks = c(seq(0,53, by=1)), xlim = c(0, 53))
plot(rdata$PIEDadults1, rdata$recruits1)
plot(rdata$BALIVE, rdata$recruits1)
plot(rdata$PPT_yr_norm, rdata$recruits1)

# look at distribution of (some) covariates
hist(rdata$BALIVE) # not too bad...Poisson-ish but with a large mean count
hist(rdata$BALIVE, breaks = c(seq(0, 750, by = 10)), xlim = c(0, 500))
# looks like there's one outlier with very high BALIVE, could be eliminated
rdata2 <- subset(rdata, BALIVE < 400) # 6720
hist(rdata2$BALIVE, breaks = c(seq(0, 400, by = 10)), xlim = c(0, 400))
hist(log(rdata$BALIVE)) # log transform has a heavy left tail

# compare BA.all to BALIVE...they should be the same
# both in sq feet, but BALIVE is on a per acre basis
# BA.all is simply the sum per plot (plot = 672.5 m^2 or 0.16617837 acre)
# expansion factor is 6.01763*BA.all
rdata$BA.all <- rdata$BA.PIED + rdata$BA.notPIED
plot(rdata$BALIVE, 6.01763*rdata$BA.all) # mostly fall on/near one:one, but BALIVE < BA.all

rdata3 <- subset(rdata, BALIVE > 0) # eliminate data points with no live trees at time 1
# 6635 plots
rdata$log.BALIVE <- log(rdata$BALIVE)

# precip variables look better after log transformation
# temp/VPD variables don't need transformation
hist(rdata$PPT_c_window_15)
hist(log(rdata$PPT_c_window_15)) # log transform is an improvement
hist(rdata$T_c_window_15) # doesn't need transformation

# standardize covariates
library(dplyr)
rdata.scaled <- rdata %>% mutate_at(scale, .vars = vars(-plot, -lat, -lon, -elev, -PApied,
                                                        -state, -county, -plotID, -CONDID, 
                                                          -measYear, -plotPrev, -PREV_MEASYEAR,
                                                          -CENSUS_INTERVAL, -recruits1, -recruits12,
                                                          -AGB_intra, -BA.PIED, -PIEDadults1,
                                                          -PIEDadults4, -PIEDadults8, -cumDIA.PIED))


# AIC = 5663.4
rmodel <- glm(data = rdata.scaled, recruits1 ~ PPT_yr_norm + T_yr_norm + BALIVE +
                offset(log(CENSUS_INTERVAL)) + offset(log(PIEDadults1 + 0.0001)), 
              family = "poisson")
### Lisa_simple
# AIC = 5722.4 with VPD, AIC = 5716 with T (poisson)
# doesn't run with family = "nb", don't know why
# with normals instead of 25-yr window, AIC = 5720.9
# with 20-yr window AIC = 5716.6 (indistinguishable from 25-yr window)
# with 15-yr window AIC = 5712.6 (some support for 15-yr window)
rmodel_full <- glm(data = rdata.scaled, recruits1 ~ 1
                   + offset(log(PIEDadults1 + 0.0001))
                   + BALIVE + I(BALIVE^2)
                   + PPT_c_window_15 + I(PPT_c_window_15^2)
                   + PPT_m_window_15 + I(PPT_m_window_15^2)
                   + T_c_window_15 + I(T_c_window_15^2)
                   + T_m_window_15 + I(T_m_window_15^2),
                   family = "poisson", offset = log(CENSUS_INTERVAL))

# MAT and MAP model is worse, AIC = 5621
rmodel_simple <- glm(data = rdata2.scaled, recruits1 ~ 1
                   + offset(log(PIEDadults1 + 0.0001))
                   + BALIVE + I(BALIVE^2)
                   + PPT_yr_norm + I(PPT_yr_norm^2)
                   + T_yr_norm + I(T_yr_norm^2),
                   family = "poisson", offset = log(CENSUS_INTERVAL))
plot(allEffects(rmodel_simple))

# try a logistic regression of recruitment recoded to 0/1
rdata.scaled$recruit01 <- ifelse(rdata.scaled$recruits1 > 0, 1, 0)
# model output doesn't make a lot of sense
# BALIVE ns (negative, but ns)
# climate effects more significant than ZIP, weird looking (unimodal)
rmodel_logreg <- glm(recruit01 ~ BALIVE + I(BALIVE^2)
                 + PPT_yr_norm + I(PPT_yr_norm^2)
                 + T_yr_norm + I(T_yr_norm^2)
                 + offset(log(CENSUS_INTERVAL)), 
                 family = binomial(link = cloglog), data = rdata.scaled)
summary(rmodel_logreg)
plot(allEffects(rmodel_logreg))

res = simulateResiduals(smodel1)

plot(res)

plotResiduals(survData.scaled$BAt1, res$scaledResiduals, quantreg = T, main = "BAt1") #resid plot doesn't look good
plotResiduals(survData.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "BALIVE") #resid plots a bit off here as well
plotResiduals(survData.scaled$PPT_yr_norm, res$scaledResiduals, quantreg = T, main = "PPT_yr_norm") # better, but not exactly on the dashed lines
plotResiduals(survData.scaled$T_yr_norm, res$scaledResiduals, quantreg = T, main = "T_yr_norm") # looks good




### from LISA also
# install.packages("R2admb")
# install.packages("glmmADMB", repos=c("http://glmmadmb.r-forge.r-project.org/repos", getOption("repos")), type="source")
library(R2admb)
library(glmmADMB)

# AIC = 4358.3
rmodel_zip <- glmmadmb(recruits1 ~ 1
                             + BALIVE + I(BALIVE^2)
                             + PPT_yr_norm + I(PPT_yr_norm^2)
                             + T_yr_norm + I(T_yr_norm^2)
                             + offset(log(CENSUS_INTERVAL))
                             + offset(log(PIEDadults1)), 
                             zeroInflation = TRUE,
                             data = rdata.scaled, 
                             family = "poisson")
plot(allEffects(rmodel_zip))


library(DHARMa)
library(glmmTMB)

# http://glmmadmb.r-forge.r-project.org/glmmADMB.html

glmmTMB()

# first glmmTMB model
# climate normals in Poisson
rmodel_zip <- glmmTMB(recruits1 ~ 1
                      + BALIVE + I(BALIVE^2)
                      + PPT_yr_norm + I(PPT_yr_norm^2)
                      + T_yr_norm + I(T_yr_norm^2)
                      + offset(log(CENSUS_INTERVAL))
                      + offset(log(PIEDadults1)), # various alternatives for this offset
                      ziformula = ~ 1 
                      + PPT_yr_norm + I(PPT_yr_norm^2)  
                      + T_yr_norm + I(T_yr_norm^2)
                      + BALIVE + I(BALIVE^2),
                      data = rdata.scaled, 
                      family = "poisson")
summary(rmodel_zip) # AIC = 4346.1
res = simulateResiduals(rmodel_zip)
plot(res, quantreg = T) # ns: p = 0.09159
# must read in Effect.glmmTMB function
# see https://github.com/glmmTMB/glmmTMB/blob/7ba86a972ddb13226a8de9eab0e113e6156fccf4/glmmTMB/R/effects.R
plot(Effect.glmmTMB("BALIVE", rmodel_zip))
plot(Effect.glmmTMB("PPT_yr_norm", rmodel_zip))

# with offset = PIEDadults4, can't be run without adding arbitrary non-zero value, since ~400 plots have zero PIED trees >4" DRC
# with offset = PIEDadults8, same problem as PIEDadults4
# with offset = cumDIA.PIED, AIC = 4573.6
# with offset = BA.PIED, AIC = 4893.0
rmodel_zipcumDIA <- glmmTMB(recruits1 ~ 1
                             + BALIVE + I(BALIVE^2)
                             + PPT_yr_norm + I(PPT_yr_norm^2)
                             + T_yr_norm + I(T_yr_norm^2)
                             + offset(log(CENSUS_INTERVAL))
                             + offset(log(cumDIA.PIED)), # various alternatives for this offset
                             ziformula = ~ 1 
                             + PPT_yr_norm + I(PPT_yr_norm^2)  
                             + T_yr_norm + I(T_yr_norm^2)
                             + BALIVE + I(BALIVE^2),
                             data = rdata.scaled, 
                             family = "poisson")
summary(rmodel_zipcumDIA) # AIC = 4573.6
res = simulateResiduals(rmodel_zipcumDIA)
plot(res, quantreg = T) # significant deviation: p = 0.02349

rmodel_zip_PIEDBA <- glmmTMB(recruits1 ~ 1
                            + BALIVE + I(BALIVE^2)
                            + PPT_yr_norm + I(PPT_yr_norm^2)
                            + T_yr_norm + I(T_yr_norm^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(BA.PIED)), # various alternatives for this offset
                            ziformula = ~ 1 
                            + PPT_yr_norm + I(PPT_yr_norm^2)  
                            + T_yr_norm + I(T_yr_norm^2)
                            + BALIVE + I(BALIVE^2),
                            data = rdata.scaled, 
                            family = "poisson")
summary(rmodel_zip_PIEDBA) # AIC = 4893.0
res = simulateResiduals(rmodel_zip_PIEDBA)
plot(res, quantreg = T) # significant deviation: p = 0.00022

# instead of climate normals
# best guesses as to time interval relevant for actual recruitment
rmodel_zip25 <- glmmTMB(recruits1 ~ 1
                           + BALIVE + I(BALIVE^2)
                           + PPT_yr_window_25 + I(PPT_yr_window_25^2)
                           + T_yr_window_25 + I(T_yr_window_25^2)
                           + offset(log(CENSUS_INTERVAL))
                           + offset(log(PIEDadults1)), 
                           ziformula = ~ 1 
                           + PPT_yr_norm + I(PPT_yr_norm^2)  
                           + T_yr_norm + I(T_yr_norm^2)
                           + BALIVE + I(BALIVE^2),
                           data = rdata.scaled, 
                           family = "poisson")
summary(rmodel_zip25) # AIC = 4343.9
res = simulateResiduals(rmodel_zip25)
plot(res, quantreg = T) # deviation of the residuals is NOT significant (p = 0.08891)


# time frame for Poisson = 20 yrs
rmodel_zip20 <- glmmTMB(recruits1 ~ 1
                          + BALIVE + I(BALIVE^2)
                          + PPT_yr_window_20 + I(PPT_yr_window_20^2)
                          + T_yr_window_20 + I(T_yr_window_20^2)
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)), 
                          ziformula = ~ 1 
                          + PPT_yr_norm + I(PPT_yr_norm^2)  
                          + T_yr_norm + I(T_yr_norm^2)
                          + BALIVE + I(BALIVE^2),
                          data = rdata.scaled, 
                          family = "poisson")
summary(rmodel_zip20) # AIC = 4342.1
res = simulateResiduals(rmodel_zip20)
plot(res, quantreg = T) # deviation is NS: p = 0.09852

# 15 yrs (Poisson)
rmodel_zip15 <- glmmTMB(recruits1 ~ 1
                        + BALIVE + I(BALIVE^2)
                        + PPT_yr_window_15 + I(PPT_yr_window_15^2)
                        + T_yr_window_15 + I(T_yr_window_15^2)
                        + offset(log(CENSUS_INTERVAL))
                        + offset(log(PIEDadults1)), 
                        ziformula = ~ 1 
                        + PPT_yr_norm + I(PPT_yr_norm^2)  
                        + T_yr_norm + I(T_yr_norm^2)
                        + BALIVE + I(BALIVE^2),
                        data = rdata.scaled, 
                        family = "poisson")
summary(rmodel_zip15) # AIC = 4344.8
res = simulateResiduals(rmodel_zip15)
plot(res, quantreg = T) # deviation is NS: p = 0.05795


# use PIPO basal area, AIC = 5615.8
rmodel_zipoiss <- glmmTMB(recruits1 ~ 1
                          + BA.PIPO + I(BA.PIPO^2)
                          + PPT_yr_window_20 + I(PPT_yr_window_20^2)
                          + T_yr_window_20 + I(T_yr_window_20^2)
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)), 
                          ziformula = ~ 1 
                          + PPT_yr_norm + I(PPT_yr_norm^2)  
                          + T_yr_norm + I(T_yr_norm^2)
                          + BALIVE + I(BALIVE^2),
                          data = rdata.scaled, 
                          family = "poisson")

# use non-PIED basal area, AIC = 5640.7
rmodel_zipoiss <- glmmTMB(recruits1 ~ 1
                          + BA.notPIED + I(BA.notPIED^2)
                          + PPT_yr_window_20 + I(PPT_yr_window_20^2)
                          + T_yr_window_20 + I(T_yr_window_20^2)
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)), 
                          ziformula = ~ 1 
                          + PPT_yr_norm + I(PPT_yr_norm^2)  
                          + T_yr_norm + I(T_yr_norm^2)
                          + BA.notPIED + I(BA.notPIED^2), # ns
                          data = rdata.scaled, 
                          family = "poisson")

# used basal area of juniper, AIC = 5643.2
rmodel_zipoiss <- glmmTMB(recruits1 ~ 1
                          + BA.juniper + I(BA.juniper^2)
                          + PPT_yr_window_20 + I(PPT_yr_window_20^2)
                          + T_yr_window_20 + I(T_yr_window_20^2)
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)), 
                          ziformula = ~ 1 
                          + PPT_yr_norm + I(PPT_yr_norm^2)  
                          + T_yr_norm + I(T_yr_norm^2)
                          + BA.juniper + I(BA.juniper^2),
                          data = rdata.scaled, 
                          family = "poisson")


# model reduction...can remove PPT, then T, then BALIVE from Bernoulli, with little change in AIC
# AIC = 4343.7, 4340.9, 5439.1
# remove PPT from Bernoulli
rmodel_zipoiss <- glmmTMB(recruits1 ~ 1
                          + BALIVE + I(BALIVE^2)
                          + PPT_yr_norm + I(PPT_yr_norm^2)
                          + T_yr_norm + I(T_yr_norm^2)
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)),
                          ziformula = ~ 1   
                          + T_yr_norm + I(T_yr_norm^2)
                          + BALIVE + I(BALIVE^2),
                          data = rdata.scaled, 
                          family = "poisson")
summary(rmodel_zipoiss)
res = simulateResiduals(rmodel_zipoiss) 
# deviation ns: p = 0.09159 (without PPT)
# deviation ns: p = 0.09159 (without PPT or T)

# this model is no different from models with predictors on the Bernoulli
# AIC = 4343.4
# residuls do not deviate significanty from expectation (p = 0.0801)
rmodel_zipoiss <- glmmTMB(recruits1 ~ 1
                          + BALIVE + I(BALIVE^2)
                          + PPT_yr_norm + I(PPT_yr_norm^2)
                          + T_yr_norm + I(T_yr_norm^2)
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)),
                          ziformula = ~ 1,
                          data = rdata.scaled, 
                          family = "poisson")
res = simulateResiduals(rmodel_zipoiss)
plot(res, quantreg = T)

# try model with 2-way interacions
# AIC = 4335.0
rmodel_zip.int <- glmmTMB(recruits1 ~ 1
                          + (BALIVE + PPT_yr_norm + T_yr_norm)^2
                          + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) 
                          + offset(log(CENSUS_INTERVAL))
                          + offset(log(PIEDadults1)),
                          ziformula = ~ 1,
                          data = rdata.scaled, 
                          family = "poisson")
res = simulateResiduals(rmodel_zip.int)
plot(res, quantreg = T) # p = 0.07203


# try models with seasonal climate normals
# THIS IS THE BEST MODEL AS OF 9/19/18
# AIC = 4317.7; support for breaking climate into 3 seasons
rmodel_zip3seas <- glmmTMB(recruits1 ~ 1
                       + BALIVE + I(BALIVE^2)
                       + PPT_wd_norm + I(PPT_wd_norm^2)
                       + PPT_m_norm + I(PPT_m_norm^2)
                       + PPT_c_norm + I(PPT_c_norm^2)
                       + T_wd_norm + I(T_wd_norm^2)
                       + T_m_norm + I(T_m_norm^2)
                       + T_c_norm + I(T_c_norm^2)
                       + offset(log(CENSUS_INTERVAL))
                       + offset(log(PIEDadults1)), 
                       ziformula = ~ 1,
                       data = rdata.scaled, 
                       family = "poisson")
res = simulateResiduals(rmodel_zip3seas)
plot(res, quantreg = T) # second best p value so far: p = 0.10428
plot(Effect.glmmTMB("BALIVE", rmodel_zip3seas))
plot(Effect.glmmTMB("PPT_c_norm", rmodel_zip3seas))
plot(Effect.glmmTMB("T_m_norm", rmodel_zip3seas))

#AIC 4319.0; no loss of fit to data when PPT data are 12-month rather than 3 seasons
rmodel_zip3seasT <- glmmTMB(recruits1 ~ 1
                             + BALIVE + I(BALIVE^2)
                             + PPT_yr_norm + I(PPT_yr_norm^2)
                             + T_wd_norm + I(T_wd_norm^2)
                             + T_m_norm + I(T_m_norm^2)
                             + T_c_norm + I(T_c_norm^2)
                             + offset(log(CENSUS_INTERVAL))
                             + offset(log(PIEDadults1)), 
                             ziformula = ~ 1,
                             data = rdata.scaled, 
                             family = "poisson")
res = simulateResiduals(rmodel_zip3seasT)
plot(res, quantreg = T)
# residuals of this model do not differ significantly from expectation (p = 0.07646)
plot(Effect.glmmTMB("BALIVE", rmodel_zip3seasT))
plot(Effect.glmmTMB("PPT_yr_norm", rmodel_zip3seasT))

# this model is worse, AIC = 4334.9
rmodel_zip3seasPPT <- glmmTMB(recruits1 ~ 1
                            + BALIVE + I(BALIVE^2)
                            + T_yr_norm + I(T_yr_norm^2)
                            + PPT_wd_norm + I(PPT_wd_norm^2)
                            + PPT_m_norm + I(PPT_m_norm^2)
                            + PPT_c_norm + I(PPT_c_norm^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(PIEDadults1)), 
                            ziformula = ~ 1,
                            data = rdata.scaled, 
                            family = "poisson")
res = simulateResiduals(rmodel_zip3seasPPT)
plot(res, quantreg = T)
# residuals of this model do not differ significantly from expectation (p = 0.08494)
plot(Effect.glmmTMB("BALIVE", rmodel_zip3seasPPT))
plot(Effect.glmmTMB("PPT_wd_norm", rmodel_zip3seasPPT))
plot(Effect.glmmTMB("PPT_m_norm", rmodel_zip3seasPPT))


# try this model with interactions
# THIS IS THE BEST MODEL AS OF 10/04/18
# AIC = 4297.9
rmodel_zip3seasTint <- glmmTMB(recruits1 ~ 1
                            + (BALIVE + PPT_yr_norm + T_wd_norm + T_c_norm + T_m_norm)^2
                            + I(BALIVE^2) + I(PPT_yr_norm^2)
                            + I(T_wd_norm^2) + I(T_m_norm^2) + I(T_c_norm^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(PIEDadults1)), 
                            ziformula = ~ 1,
                            data = rdata.scaled, 
                            family = "poisson")
res = simulateResiduals(rmodel_zip3seasTint)
plot(res, quantreg = T) # best p value I've seen so far: p=0.14391
plot(Effect.glmmTMB("BALIVE", rmodel_zip3seasTint))
plot(Effect.glmmTMB("PPT_yr_norm", rmodel_zip3seasTint))
plot(Effect.glmmTMB("T_m_norm", rmodel_zip3seasTint))
plot(Effect.glmmTMB("BALIVE:PPT_yr_norm", rmodel_zip3seasTint)) # doesn't work


#plotResiduals(rdata.scaled$BALIVE, res$scaledResiduals, quantreg = T, main = "PREVDIA")
plotResiduals(rdata.scaled$PPT_c_window_25, res$scaledResiduals, quantreg = T, main = "PPT_yr")
plotResiduals(rdata.scaled$T_c_window_25, res$scaledResiduals, quantreg = T, main = "BALIVE")


#### NEW STUFF 11/12/18
### dealing with std'ized covariates

# specify the predictors in the "best" model (or candidate best)
r.predictors <- c("T_yr_norm", "PPT_yr_norm", "BALIVE") # rmodel_zipoiss
# eventually rewrite this so that it can handle alternative "best" models
#r.predictors <- c("T_wd_norm", "T_c_norm", "T_m_norm", "PPT_wd_norm", "PPT_c_norm", "PPT_m_norm", "BALIVE") # rmodel_zip3seas
#r.predictors <- c("T_wd_norm", "T_c_norm", "T_m_norm", "PPT_yr_norm", "BALIVE") # rmodel_zip3seasTint

get_scale = function(data, predictors) {
  sc = list("scale" = NULL, "center"  = NULL)
  for (i in predictors) {
    sc$scale[i] = attributes(data[, i])$"scaled:scale"
    sc$center[i] = attributes(data[, i])$"scaled:center"
  }
  return(sc)
}

r.scaling = get_scale(rdata.scaled, r.predictors)

# remove scaling information from the dataset so that the model doesnt expect scaled data in predict()
for (i in r.predictors) {
  attributes(rdata.scaled[, i]) = NULL
}

# export model for coefficients and scaling information -------------------
save(rmodel_zipoiss, r.scaling, file = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Code/IPM/RecruitRescaling.Rdata")



# ------------------------------------------------------------------------------
# OLD (as of 11/12/18) STUFF FOLLOWS



# model without zero inflation...much much worse
rmodel_nzipoiss <- glmmadmb(recruits ~ 1
                            + BALIVE + I(BALIVE^2)
                            + PPT_c_window_25 + I(PPT_c_window_25^2)
                            + PPT_m_window_25 + I(PPT_m_window_25^2)
                            + T_c_window_25 + I(T_c_window_25^2)
                            + T_m_window_25 + I(T_m_window_25^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(AGB_intra)), 
                            zeroInflation = F,
                            data = rdata.scaled, 
                            family = "poisson")
summary(rmodel_nzipoiss) # Loglike = -2835.95 (with VPD), -2805.03 (with T)

# logLike = -2651.83 (worse), AIC = 5319.7
rmodel_zipoiss.2 <- glmmadmb(recruits ~ 1
                            + BALIVE + I(BALIVE^2)
                            + PPT_yr_norm + I(PPT_yr_norm^2)
                            + T_yr_norm + I(T_yr_norm^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(AGB_intra)), 
                            zeroInflation = TRUE,
                            data = rdata.scaled, 
                            family = "poisson")

# 3 seasons is better than year or cool + monsoon; Loglike = -2613.14 (with T), AIC = 5258.3
rmodel_zipoiss.3 <- glmmadmb(recruits ~ 1
                            + BALIVE + I(BALIVE^2)
                            + PPT_c_window_25 + I(PPT_c_window_25^2)
                            + PPT_wd_window_25 + I(PPT_wd_window_25^2)
                            + PPT_m_window_25 + I(PPT_m_window_25^2)
                            + T_c_window_25 + I(T_c_window_25^2)
                            + T_wd_window_25 + I(T_wd_window_25^2)
                            + T_m_window_25 + I(T_m_window_25^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(AGB_intra)), 
                            zeroInflation = T,
                            data = rdata.scaled, 
                            family = "poisson")


# double-checking that annual climate variables not better - they're not
# AIC = 4687.9
rmodel_zipoiss.4a <- glmmadmb(recruits ~ 1
                             + BALIVE + PPT_yr_norm + T_yr_norm
                             + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2), 
                             zeroInflation = T,
                             data = rdata.scaled, 
                             family = "poisson")
# even with interactions + quadratics, not better; AIC = 4689.4
rmodel_zipoiss.4b <- glmmadmb(recruits ~ 1
                              + (BALIVE + PPT_yr_norm + T_yr_norm)^2
                              + I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2), 
                              zeroInflation = T,
                              data = rdata.scaled, 
                              family = "poisson")

# temp effects are weak - reduce temporal resolution for T variables?
# nope: AIC = 4672.7
rmodel_zipoiss.5 <- glmmadmb(recruits ~ 1
                                 + BALIVE + I(BALIVE^2)
                                 + PPT_c_window_25 + I(PPT_c_window_25^2)
                                 + PPT_wd_window_25 + I(PPT_wd_window_25^2)
                                 + PPT_m_window_25 + I(PPT_m_window_25^2)
                                 + T_yr_norm + I(T_yr_norm^2), 
                                 zeroInflation = T,
                                 data = rdata, 
                                 family = "poisson")
# compare 25-yr vs. 20-yr vs. 15-yr windows
# 25-yr AIC = 4684...now 4651
# 20-yr AIC = 4680...now 4647.9
# 15-yr AIC = 4687.2...now 4653.9
rmodel_zipoiss.5a <- glmmadmb(recruits ~ 1
                                 + BALIVE + I(BALIVE^2)
                                 + PPT_c_window_20 + I(PPT_c_window_20^2)
                                 + PPT_wd_window_20 + I(PPT_wd_window_20^2)
                                 + PPT_m_window_20 + I(PPT_m_window_20^2)
                                 + T_c_window_20 + I(T_c_window_20^2)
                                 + T_wd_window_20 + I(T_wd_window_20^2)
                                 + T_m_window_20 + I(T_m_window_20^2), 
                                 zeroInflation = T,
                                 data = rdata.scaled, 
                                 family = "poisson")
plot(allEffects(rmodel_zipoiss.5a))

# remove BALIVE; AIC = 4655.7
rmodel_zipoiss.5b <- glmmadmb(recruits ~ 1
                              + PPT_c_window_20 + I(PPT_c_window_20^2)
                              + PPT_wd_window_20 + I(PPT_wd_window_20^2)
                              + PPT_m_window_20 + I(PPT_m_window_20^2)
                              + T_c_window_20 + I(T_c_window_20^2)
                              + T_wd_window_20 + I(T_wd_window_20^2)
                              + T_m_window_20 + I(T_m_window_20^2), 
                              zeroInflation = T,
                              data = rdata.scaled, 
                              family = "poisson")

# remove PPT_wd; AIC = 4653.9
rmodel_zipoiss.5c <- glmmadmb(recruits ~ 1
                              + BALIVE + I(BALIVE^2)
                              + PPT_c_window_20 + I(PPT_c_window_20^2)
                              + PPT_m_window_20 + I(PPT_m_window_20^2)
                              + T_c_window_20 + I(T_c_window_20^2)
                              + T_wd_window_20 + I(T_wd_window_20^2)
                              + T_m_window_20 + I(T_m_window_20^2), 
                              zeroInflation = T,
                              data = rdata.scaled, 
                              family = "poisson")

# remove T_c; AIC = 4648.7 (within one AIC unit of best model)
rmodel_zipoiss.5d <- glmmadmb(recruits ~ 1
                              + BALIVE + I(BALIVE^2)
                              + PPT_c_window_20 + I(PPT_c_window_20^2)
                              + PPT_wd_window_20 + I(PPT_wd_window_20^2)
                              + PPT_m_window_20 + I(PPT_m_window_20^2)
                              + T_wd_window_20 + I(T_wd_window_20^2)
                              + T_m_window_20 + I(T_m_window_20^2), 
                              zeroInflation = T,
                              data = rdata.scaled, 
                              family = "poisson")
plot(allEffects(rmodel_zipoiss.5d))


########## OLD STUFF follows

# logLike = -2798.1 (with 25-yr window), -2796.63 (with 20-yr window), -2795.74 (with 15-yr window) 
rmodel_zipoiss.5 <- glmmadmb(recruits ~ 1
                            + (BALIVE + PPT_c_window_15 + PPT_m_window_15 + T_c_window_15 + T_m_window_15)^2
                            + I(BALIVE^2)+ I(PPT_c_window_15^2)+ I(PPT_m_window_15^2) + I(T_m_window_15^2) + I(T_c_window_25^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(AGB_intra)), 
                            zeroInflation = F,
                            data = rdata, 
                            family = "poisson")


rmodel_nzipoiss.6 <- glmmadmb(recruits ~ 1
                            + (BALIVE + PPT_c_window_25 + PPT_m_window_25 + T_c_window_25 + T_m_window_25)^2
                            + I(BALIVE^2)+ I(PPT_c_window_25^2)+ I(PPT_m_window_25^2) + I(T_m_window_25^2) + I(T_c_window_25^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(AGB_intra)), 
                            zeroInflation = F,
                            data = rdata, 
                            family = "poisson")







rdata$ID = as.factor(rdata$plot)
rmodel_zipoissM <- glmmadmb(recruits ~ 1
                            + BALIVE + I(BALIVE^2)
                            + PPT_c_window_25 + I(PPT_c_window_25^2)
                            + PPT_m_window_25 + I(PPT_m_window_25^2)
                            + VPD_c_window_25 + I(VPD_c_window_25^2)
                            + VPD_m_window_25 + I(VPD_m_window_25^2)
                            + offset(log(CENSUS_INTERVAL))
                            + offset(log(AGB_intra))
                            + (1|ID), 
                            zeroInflation = T,
                            data = rdata, 
                            family = "poisson")
summary(rmodel_zipoissM)


########## OLD STUFF

# log transformation of predicters barely supported (delta AIC = 3)
rmodel_zip.logtr <- glmmadmb(recruits ~ 1
                             + log.BALIVE + I(log.BALIVE^2)
                             + log.PPT_c_window_25 + I(log.PPT_c_window_25^2)
                             + log.PPT_m_window_25 + I(log.PPT_m_window_25^2)
                             + T_c_window_25 + I(T_c_window_25^2)
                             + T_m_window_25 + I(T_m_window_25^2)
                             + offset(log(CENSUS_INTERVAL))
                             + offset(log(AGB_intra)), 
                             zeroInflation = TRUE,
                             data = rdata.scaled, 
                             family = "poisson")
summary(rmodel_ziptr) # AIC = 5265.8, Loglike = -2620.88 (with T)
# what improvement there was, was due to log.BALIVE, not log.PPT
rmodel_zip.logtr2 <- glmmadmb(recruits ~ 1
                              + log.BALIVE + I(log.BALIVE^2)
                              + PPT_c_window_25 + I(PPT_c_window_25^2)
                              + PPT_m_window_25 + I(PPT_m_window_25^2)
                              + T_c_window_25 + I(T_c_window_25^2)
                              + T_m_window_25 + I(T_m_window_25^2)
                              + offset(log(CENSUS_INTERVAL))
                              + offset(log(AGB_intra)), 
                              zeroInflation = TRUE,
                              data = rdata.scaled, 
                              family = "poisson")





# Model: elevation + BA
rmodel <- glm(data = rdata, recruits ~ elev + BALIVE +
                offset(log(CENSUS_INTERVAL)) + offset(log(AGB_intra + 0.00001)), 
              family = "poisson")
summary(rmodel)
save(rmodel, file = "D:/EvansLab/Final/Models/B/recr.rda")
pdf("D:/EvansLab/Final/Output/B/RecruitmentModel.pdf")
coefplot(rmodel, intercept = F)
hist(rdata$recruits, col = rgb(1,0,0,0.75), main = "Observed (red) vs. modeled (green) recruitment count", xlab = "Number of recruits")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
png("D:/EvansLab/Final/Manuscript/FigS8.png")
hist(rdata$recruits, col = rgb(1,0,0,0.75), xlab = "Number of recruits", main = "")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()

# Model: elevation + climate
rmodel <- glm(data = rdata, recruits ~ elev + PPT_c_window_25 + PPT_w_window_25 + VPD_c_window_25 + VPD_w_window_25 +
                offset(log(measInterval)) + offset(log(AGB_intra + 0.00001)), 
              family = "poisson")
rmodel <- update(rmodel, . ~ . - VPD_c_window_25)
# rmodel <- update(rmodel, . ~ . + I(VPD_w_window_25))
rmodel <- update(rmodel, . ~ . - VPD_w_window_25)
rmodel <- update(rmodel, . ~ . - elev)
summary(rmodel)
save(rmodel, file = "D:/EvansLab/Final/Models/C/recr.rda")
pdf("D:/EvansLab/Final/Output/C/RecruitmentModel.pdf")
coefplot(rmodel, intercept = F)
hist(rdata$recruits, col = rgb(1,0,0,0.75), main = "Observed (red) vs. modeled (green) recruitment count", xlab = "Number of recruits")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
coefficientplot <- coefplot(rmodel, intercept = F) + scale_y_discrete(labels = c("Cool-season\nprecipitation", "Warm-season\nprecipitation")) +
  ggtitle("")
ggsave("D:/EvansLab/Final/Manuscript/FigS7.png", coefficientplot, scale = 1.5)
# Model: elevation + BA + climate
rmodel <- glm(data = rdata, recruits ~ elev + baLive + PPT_c_window_25 + PPT_w_window_25 + VPD_c_window_25 + VPD_w_window_25 +
                offset(log(measInterval)) + offset(log(AGB_intra + 0.00001)), 
              family = "poisson")
rmodel <- update(rmodel, . ~ . - VPD_w_window_25)
rmodel <- update(rmodel, . ~ . - VPD_c_window_25)
rmodel <- update(rmodel, . ~ . - elev)
rmodel <- update(rmodel, . ~ . - baLive)
summary(rmodel)
save(rmodel, file = "D:/EvansLab/Final/Models/BC/recr.rda")
pdf("D:/EvansLab/Final/Output/BC/RecruitmentModel.pdf")
coefplot(rmodel, intercept = F)
hist(rdata$recruits, col = rgb(1,0,0,0.75), main = "Observed (red) vs. modeled (green) recruitment count", xlab = "Number of recruits")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
png("D:/EvansLab/Final/Manuscript/FigS9.png")
hist(rdata$recruits, col = rgb(1,0,0,0.75), xlab = "Number of recruits", main = "")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
