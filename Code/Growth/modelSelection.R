library(ggeffects)
library(ggplot2)
library(coefplot)
library(cowplot)
library(effects)

# Read data
# start with the same data file that is used for analysis of survival/mortality
grdata <- read.csv("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

# Only keep trees that didn't die
grdata <- subset(grdata, STATUSCD == 1) #15742

# Create increment columns
# note that growth increments need to be moved to the positive realm (by adding a constant)
# IF log transform is used
grdata$AGB_INCR <- grdata$DRYBIO_AG_DIFF / grdata$CENSUS_INTERVAL
grdata$DIA_INCR <- grdata$DIA_DIFF / grdata$CENSUS_INTERVAL
grdata$BA_INCR <- grdata$BA_DIFF / grdata$CENSUS_INTERVAL

# distribution of size
#hist(grdata$BAt1) # size is lognormal-ish, with a wonky bit at the small end of the scale (due to min size threshold?)
#hist(grdata$BAt1, breaks = c(seq(0, 1220, by = 10)), xlim = c(0, 500))
#hist(log(grdata$BAt1)) # not quite lognormal, but worth trying in the models below (heavy in the left tail)
hist(grdata$PREVDIA)

# distribution of other predictors
hist(grdata$PPT_yr) # not too bad...Poisson-ish but with a large mean count
hist(log(grdata$PPT_yr)) # more normal-looking
hist(grdata$T_yr) # looks normal
hist(grdata$BALIVE) # not too bad...Poisson-ish but with a large mean count
hist(log(grdata$BALIVE)) # log transform has a heavy left tail

# examine distribution of response(s)
#hist(grdata$AGB_INCR)
#hist(grdata$AGB_INCR, breaks = c(seq(-106, 165, by = 0.5)), xlim = c(-10, 15))
hist(grdata$DIA_INCR)
summary(grdata$DIA_INCR)
hist(grdata$DIA_INCR, breaks = c(seq(-2.5, 1.75, by = 0.01)), xlim = c(-0.5, 0.5))
# should the most extreme diameter increments be filtered out of the data set?
# are the most extreme values biologically plausible?

#hist(grdata$BA_INCR)
#hist(grdata$BA_INCR, breaks = c(seq(-76, 53, by = 0.1)), xlim = c(-5, 5))

# transformed predictors
grdata$log.PPT_yr <- log(grdata$PPT_yr)
grdata$log.size <- log(grdata$BAt1)
grdata$log.BALIVE <- log(grdata$BALIVE)

# examine largest trees
#grdata.2 <- subset(grdata, BAt1 < 500) # 15704, removed 38 cases
lgst.trees <- subset(grdata, BAt1 >500)

# eliminate 69 extreme cases (out of 15,742)
grdata.2 <- subset(grdata, DIA_INCR < 0.4 & DIA_INCR > -0.4) # 15639; 65 cases eliminated
extreme.incr <- subset(grdata, DIA_INCR > 0.4 | DIA_INCR < -0.4)
# slightly more conservative would be the following:
# grdata <- subset(grdata, BA_INCR < 5.0 & BA_INCR >-5.0) 

# standardize covariates
library(dplyr)
grdata.scaled <- grdata %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                    -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                    -CENSUS_INTERVAL,
                                                    -AGB_INCR, -DIA_INCR, -BA_INCR))

grdata2.scaled <- grdata.2 %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                          -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                          -CENSUS_INTERVAL,
                                                          -AGB_INCR, -DIA_INCR, -BA_INCR))

library(lme4)
library(lmerTest)
library(MuMIn) # use MuMin to choose between models (AICc)
library(DHARMa) # use DHARMa to check residuals

# compare different response variables
# DIA, AGB, BA
gmodel.AGB <- lmer(AGB_INCR ~ PREV_DRYBIO_AG + I(PREV_DRYBIO_AG^2) + BALIVE + PPT_c + PPT_wd + PPT_m + VPD_c + VPD_wd + VPD_m + (1|PLT_CN), data = grdata.scaled)
gmodel.DIA <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c + PPT_wd + PPT_m + VPD_c + VPD_wd + VPD_m + (1|PLT_CN), data = grdata.scaled)
gmodel.BA <- lmer(BA_INCR ~ BAt1 + I(BAt1^2) + BALIVE + PPT_c + PPT_wd + PPT_m + VPD_c + VPD_wd + VPD_m + (1|PLT_CN), data = grdata.scaled)
mod.comp <- model.sel(gmodel.AGB, gmodel.DIA, gmodel.BA)

# check residuals
class(gmodel.AGB) <- "lmerMod"
plot(simulateResiduals(gmodel.AGB, integerResponse = F), quantreg = T)
plot(gmodel.AGB)

class(gmodel.DIA) <- "lmerMod"
plot(simulateResiduals(gmodel.DIA, integerResponse = F), quantreg = T)
plot(gmodel.DIA)

class(gmodel.BA) <- "lmerMod"
plot(simulateResiduals(gmodel.BA, integerResponse = F), quantreg = T)
plot(gmodel.BA)

### LISA: only the residuals of the DIA model look more or less acceptable, for BA and AGB we would require a non-linear model, I suppose. 

gmodel.1a <- lmer(DIA_INCR ~ BAt1 + I(BAt1^2) + BALIVE + PPT_yr + VPD_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.1b <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + VPD_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.1c <- lmer(DIA_INCR ~ PREV_DRYBIO_AG + I(PREV_DRYBIO_AG^2) + BALIVE + PPT_yr + VPD_yr + (1|PLT_CN), data = grdata.scaled)

class(gmodel.1a) <- "lmerMod"
plot(simulateResiduals(gmodel.1a, integerResponse = F), quantreg = T)

# residuals looks good for model using PREVDIA as a (size) predictor
class(gmodel.1b) <- "lmerMod"
plot(simulateResiduals(gmodel.1b, integerResponse = F), quantreg = T)

# compare annual vs. 3 vs. 4 seasons...likes annual better
gmodel.1a <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + VPD_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.1b <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c + PPT_wd + PPT_m + VPD_c + VPD_wd + VPD_m + (1|PLT_CN), data = grdata.scaled)
gmodel.1c <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c + PPT_pf + PPT_fs + PPT_m + VPD_c + VPD_pf + VPD_fs + VPD_m + (1|PLT_CN), data = grdata.scaled)
mod.comp1 <- model.sel(gmodel.1a, gmodel.1b, gmodel.1c)

# compare T vs. VPD...delta AIC is 2.94 (slight preference for T)
gmodel.2a <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + VPD_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.2b <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + T_yr + (1|PLT_CN), data = grdata.scaled)
mod.comp2 <- model.sel(gmodel.2a, gmodel.2b)

# compare normals vs. census interval vs. anomalies...census interval is best
gmodel.3a <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + T_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.3b <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c + PPT_wd + PPT_m + T_c + T_wd + T_m + (1|PLT_CN), data = grdata.scaled)
gmodel.3c <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr_norm + T_yr_norm + (1|PLT_CN), data = grdata.scaled)
gmodel.3d <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + T_c_norm + T_wd_norm + T_m_norm + (1|PLT_CN), data = grdata.scaled)
gmodel.3e <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr_anom + T_yr_anom + (1|PLT_CN), data = grdata.scaled)
gmodel.3f <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c_anom + PPT_wd_anom + PPT_m_anom + T_c_anom + T_wd_anom + T_m_anom + (1|PLT_CN), data = grdata.scaled)
gmodel.3g <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr_anom + T_yr_anom + PPT_yr_norm + T_yr_norm + (1|PLT_CN), data = grdata.scaled)
gmodel.3h <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + T_c_norm + T_wd_norm + T_m_norm + PPT_c_anom + PPT_wd_anom + PPT_m_anom + T_c_anom + T_wd_anom + T_m_anom + (1|PLT_CN), data = grdata.scaled)
gmodel.3i <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPTex_yr_anom + Tex_yr_anom + (1|PLT_CN), data = grdata.scaled)
gmodel.3j <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + Tex_c_anom + Tex_wd_anom + Tex_m_anom + (1|PLT_CN), data = grdata.scaled)
gmodel.3k <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPTex_yr_anom + Tex_yr_anom + PPT_yr_norm + T_yr_norm + (1|PLT_CN), data = grdata.scaled)
gmodel.3l <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + T_c_norm + T_wd_norm + T_m_norm + PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + Tex_c_anom + Tex_wd_anom + Tex_m_anom + (1|PLT_CN), data = grdata.scaled)

mod.comp3 <- model.sel(gmodel.3a, gmodel.3b, gmodel.3c, gmodel.3d, gmodel.3e, gmodel.3f, 
                       gmodel.3g, gmodel.3h, gmodel.3i, gmodel.3j, gmodel.3k, gmodel.3l)
# gmodel.3a is best among these (true whether BAt1 or PREVDIA is the size predictor)

# add 2-way interactions, excluding size^2
# clear support for 2-way interactions
gmodel.4a <- lmer(DIA_INCR ~ (PREVDIA + PREVDIA + PPT_yr + T_yr)^2 + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.4b <- lmer(DIA_INCR ~ (PREVDIA + PREVDIA + PPT_c + PPT_wd + PPT_m + T_c + T_wd + T_m)^2 + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.4c <- lmer(DIA_INCR ~ (PREVDIA + PREVDIA + PPT_c + PPT_pf + PPT_fs + PPT_m + T_c + T_pf + T_fs + T_m)^2 + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.4d <- lmer(DIA_INCR ~ (PREVDIA + PREVDIA + PPT_yr + VPD_yr)^2 + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.4e <- lmer(DIA_INCR ~ (PREVDIA + PREVDIA + PPT_c + PPT_wd + PPT_m + VPD_c + VPD_wd + VPD_m)^2 + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.4f <- lmer(DIA_INCR ~ (PREVDIA + PREVDIA + PPT_c + PPT_pf + PPT_fs + PPT_m + VPD_c + VPD_pf + VPD_fs + VPD_m)^2 + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.4g <- lmer(DIA_INCR ~ (PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + VPD_yr)^2 + (1|PLT_CN), data = grdata.scaled)
gmodel.4h <- lmer(DIA_INCR ~ (PREVDIA + I(PREVDIA^2) + BALIVE + PPT_yr + T_yr)^2 + (1|PLT_CN), data = grdata.scaled)
mod.comp4 <- model.sel(gmodel.3a, gmodel.4a, gmodel.4b, gmodel.4c,
                       gmodel.4d, gmodel.4e, gmodel.4f, gmodel.4g, gmodel.4h)
# when no data are filtered (DIA_INCR)...
# likes models 4g and 4h best, after 3a

class(gmodel.3a) <- "lmerMod"
plot(simulateResiduals(gmodel.3a, integerResponse = F), quantreg = T) # residuals look good

# try model with all quadratics instead of 2-way interactions
gmodel.5 <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                            PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                            (1|PLT_CN), data = grdata.scaled)
growSD <- sd(resid(gmodel.5))
mod.comp5 <- model.sel(gmodel.3a, gmodel.4g, gmodel.4h, gmodel.5) # model 5 best, by 15.57 (model 3)
gmodel2.5 <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                   PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                   (1|PLT_CN), data = grdata2.scaled)

gmodel3.5 <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + BALIVE + I(BALIVE^2) +
                   PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                   (1|PLT_CN), data = grdata3.scaled)
plot(allEffects(gmodel.5)) # effects look the same as they did with BAt1 as size predictor
plot(allEffects(gmodel2.5))
plot(allEffects(gmodel3.5))
 
# log transform of BALIVE?
gmodel.5logBALIVE <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + log.BALIVE + I(log.BALIVE^2) +
                          PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                          (1|PLT_CN), data = grdata.scaled)

gmodel2.5logBALIVE <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + log.BALIVE + I(log.BALIVE^2) +
                           PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                           (1|PLT_CN), data = grdata2.scaled)

gmodel3.5logBALIVE <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + log.BALIVE + I(log.BALIVE^2) +
                           PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                           (1|PLT_CN), data = grdata3.scaled)
plot(allEffects(gmodel.5logBALIVE))
plot(allEffects(gmodel2.5logBALIVE))
plot(allEffects(gmodel3.5logBALIVE))


mod.comp5 <- model.sel(gmodel.5, gmodel.5logBALIVE) # equal with log.BALIVE
mod2.comp5 <- model.sel(gmodel2.5, gmodel2.5logBALIVE) # log transformation of BALIVE is worse (19.71)
mod3.comp5 <- model.sel(gmodel3.5, gmodel3.5logBALIVE) # 33.42 AIC units from model with no log transformations
plot(allEffects(gmodel.5logsize))
plot(allEffects(gmodel2.5logsize))
plot(allEffects(gmodel3.5logsize))


mod.comp5 <- model.sel(gmodel.3a, gmodel.4g, gmodel.4h, gmodel.5) 
# models 4g and 4h better with unfiltered data
# model 5 best with filtered data
# most support for quadratic terms on size and BALIVE

res = simulateResiduals(gmodel.5, integerResponse = F) #doesn't work??
plotResiduals(grdata.scaled$PREVDIA, res$scaledResiduals, quantreg = T, main = "PREVDIA")

# 2 quadratic terms (size and BALIVE), but not in interactions
gmodel.5a <- lmer(DIA_INCR ~ (PREVDIA + BALIVE + PPT_yr + VPD_yr)^2 + I(PREVDIA^2) + I(BALIVE^2) + (1|PLT_CN), data = grdata.scaled)
#plot(allEffects(gmodel.5b))
gmodel.5c <- lmer(DIA_INCR ~ (PREVDIA + BALIVE + PPT_yr + VPD_yr)^2 + I(PREVDIA^2) + I(PREVDIA^3) + I(BALIVE^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.5e <- lmer(DIA_INCR ~ (PREVDIA + BALIVE + PPT_yr + VPD_yr)^2 + I(PREVDIA^2) + I(BALIVE^2) + 
                    I(PREVDIA^2):BALIVE + I(PREVDIA^2):PPT_yr + I(PREVDIA^2):VPD_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.5f <- lmer(DIA_INCR ~ (PREVDIA + BALIVE + PPT_yr + VPD_yr)^2 + I(PREVDIA^2) + I(BALIVE^2) + 
                    I(BALIVE^2):BAt1 + I(BALIVE^2):PPT_yr + I(BALIVE^2):VPD_yr + (1|PLT_CN), data = grdata.scaled)
gmodel.5g <- lmer(DIA_INCR ~ (PREVDIA + BALIVE + PPT_yr + VPD_yr)^2 + I(PREVDIA^2) + I(PREVDIA^3) + I(BALIVE^2) + 
                    I(PREVDIA^2):BALIVE + I(PREVDIA^2):PPT_yr + I(PREVDIA^2):VPD_yr + (1|PLT_CN), data = grdata.scaled)

mod.comp5 <- model.sel(gmodel.3a, gmodel.4g, gmodel.4h, gmodel.5, gmodel.5a, 
                       gmodel.5c, gmodel.5e, gmodel.5f, 
                       gmodel.5g) 
# model 5g best with unfiltered data (8/7/18)
# model 5 best with filtered data (DIA_INCR extremes removed)

gmodel.6 <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + I(PREVDIA^3) + BALIVE + I(BALIVE^2) +
                   PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                   (1|PLT_CN), data = grdata.scaled)
mod.comp6 <- model.sel(gmodel.3a, gmodel.4g, gmodel.4h, gmodel.5, gmodel.5a, 
                       gmodel.5c, gmodel.5e, gmodel.5f, 
                       gmodel.5g, gmodel.6) 
plot(allEffects(gmodel.6), cex = 0.5) 

gmodel.7 <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + I(PREVDIA^3) + BALIVE + I(BALIVE^2) +
                   PPT_yr + I(PPT_yr^2) + T_yr + I(T_yr^2) + 
                   (1|PLT_CN), data = grdata.scaled)
mod.comp7 <- model.sel(gmodel.5, gmodel.6, gmodel.7)
gmodel.8 <- lmer(DIA_INCR ~ PREVDIA + I(PREVDIA^2) + I(PREVDIA^3) + BALIVE + I(BALIVE^2) +
                   PPT_yr + I(PPT_yr^2) + 
                   (1|PLT_CN), data = grdata.scaled)
mod.comp8 <- model.sel(gmodel.5, gmodel.6, gmodel.7, gmodel.8)
gmodel.9 <- lmer(DIA_INCR ~ (PREVDIA + BALIVE + PPT_yr)^2 +
                            I(PREVDIA^2) + I(PREVDIA^3) +  I(BALIVE^2) +
                            I(PPT_yr^2) + 
                            (1|PLT_CN), data = grdata.scaled)
mod.comp9 <- model.sel(gmodel.5, gmodel.6, gmodel.7, gmodel.8, gmodel.9)

# after log-transform of BALIVE, model 5i is best (simplest model: main effects + size^2)
gmodel.5h <- lmer(DIA_INCR ~ (BAt1 + log.BALIVE + PPT_yr + VPD_yr)^2 + I(BAt1^2) + I(log.BALIVE^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.5i <- lmer(DIA_INCR ~ BAt1 + log.BALIVE + PPT_yr + VPD_yr + I(BAt1^2) + (1|PLT_CN), data = grdata.scaled)
gmodel.5j <- lmer(DIA_INCR ~ BAt1 + I(BAt1^2) + log.BALIVE + I(log.BALIVE^2) +
                   PPT_yr + I(PPT_yr^2) + VPD_yr + I(VPD_yr^2) + 
                   (1|PLT_CN), data = grdata.scaled)
mod.comp5 <- model.sel(gmodel.3a, gmodel.4g, gmodel.4h, gmodel.5, gmodel.5a, 
                       gmodel.5b, gmodel.5c, gmodel.5d, gmodel.5e, gmodel.5f, 
                       gmodel.5g, gmodel.5h, gmodel.5i, gmodel.5j) 


#### NEW STUFF 11/07/18
### dealing with std'ized covariates

# specify the predictors in the "best" model (or candidate best)
gr.predictors <- c("PREVDIA", "VPD_yr", "PPT_yr", "BALIVE") 
# eventually rewrite this so that it can handle alternative "best" models

get_scale = function(data, predictors) {
  sc = list("scale" = NULL, "center"  = NULL)
  for (i in predictors) {
    sc$scale[i] = attributes(data[, i])$"scaled:scale"
    sc$center[i] = attributes(data[, i])$"scaled:center"
  }
  return(sc)
}

gr.scaling = get_scale(grdata.scaled, gr.predictors)

# remove scaling information from the dataset so that the model doesnt expect scaled data in predict()
for (i in gr.predictors) {
  attributes(grdata.scaled[, i]) = NULL
}

# export model for coefficients and scaling information -------------------
save(gmodel.5, gr.scaling, growSD, file = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Code/IPM/GrRescaling.Rdata")





#### OLD STUFF...Michiel
gmodel <- update(gmodel, . ~ . - PPT_c)
gmodel <- update(gmodel, . ~ . - VPD_w)
gmodel <- update(gmodel, . ~ . - VPD_c)
gmodel <- update(gmodel, . ~ . + I(DRYBIO_AG^2))
# gmodel <- update(gmodel, . ~ . + I(PPT_w^2))
summary(gmodel)
growSD <- sd(resid(gmodel))
save(gmodel, growSD, file = "D:/EvansLab/Final/Models/BC/grow.rda")
pdf("D:/EvansLab/Final/Output/BC/GrowthModel.pdf")
coefplot(gmodel, intercept = F)
preds <- ggpredict(gmodel, terms = c("DRYBIO_AG"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. biomass")
preds <- ggpredict(gmodel, terms = c("elev"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. elevation")
preds <- ggpredict(gmodel, terms = c("baLive"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. plot BA")
preds <- ggpredict(gmodel, terms = c("PPT_w"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. PPT warm")
dev.off()
# Publication figures
preds <- ggpredict(gmodel, terms = c("DRYBIO_AG"))
A <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Biomass")
preds <- ggpredict(gmodel, terms = c("elev"))
B <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Elevation")
preds <- ggpredict(gmodel, terms = c("baLive"))
C <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Plot basal area")
preds <- ggpredict(gmodel, terms = c("PPT_w"))
D <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Warm-season precipitation")
all <- plot_grid(A, B, C, D, labels = c("A", "B", "C", "D"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/FigS3.png", all, base_aspect_ratio = 1.5)

### reduced models follow (without climate, without competition)
# Model: size + elevation + BA
gmodel <- lm(data = grdata, AGB_INCR ~ DRYBIO_AG + I(DRYBIO_AG^2) + BALIVE)
summary(gmodel)
growSD <- sd(resid(gmodel))
save(gmodel, growSD, file = "D:/EvansLab/Final/Models/B/grow.rda")
pdf("D:/EvansLab/Final/Output/B/GrowthModel.pdf")
coefplot(gmodel, intercept = F)
preds <- ggpredict(gmodel, terms = c("DRYBIO_AG"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. biomass")
preds <- ggpredict(gmodel, terms = c("elev"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. elevation")
preds <- ggpredict(gmodel, terms = c("baLive"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. plot BA")
dev.off()
preds <- ggpredict(gmodel, terms = c("DRYBIO_AG"))
A <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Biomass")
preds <- ggpredict(gmodel, terms = c("elev"))
B <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Elevation")
preds <- ggpredict(gmodel, terms = c("baLive"))
C <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Plot basal area")
all <- plot_grid(A, B, C, labels = c("A", "B", "C"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/FigS1.png", all, base_aspect_ratio = 1.5)

# Model: size + elevation + climate
gmodel <- lm(data = grdata, GROWTH_INCR ~ DRYBIO_AG + elev + PPT_c + PPT_w + VPD_c + VPD_w)
gmodel <- update(gmodel, . ~ . - PPT_c)
gmodel <- update(gmodel, . ~ . - VPD_w)
gmodel <- update(gmodel, . ~ . + I(DRYBIO_AG^2))
# gmodel <- update(gmodel, . ~ . + I(VPD_c^2))
# gmodel <- update(gmodel, . ~ . + I(VPD_c^3))
# gmodel <- update(gmodel, . ~ . + I(VPD_c^4))
gmodel <- update(gmodel, . ~ . - VPD_c) # VPD_c was finally removed because couldn't find a good shape where growth decreased with increasing VPD
# gmodel <- update(gmodel, . ~ . + VPD_w) # No improvement
# gmodel <- update(gmodel, . ~ . + I(PPT_w^2)) # Not significant
gmodel <- update(gmodel, . ~ . - elev)
summary(gmodel)
growSD <- sd(resid(gmodel))
save(gmodel, growSD, file = "D:/EvansLab/Final/Models/C/grow.rda")
pdf("D:/EvansLab/Final/Output/C/GrowthModel.pdf")
coefplot(gmodel, intercept = F)
preds <- ggpredict(gmodel, terms = c("DRYBIO_AG"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. biomass")
# preds <- ggpredict(gmodel, terms = c("elev"))
# ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. elevation")
preds <- ggpredict(gmodel, terms = c("PPT_w"))
ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. PPT warm")
# preds <- ggpredict(gmodel, terms = c("VPD_c"))
# ggplot(preds, aes(x, predicted)) + geom_line() + ggtitle("Growth vs. VPD cool")
dev.off()
preds <- ggpredict(gmodel, terms = c("DRYBIO_AG"))
A <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Biomass")
preds <- ggpredict(gmodel, terms = c("PPT_w"))
B <- ggplot(preds, aes(x, predicted)) + geom_line() + ylab("Growth") + xlab("Warm-season precipitation")
all <- plot_grid(A, B, labels = c("A", "B"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/FigS2.png", all, base_aspect_ratio = 1.5)

