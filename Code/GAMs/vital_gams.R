#### GAMs for PIED vital rate models

library(mgcv)
library(DHARMa) # use DHARMa to check residuals
library(tidyverse)

## Growth data
# Read data
# start with the same data file that is used for analysis of survival/mortality
grdata <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

# Only keep trees that didn't die
grdata <- subset(grdata, STATUSCD == 1) #18204

# Create increment columns
# note that growth increments need to be moved to the positive realm (by adding a constant)
# IF log transform is used
grdata$AGB_INCR <- grdata$DRYBIO_AG_DIFF / grdata$CENSUS_INTERVAL
grdata$DIA_INCR <- grdata$DIA_DIFF / grdata$CENSUS_INTERVAL
grdata$BA_INCR <- grdata$BA_DIFF / grdata$CENSUS_INTERVAL

grdata=grdata[-which(grdata$DIA_INCR<(-1)),]

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

# standardize covariates
library(dplyr)
grdata.scaled <- grdata %>% mutate_at(scale, .vars = vars(-CN, -PREV_TRE_CN, -PLT_CN, -PREV_PLT_CN, -CONDID,
                                                          -STATUSCD, -MEASYEAR, -PREV_MEASYEAR, 
                                                          -CENSUS_INTERVAL,
                                                          -AGB_INCR, -DIA_INCR, -BA_INCR))

## Survival data
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

## Recruitment data
# Read and process data
rdata <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)

# looking at the data
recruit.table <- table(rdata[, c("PIEDadults1", "recruits1")])
colSums(recruit.table)
table(rdata[, c("PIEDadults4", "recruits1")])

plot.info <- c("plot", "lat", "lon", "elev", "state", "county", "plotID", "PIEDadults1", "recruits1")
colonize.plots <- rdata[rdata$PIEDadults1 == 0 & rdata$recruits1 > 0, plot.info]

# Subset data: only plots that contained PIED trees at time 1 census
# AGB: note that it was previously incorrectly time 2
rdata <- subset(rdata, PIEDadults1 > 0)
#rdata <- subset(rdata, PIEDadults8 > 0)

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

# compare BA.all to BALIVE...they should be the same
# both in sq feet, but BALIVE is on a per acre basis
# BA.all is simply the sum per plot (plot = 672.5 m^2 or 0.16617837 acre)
# expansion factor is 6.01763*BA.all
rdata$BA.all <- rdata$BA.PIED + rdata$BA.notPIED
plot(rdata$BALIVE, 6.01763*rdata$BA.all) # mostly fall on/near one:one, but BALIVE < BA.all

#rdata3 <- subset(rdata, BALIVE > 0) # eliminate data points with no live trees at time 1
# 6635 plots

# standardize covariates
rdata.scaled <- rdata %>% mutate_at(scale, .vars = vars(-plot, -lat, -lon, -elev, -PApied,
                                                        -state, -county, -plotID, -CONDID, 
                                                        -measYear, -plotPrev, -PREV_MEASYEAR,
                                                        -CENSUS_INTERVAL, -recruits1, -recruits12,
                                                        -AGB_intra, -BA.PIED, -PIEDadults1,
                                                        -PIEDadults4, -PIEDadults8, -cumDIA.PIED))

## Growth models
grdata.scaled$PLT_CN_factor<-as.factor(grdata.scaled$PLT_CN)
k=5
gmodel.clim.gam <- bam(DIA_INCR ~ s(PREVDIA,k=k) + s(PPT_yr_norm,k=k) + s(T_yr_norm,k=k) + 
                          s(PLT_CN_factor,bs = "re"), data = grdata.scaled)
gmodel.clim.int.gam <- bam(DIA_INCR ~ s(PREVDIA,k=k) + te(PPT_yr_norm,T_yr_norm,k=k) + 
                         s(PLT_CN_factor,bs = "re"), data = grdata.scaled)
gmodel.clim.comp.gam <- bam(DIA_INCR ~ s(PREVDIA,k=k) + s(BALIVE,k=k) + s(PPT_yr_norm,k=k) + s(T_yr_norm,k=k) + 
                               s(PLT_CN_factor,bs = "re"), data = grdata.scaled)
gmodel.int.gam <- bam(DIA_INCR ~ s(PREVDIA,k=k) + s(BALIVE,k=k) + te(PPT_yr_norm,T_yr_norm,k=k) +  
                        s(PLT_CN_factor,bs = "re"), data = grdata.scaled)

res = simulateResiduals(gmodel.clim.gam, integerResponse = F)
plot(res, quantreg = T) #p = 0

res = simulateResiduals(gmodel.clim.comp.gam, integerResponse = F)
plot(res, quantreg = T) #p = 0

res = simulateResiduals(gmodel.int.gam, integerResponse = F)
plot(res, quantreg = T) #p = 0

## Survival models
survData3.scaled$PLT_CN_factor<-as.factor(survData3.scaled$PLT_CN)
survData.scaled$PLT_CN_factor<-as.factor(survData.scaled$PLT_CN)
k=8
smodel.clim.gam <- bam(mort ~ s(PREVDIA,k=3) + s(PPT_yr_norm,k=k) + s(T_yr_norm,k=5)+ 
                         s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                       family = binomial(link = cloglog), data = survData3.scaled)
smodel.clim.int.gam <- bam(mort ~ s(PREVDIA,k=3) + te(PPT_yr_norm,T_yr_norm,k=5)+ 
                             s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                           family = binomial(link = cloglog), data = survData3.scaled)
smodel.clim.comp.gam <- bam(mort ~ s(PREVDIA,k=3) + s(BALIVE,k=k) + s(PPT_yr_norm,k=k) + s(T_yr_norm,k=5)+ 
                              s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                            family = binomial(link = cloglog), data = survData3.scaled)
smodel.int.gam <- bam(mort ~ s(PREVDIA,k=3) + s(BALIVE,k=k) + te(PPT_yr_norm,T_yr_norm,k=5)+ 
                        s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                      family = binomial(link = cloglog), data = survData3.scaled)

smodel.clim.gam <- bam(mort ~ s(PREVDIA,k=3) + s(PPT_yr_norm,k=5) +
                             T_yr_norm + 
                             s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                           family = binomial(link = cloglog), data = survData3.scaled)
smodel.clim.int.gam <- bam(mort ~ s(PREVDIA,k=3) + s(PPT_yr_norm,k=5) +
                        T_yr_norm + T_yr_norm:PPT_yr_norm +
                        s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                      family = binomial(link = cloglog), data = survData3.scaled)
smodel.clim.comp.gam <- bam(mort ~ s(PREVDIA,k=3) + s(BALIVE,k=k) + s(PPT_yr_norm,k=5) +
                        T_yr_norm + 
                        s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                      family = binomial(link = cloglog), data = survData3.scaled)
smodel.int.gam <- bam(mort ~ s(PREVDIA,k=3) + s(BALIVE,k=k) + s(PPT_yr_norm,k=5) +
                      T_yr_norm + T_yr_norm:PPT_yr_norm +
                        s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                      family = binomial(link = cloglog), data = survData3.scaled)

smodel.clim.fire.gam <- bam(mort ~ s(PREVDIA,k=3) + s(PPT_yr_norm,k=k) + T_yr_norm + 
                              s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                            family = binomial(link = cloglog), data = survData.scaled)
smodel.clim.int.fire.gam <- bam(mort ~ s(PREVDIA,k=3) + s(PPT_yr_norm) + T_yr_norm + T_yr_norm:PPT_yr_norm +
                                  s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                                family = binomial(link = cloglog), data = survData.scaled)
smodel.clim.comp.fire.gam <- bam(mort ~ s(PREVDIA,k=3) + s(BALIVE,k=k) + s(PPT_yr_norm,k=k) + T_yr_norm + 
                                   s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                                 family = binomial(link = cloglog), data = survData.scaled)
smodel.int.fire.gam <- bam(mort ~ s(PREVDIA,k=3) + s(BALIVE,k=k) + s(PPT_yr_norm) + T_yr_norm + T_yr_norm:PPT_yr_norm +
                             s(PLT_CN_factor,bs = "re") + offset(log(CENSUS_INTERVAL)),
                           family = binomial(link = cloglog), data = survData.scaled)

res = simulateResiduals(smodel.clim.gam10)
plot(res, quantreg = T) #p = 0

res = simulateResiduals(smodel.clim.comp.gam)
plot(res, quantreg = T) #p = 0

res = simulateResiduals(smodel.int.gam)
plot(res, quantreg = T) #p = 0

## Recriutment models
rdata.scaled$off<-log(rdata.scaled$CENSUS_INTERVAL) + log(rdata.scaled$PIEDadults1)
k=5
rmodel.clim.gam <- gam(list(recruits1 ~ s(PPT_yr_norm,k=k) + s(T_yr_norm,k=k)+ 
                              offset(off),~offset(off)),family = ziplss(), data = rdata.scaled)
rmodel.clim.int.gam <- gam(list(recruits1 ~ te(PPT_yr_norm,T_yr_norm,k=k)+
                              offset(off),~offset(off)),family = ziplss(), data = rdata.scaled)
rmodel.clim.comp.gam <- gam(list(recruits1 ~ s(BALIVE,k=k) + s(PPT_yr_norm,k=k) + s(T_yr_norm,k=k)+ 
                              offset(off),~offset(off)),family = ziplss(), data = rdata.scaled)
rmodel.int.gam <- gam(list(recruits1 ~ s(BALIVE,k=k) + te(PPT_yr_norm,T_yr_norm,k=k)+ 
                              offset(off),~offset(off)),family = ziplss(), data = rdata.scaled)

res = simulateResiduals(rmodel.clim.gam)
plot(res, quantreg = T) #p = 0

res = simulateResiduals(rmodel.clim.comp.gam)
plot(res, quantreg = T) #p = 0

res = simulateResiduals(rmodel.int.gam)
plot(res, quantreg = T) #p = 0

## Save models and scaling
# growSD is used for building IPM (see BuildIPM.R)
growSD.clim.gam <- sd(resid(gmodel.clim.gam))
growSD.clim.int.gam <- sd(resid(gmodel.clim.int.gam))
growSD.clim.comp.gam <- sd(resid(gmodel.clim.comp.gam))
growSD.int.gam <- sd(resid(gmodel.int.gam))

# specify the predictors in the exported models
gr.predictors <- c("PREVDIA", "T_yr_norm", "PPT_yr_norm", "BALIVE") 

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
save(gmodel.clim.gam,gmodel.clim.int.gam,gmodel.clim.comp.gam,gmodel.int.gam, 
     gr.scaling, growSD.clim.gam,growSD.clim.int.gam,growSD.clim.comp.gam,
     growSD.int.gam,file = "./Code/IPM/GrRescaling_gam.Rdata")

# specify the predictors in the "best" model (or candidate best)
#surv.predictors <- c("PREVDIA", "T_yr_norm", "PPT_yr_norm")
surv.predictors <- c("PREVDIA", "T_yr_norm", "PPT_yr_norm", "BALIVE")
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
surv.scaling.fire = get_scale(survData.scaled, surv.predictors)
# remove scaling information from the dataset so that the model doesnt expect scaled data in predict()
for (i in surv.predictors) {
  attributes(survData.scaled[, i]) = NULL
}

# export model for coefficients and scaling information -------------------
save(smodel.clim.gam,smodel.clim.int.gam,smodel.clim.comp.gam,smodel.int.gam,
     surv.scaling,  file = "./Code/IPM/SurvRescaling_gam.Rdata")

save(smodel.clim.fire.gam,smodel.clim.int.fire.gam,smodel.clim.comp.fire.gam,smodel.int.fire.gam,
     surv.scaling.fire,  file = "./Code/IPM/SurvRescalingFire_gam.Rdata")

# specify the predictors in the models
r.predictors <- c("T_yr_norm", "PPT_yr_norm", "BALIVE") 

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
save(rmodel.clim.gam,rmodel.clim.int.gam,rmodel.clim.comp.gam,rmodel.int.gam,
     r.scaling, file = "./Code/IPM/RecruitRescaling_gam.Rdata")

