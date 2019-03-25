library(ggeffects)
library(ggplot2)
library(coefplot)
library(effects)

# Read and process data
rdata <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)

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

# compare BA.all to BALIVE...they should be the same
# both in sq feet, but BALIVE is on a per acre basis
# BA.all is simply the sum per plot (plot = 672.5 m^2 or 0.16617837 acre)
# expansion factor is 6.01763*BA.all
rdata$BA.all <- rdata$BA.PIED + rdata$BA.notPIED
plot(rdata$BALIVE, 6.01763*rdata$BA.all) # mostly fall on/near one:one, but BALIVE < BA.all

#rdata3 <- subset(rdata, BALIVE > 0) # eliminate data points with no live trees at time 1
# 6635 plots

# standardize covariates
library(dplyr)
rdata.scaled <- rdata %>% mutate_at(scale, .vars = vars(-plot, -lat, -lon, -elev, -PApied,
                                                        -state, -county, -plotID, -CONDID, 
                                                          -measYear, -plotPrev, -PREV_MEASYEAR,
                                                          -CENSUS_INTERVAL, -recruits1, -recruits12,
                                                          -AGB_intra, -BA.PIED, -PIEDadults1,
                                                          -PIEDadults4, -PIEDadults8, -cumDIA.PIED))


### from LISA also

library(DHARMa)
library(glmmTMB)
glmmTMB()

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
#ELS update: AIC = 4979.1
res = simulateResiduals(rmodel_zip)
plot(res, quantreg = T) # ns: p = 0.09159
#ELs update: p = 0.160998
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
#ELS update: AIC = 5250.3
res = simulateResiduals(rmodel_zipcumDIA)
plot(res, quantreg = T) # significant deviation: p = 0.02349
#ELS update: p = 0.01968

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
#ELS update: AIC = 5626.9
res = simulateResiduals(rmodel_zip_PIEDBA)
plot(res, quantreg = T) # significant deviation: p = 0.00022
#ELS update: p = 0.00028

# Proceeding with PIEDadults1 as the offset
# (don't need to choose size predictor - this is independent of adult size)

# Demonstrate the effect of quadratics
rmodel.lin <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_yr_norm + T_yr_norm + 
                      offset(log(CENSUS_INTERVAL)) + 
                      offset(log(PIEDadults1)), # various alternatives for this offset
                    ziformula = ~ 1,
                    data = rdata.scaled, family = "poisson")
rmodel.q <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_yr_norm + T_yr_norm + 
        I(BALIVE^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
        offset(log(CENSUS_INTERVAL)) + 
        offset(log(PIEDadults1)), # various alternatives for this offset
        ziformula = ~ 1,
        data = rdata.scaled, family = "poisson")

mod.comp0<-model.sel(rmodel.lin,rmodel.q)
# strong preference for quadratics (delta AIC = 26.81)

#Compare norms to different time windows
rmodel.1a <- glmmTMB(recruits1 ~ 1
                        + BALIVE + I(BALIVE^2)
                        + PPT_yr_window_25 + I(PPT_yr_window_25^2)
                        + T_yr_window_25 + I(T_yr_window_25^2)
                        + offset(log(CENSUS_INTERVAL))
                        + offset(log(PIEDadults1)), 
                        ziformula = ~ 1,
                        data = rdata.scaled, 
                        family = "poisson")
rmodel.1b <- glmmTMB(recruits1 ~ 1
                        + BALIVE + I(BALIVE^2)
                        + PPT_yr_window_20 + I(PPT_yr_window_20^2)
                        + T_yr_window_20 + I(T_yr_window_20^2)
                        + offset(log(CENSUS_INTERVAL))
                        + offset(log(PIEDadults1)), 
                        ziformula = ~ 1,
                        data = rdata.scaled, 
                        family = "poisson")
# 15 yrs (Poisson)
rmodel.1c <- glmmTMB(recruits1 ~ 1
                        + BALIVE + I(BALIVE^2)
                        + PPT_yr_window_15 + I(PPT_yr_window_15^2)
                        + T_yr_window_15 + I(T_yr_window_15^2)
                        + offset(log(CENSUS_INTERVAL))
                        + offset(log(PIEDadults1)), 
                        ziformula = ~ 1,
                        data = rdata.scaled, 
                        family = "poisson")
mod.comp1 <- model.sel(rmodel.q,rmodel.1a, rmodel.1b,rmodel.1c)
#No difference, will proceed with norms for consistency with other vital rate models

# compare annual vs. 3 vs. 4 seasons...likes 3-season better 
rmodel.2a <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + 
                      I(BALIVE^2) + 
                      I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                      I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) +
                      offset(log(CENSUS_INTERVAL)) + 
                      offset(log(PIEDadults1)), # various alternatives for this offset
                    ziformula = ~ 1,
                    data = rdata.scaled, family = "poisson")
rmodel.2b <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + + PPT_pf_norm + PPT_fs_norm + PPT_m_norm + 
                       T_c_norm + T_pf_norm + T_fs_norm + T_m_norm + 
                       I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_pf_norm^2) + I(PPT_fs_norm^2) + I(PPT_m_norm^2) + 
                      I(T_c_norm^2) + I(T_pf_norm^2) + I(T_fs_norm^2) + (T_m_norm^2) + 
                      offset(log(CENSUS_INTERVAL)) + 
                      offset(log(PIEDadults1)), # various alternatives for this offset
                    ziformula = ~ 1,
                    data = rdata.scaled, family = "poisson")

mod.comp2 <- model.sel(rmodel.q, rmodel.2a, rmodel.2b)

# compare normals vs. anomalies...census interval is best
rmodel.3a <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_yr_anom + T_yr_anom + I(BALIVE^2) + 
                       I(PPT_yr_anom^2) + I(T_yr_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3b <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_anom + PPT_wd_anom + PPT_m_anom + 
                       T_c_anom + T_wd_anom + T_m_anom + I(BALIVE^2) + 
                       I(PPT_c_anom^2) + I(PPT_wd_anom^2) + I(PPT_m_anom^2) + 
                       I(T_c_anom^2) + I(T_wd_anom^2) + (T_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3c <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_yr_anom + T_yr_anom + PPT_yr_norm + T_yr_norm + 
                       I(BALIVE^2) + 
                       I(PPT_yr_anom^2) + I(T_yr_anom^2) + I(PPT_yr_norm^2) + I(T_yr_norm^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3d <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPT_c_anom + PPT_wd_anom + PPT_m_anom + 
                       T_c_anom + T_wd_anom + T_m_anom + 
                       I(PPT_c_anom^2) + I(PPT_wd_anom^2) + I(PPT_m_anom^2) + 
                       I(T_c_anom^2) + I(T_wd_anom^2) + (T_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3e <- glmmTMB(recruits1 ~ 1 + BALIVE + PPTex_yr_anom + Tex_yr_anom + I(BALIVE^2) + 
                       I(PPTex_yr_anom^2) + I(Tex_yr_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3f <- glmmTMB(recruits1 ~ 1 + BALIVE + PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + I(BALIVE^2) + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3g <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_yr_norm + T_yr_norm + I(BALIVE^2) + 
                       I(PPT_yr_norm^2) + I(T_yr_norm^2) + PPTex_yr_anom + Tex_yr_anom + 
                       I(PPTex_yr_anom^2) + I(Tex_yr_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.3h <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")

mod.comp3 <- model.sel(rmodel.2a, rmodel.3a, rmodel.3b, rmodel.3c, rmodel.3d, rmodel.3e, rmodel.3f,
                       rmodel.3g, rmodel.3h)
# rmodel.3h is best among these 

# add drought anomalies
rmodel.4a <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + 
                       PPT_drought + Tmean_drought + I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       I(PPT_drought^2) + I(Tmean_drought^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.4b <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + 
                       PPT_pf_dr + PPT_c_dr + PPT_fs_dr + PPT_m_dr + Tmean_drought + I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       I(PPT_pf_dr_anom^2) + I(PPT_c_dr^2) + I(PPT_fs_dr^2) + I(PPT_m_dr^2) + I(Tmean_drought^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.4c <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + 
                       PPT_dr_anom + T_dr_anom + I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       I(PPT_dr_anom^2) + I(T_dr_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.4d <- glmmTMB(recruits1 ~ 1 + BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                       T_c_norm + T_wd_norm + T_m_norm + 
                       PPT_pf_dr_anom + PPT_c_dr_anom + PPT_fs_dr_anom + PPT_m_dr_anom + T_dr_anom + 
                       I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       I(PPT_pf_dr_anom^2) + I(PPT_c_dr_anom^2) + I(PPT_fs_dr_anom^2) + I(PPT_m_dr_anom^2) + I(T_dr_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")

mod.comp4<-model.sel(rmodel.3h,rmodel.4a,rmodel.4b,rmodel.4c,rmodel.4d)
# rmodel.3h and 4c are the best - proceeding with 3h because it is simpler

# add 2-way interactions, excluding quadratics
rmodel.5a <- glmmTMB(recruits1 ~ 1 + (BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                                        T_c_norm + T_wd_norm + T_m_norm + 
                                        PPT_dr_anom + T_dr_anom)^2 + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.5b <- glmmTMB(recruits1 ~ 1 + (BALIVE + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                                        T_c_norm + T_wd_norm + T_m_norm)^2 + I(BALIVE^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")

mod.comp5 <- model.sel(rmodel.3h, rmodel.5a, rmodel.5b)
# rmodel.5b is best

# Try different basal area predictors
rmodel.6a <- glmmTMB(recruits1 ~ 1 + (BA.PIPO + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                                        T_c_norm + T_wd_norm + T_m_norm)^2 + I(BA.PIPO^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.6b <- glmmTMB(recruits1 ~ 1 + (BA.notPIED + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                                        T_c_norm + T_wd_norm + T_m_norm)^2 + I(BA.notPIED^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
rmodel.6c <- glmmTMB(recruits1 ~ 1 + (BA.juniper + PPT_c_norm + PPT_wd_norm + PPT_m_norm + 
                                        T_c_norm + T_wd_norm + T_m_norm)^2 + I(BA.juniper^2) + 
                       I(PPT_c_norm^2) + I(PPT_wd_norm^2) + I(PPT_m_norm^2) + 
                       I(T_c_norm^2) + I(T_wd_norm^2) + (T_m_norm^2) + 
                       PPTex_c_anom + PPTex_wd_anom + PPTex_m_anom + 
                       Tex_c_anom + Tex_wd_anom + Tex_m_anom + 
                       I(PPTex_c_anom^2) + I(PPTex_wd_anom^2) + I(PPTex_m_anom^2) + 
                       I(Tex_c_anom^2) + I(Tex_wd_anom^2) + (Tex_m_anom^2) + 
                       offset(log(CENSUS_INTERVAL)) + 
                       offset(log(PIEDadults1)), # various alternatives for this offset
                     ziformula = ~ 1,
                     data = rdata.scaled, family = "poisson")
mod.comp6<-model.sel(rmodel.5b,rmodel.6a,rmodel.6b,rmodel.6c)
# 5b is by far the best

res = simulateResiduals(rmodel.5b)
plot(res, quantreg = T) # ns: p = 0.44424

### Models to export
rmodel.clim<-glmmTMB(recruits1 ~ 1
                       + PPT_yr_norm + T_yr_norm 
                       + I(PPT_yr_norm^2) + I(T_yr_norm^2) 
                       + offset(log(CENSUS_INTERVAL))
                       + offset(log(PIEDadults1)), # various alternatives for this offset
                       ziformula = ~ 1,
                       data = rdata.scaled, 
                       family = "poisson")
rmodel.clim.comp<-glmmTMB(recruits1 ~ 1
                       + BALIVE + PPT_yr_norm + T_yr_norm 
                       + offset(log(CENSUS_INTERVAL))
                       + offset(log(PIEDadults1)), # various alternatives for this offset
                       ziformula = ~ 1,
                       data = rdata.scaled, 
                       family = "poisson")
rmodel.int<-glmmTMB(recruits1 ~ 1
               + (BALIVE + PPT_yr_norm + T_yr_norm)^2 
               + offset(log(CENSUS_INTERVAL))
               + offset(log(PIEDadults1)), # various alternatives for this offset
               ziformula = ~ 1,
               data = rdata.scaled, 
               family = "poisson")
rmodel.best<-rmodel.5b

### dealing with std'ized covariates

# specify the predictors in the "best" model (or candidate best)
r.predictors <- c("T_yr_norm", "PPT_yr_norm", "BALIVE", 
                  "PPT_wd_norm", "PPT_c_norm", "PPT_m_norm", "T_wd_norm", "T_c_norm", "T_m_norm",
                  "PPTex_wd_anom", "PPTex_c_anom", "PPTex_m_anom", "Tex_wd_anom", "Tex_c_anom", "Tex_m_anom") # rmodel_zipoiss
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
save(rmodel.clim,rmodel.clim.comp,rmodel.int,rmodel.best, r.scaling, file = "./Code/IPM/RecruitRescaling.Rdata")
