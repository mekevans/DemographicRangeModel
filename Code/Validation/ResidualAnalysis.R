########## Residual analysis for PIED############

## Created by Emily Schultz
## Created on 29 Aug 2019

library(mgcv)
library(scam)
library(DHARMa) # use DHARMa to check residuals
library(tidyverse)
library(lattice)
library(cowplot)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(penalized)
library(scales)
library(akima)
library(ggalt)

invlogit<-function(x){exp(x)/(1+exp(x))}

mytheme2<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.text=element_text(size=11),legend.title=element_text(size=11),
                legend.key = element_rect(fill = "white"),axis.text=element_text(size=11),
                axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
                axis.line.x = element_line(color="black", size = 0.3),
                axis.line.y = element_line(color="black", size = 0.3))

PApied<-raster("./Processed/Validation/presenceAbsenceRaster2.tif")

PRISM.norm.path <-  "./ClimateData/PRISM/Normals/"
ppt_yr_raster <- raster(paste0(PRISM.norm.path, "PPT_year.tif"))
t_yr_raster <- raster(paste0(PRISM.norm.path, "T_year.tif"))
ba_raster <- raster("./BA/balive_RF.tif")
ppt_yr_raster <- aggregate(ppt_yr_raster, 4)
t_yr_raster <- aggregate(t_yr_raster, 4)
ba_raster <- aggregate(ba_raster, 4)

lambda_c<-raster("./Output/tifs/PIED.clim_lambda_gam.tif")
#lambda_ci<-raster("./Output/tifs/PIED.climint_lambda_gam.tif")
lambda_cc<-raster("./Output/tifs/PIED.climcomp_lambda_gam.tif")
lambda_i<-raster("./Output/tifs/PIED.int_lambda_gam.tif")

FIA_lambda_c<-read.csv("./Output/FIA_lambda_gam_c.csv")
FIA_lambda_cc<-read.csv("./Output/FIA_lambda_gam_cc.csv")
FIA_lambda_i<-read.csv("./Output/FIA_lambda_gam_i.csv")

FIA_lambda<-merge(FIA_lambda_cc,FIA_lambda_i)

FIA_lambda<-as.data.frame(PApied,xy=TRUE)
names(FIA_lambda)<-c("lon","lat","PApied")
FIA_lambda$lambda_c<-values(lambda_c)
#FIA_lambda$lambda_ci<-values(lambda_ci)
FIA_lambda$lambda_cc<-values(lambda_cc)
FIA_lambda$lambda_i<-values(lambda_i)
FIA_lambda$BALIVE<-values(ba_raster)
FIA_lambda$PPT_yr<-values(ppt_yr_raster)
FIA_lambda$T_yr<-values(t_yr_raster)

FIA_lambda<-na.omit(FIA_lambda)

FIA_pied_pres<-subset(FIA_lambda,PApied==1)

min_ba<-min(FIA_pied_pres$BALIVE, na.rm=T)
max_ba<-max(FIA_pied_pres$BALIVE, na.rm=T)
min_ppt<-min(FIA_pied_pres$PPT_yr, na.rm=T)
max_ppt<-max(FIA_pied_pres$PPT_yr, na.rm=T)
min_t<-min(FIA_pied_pres$T_yr, na.rm=T)
max_t<-max(FIA_pied_pres$T_yr, na.rm=T)

FIA_lambda_noex<-subset(FIA_lambda,FIA_lambda$PPT_yr<max_ppt & FIA_lambda$PPT_yr>min_ppt &
                          FIA_lambda$T_yr<max_t & FIA_lambda$T_yr>min_t &
                          FIA_lambda$BALIVE<max_ba & FIA_lambda$BALIVE>min_ba)
##Residual Analysis
k=5
pa_c <- gam(PApied~s(lambda_c,k=k),family=binomial(),
          data=FIA_lambda_noex)
pa_ci <- gam(PApied~s(lambda_ci,k=k),family=binomial(),
            data=FIA_lambda)
pa_cc <- gam(PApied~s(lambda_cc,k=k),family=binomial(),
            data=FIA_lambda_noex)
pa_i <- gam(PApied~s(lambda_i),family=binomial(),
             data=FIA_lambda_noex)

pa_c <- glm(PApied~lambda_c,family=binomial(),
            data=FIA_lambda)
#pa_ci <- glm(PApied~lambda_ci,family=binomial(),
#            data=FIA_lambda)
pa_cc <- glm(PApied~lambda_cc,family=binomial(),
            data=FIA_lambda)
pa_i <- glm(PApied~lambda_i,family=binomial(),
            data=FIA_lambda)
#,bs="mpi",m=2

#Plot models
ncuts=50
chopsize_lam_c<-cut(FIA_lambda$lambda_c,ncuts)
#chopsize_lam_ci<-cut(FIA_lambda$lambda_ci,ncuts)
chopsize_lam_cc<-cut(FIA_lambda$lambda_cc,ncuts)
chopsize_lam_i<-cut(FIA_lambda$lambda_i,ncuts)

count_binned_lam_c<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_c),length))
lam_c_binned<-as.vector(sapply(split(FIA_lambda$lambda_c,chopsize_lam_c),mean,na.rm=T))
pres_c_binned<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_c),mean,na.rm=T))

#count_binned_lam_ci<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_ci),length))
#lam_ci_binned<-as.vector(sapply(split(FIA_lambda$lambda_ci,chopsize_lam_ci),mean,na.rm=T))
#pres_ci_binned<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_ci),mean,na.rm=T))

count_binned_lam_cc<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_cc),length))
lam_cc_binned<-as.vector(sapply(split(FIA_lambda$lambda_cc,chopsize_lam_cc),mean,na.rm=T))
pres_cc_binned<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_cc),mean,na.rm=T))

count_binned_lam_i<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_i),length))
lam_i_binned<-as.vector(sapply(split(FIA_lambda$lambda_i,chopsize_lam_i),mean,na.rm=T))
pres_i_binned<-as.vector(sapply(split(FIA_lambda$PApied,chopsize_lam_i),mean,na.rm=T))

pres_binned<-data.frame(count_lam=c(count_binned_lam_c,count_binned_lam_cc,count_binned_lam_i),
                        lam=c(lam_c_binned,lam_cc_binned,lam_i_binned),
                        pres=c(pres_c_binned,pres_cc_binned,pres_i_binned),
                        pred=c(invlogit(predict(pa_c,newdata=data.frame(lambda_c=lam_c_binned))),
                               invlogit(predict(pa_cc,newdata=data.frame(lambda_cc=lam_cc_binned))),
                               invlogit(predict(pa_i,newdata=data.frame(lambda_i=lam_i_binned)))),
                        model=c(rep("c",ncuts),rep("cc",ncuts),rep("i",ncuts)))


#No extrapolation
#Plot models
ncuts=50
chopsize_lam_c<-cut(FIA_lambda_noex$lambda_c,ncuts)
#chopsize_lam_ci<-cut(FIA_lambda_noex$lambda_ci,ncuts)
chopsize_lam_cc<-cut(FIA_lambda_noex$lambda_cc,ncuts)
chopsize_lam_i<-cut(FIA_lambda_noex$lambda_i,ncuts)

count_binned_lam_c<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_c),length))
lam_c_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_c,chopsize_lam_c),mean,na.rm=T))
pres_c_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_c),mean,na.rm=T))

#count_binned_lam_ci<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_ci),length))
#lam_ci_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_ci,chopsize_lam_ci),mean,na.rm=T))
#pres_ci_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_ci),mean,na.rm=T))

count_binned_lam_cc<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_cc),length))
lam_cc_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_cc,chopsize_lam_cc),mean,na.rm=T))
pres_cc_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_cc),mean,na.rm=T))

count_binned_lam_i<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_i),length))
lam_i_binned<-as.vector(sapply(split(FIA_lambda_noex$lambda_i,chopsize_lam_i),mean,na.rm=T))
pres_i_binned<-as.vector(sapply(split(FIA_lambda_noex$PApied,chopsize_lam_i),mean,na.rm=T))

pres_binned<-data.frame(count_lam=c(count_binned_lam_c,count_binned_lam_cc,count_binned_lam_i),
                        lam=c(lam_c_binned,lam_cc_binned,lam_i_binned),
                        pres=c(pres_c_binned,pres_cc_binned,pres_i_binned),
                        pred=c(invlogit(predict(pa_c,newdata=data.frame(lambda_c=lam_c_binned))),
                               invlogit(predict(pa_cc,newdata=data.frame(lambda_cc=lam_cc_binned))),
                               invlogit(predict(pa_i,newdata=data.frame(lambda_i=lam_i_binned)))),
                        model=c(rep("c",ncuts),rep("cc",ncuts),rep("i",ncuts)))


## Residuals
res_c = simulateResiduals(pa_c)
#res_ci = simulateResiduals(pa_ci)
res_cc = simulateResiduals(pa_cc)
res_i = simulateResiduals(pa_i)

FIA_lambda$resid_c<-res_c$scaledResiduals
#FIA_lambda$resid_ci<-res_ci$scaledResiduals
FIA_lambda$resid_cc<-res_cc$scaledResiduals
FIA_lambda$resid_i<-res_i$scaledResiduals
FIA_lambda$pred_c<-invlogit(pa_c$coefficients[1]+pa_c$coefficients[2]*FIA_lambda$lambda_c)
#FIA_lambda$pred_ci<-invlogit(pa_ci$coefficients[1]+pa_ci$coefficients[2]*FIA_lambda$lambda_ci)
FIA_lambda$pred_cc<-invlogit(pa_cc$coefficients[1]+pa_cc$coefficients[2]*FIA_lambda$lambda_cc)
FIA_lambda$pred_i<-invlogit(pa_i$coefficients[1]+pa_i$coefficients[2]*FIA_lambda$lambda_i)

FIA_lambda$resid_c2<-FIA_lambda$pred_c-FIA_lambda$PApied
#FIA_lambda$resid_ci2<-FIA_lambda$pred_ci-FIA_lambda$PApied
FIA_lambda$resid_cc2<-FIA_lambda$pred_cc-FIA_lambda$PApied
FIA_lambda$resid_i2<-FIA_lambda$pred_i-FIA_lambda$PApied

testSpatialAutocorrelation(simulationOutput = res_c, x = FIA_lambda$lon, y= FIA_lambda$lat)
testSpatialAutocorrelation(simulationOutput = res_c)

resid_group<-recalculateResiduals(res_c , group = xxx)

ncuts=30
#chopsize_elev<-cut(FIA_lambda$elev,ncuts)
chopsize_ba<-cut(FIA_lambda$BALIVE,ncuts)
chopsize_ppt<-cut(FIA_lambda$PPT_yr,ncuts)
chopsize_t<-cut(FIA_lambda$T_yr,ncuts)

#count_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_elev),length))
#elev_binned<-as.vector(sapply(split(FIA_lambda$elev,chopsize_elev),mean,na.rm=T))
#resid_c_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_elev),mean,na.rm=T))
#resid_c2_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_c2,chopsize_elev),mean,na.rm=T))
#resid_ci_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_ci,chopsize_elev),mean,na.rm=T))
#resid_ci2_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_ci2,chopsize_elev),mean,na.rm=T))
#resid_cc_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_cc,chopsize_elev),mean,na.rm=T))
#resid_cc2_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_cc2,chopsize_elev),mean,na.rm=T))
#resid_i_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_i,chopsize_elev),mean,na.rm=T))
#resid_i2_binned_elev<-as.vector(sapply(split(FIA_lambda$resid_i2,chopsize_elev),mean,na.rm=T))

count_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_ba),length))
ba_binned<-as.vector(sapply(split(FIA_lambda$BALIVE,chopsize_ba),mean,na.rm=T))
resid_c_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_ba),mean,na.rm=T))
resid_c2_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_c2,chopsize_ba),mean,na.rm=T))
#resid_ci_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_ci,chopsize_ba),mean,na.rm=T))
#resid_ci2_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_ci2,chopsize_ba),mean,na.rm=T))
resid_cc_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_cc,chopsize_ba),mean,na.rm=T))
resid_cc2_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_cc2,chopsize_ba),mean,na.rm=T))
resid_i_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_i,chopsize_ba),mean,na.rm=T))
resid_i2_binned_ba<-as.vector(sapply(split(FIA_lambda$resid_i2,chopsize_ba),mean,na.rm=T))

count_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_ppt),length))
ppt_binned<-as.vector(sapply(split(FIA_lambda$PPT_yr,chopsize_ppt),mean,na.rm=T))
resid_c_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_ppt),mean,na.rm=T))
resid_c2_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_c2,chopsize_ppt),mean,na.rm=T))
#resid_ci_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_ci,chopsize_ppt),mean,na.rm=T))
#resid_ci2_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_ci2,chopsize_ppt),mean,na.rm=T))
resid_cc_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_cc,chopsize_ppt),mean,na.rm=T))
resid_cc2_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_cc2,chopsize_ppt),mean,na.rm=T))
resid_i_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_i,chopsize_ppt),mean,na.rm=T))
resid_i2_binned_ppt<-as.vector(sapply(split(FIA_lambda$resid_i2,chopsize_ppt),mean,na.rm=T))

count_binned_t<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_t),length))
t_binned<-as.vector(sapply(split(FIA_lambda$T_yr,chopsize_t),mean,na.rm=T))
resid_c_binned_t<-as.vector(sapply(split(FIA_lambda$resid_c,chopsize_t),mean,na.rm=T))
resid_c2_binned_t<-as.vector(sapply(split(FIA_lambda$resid_c2,chopsize_t),mean,na.rm=T))
#resid_ci_binned_t<-as.vector(sapply(split(FIA_lambda$resid_ci,chopsize_t),mean,na.rm=T))
#resid_ci2_binned_t<-as.vector(sapply(split(FIA_lambda$resid_ci2,chopsize_t),mean,na.rm=T))
resid_cc_binned_t<-as.vector(sapply(split(FIA_lambda$resid_cc,chopsize_t),mean,na.rm=T))
resid_cc2_binned_t<-as.vector(sapply(split(FIA_lambda$resid_cc2,chopsize_t),mean,na.rm=T))
resid_i_binned_t<-as.vector(sapply(split(FIA_lambda$resid_i,chopsize_t),mean,na.rm=T))
resid_i2_binned_t<-as.vector(sapply(split(FIA_lambda$resid_i2,chopsize_t),mean,na.rm=T))

res_binned<-data.frame(#count_elev=count_binned_elev,elev=elev_binned,
                       #resid_c_elev=resid_c_binned_elev,resid_ci_elev=resid_ci_binned_elev,
                       #resid_cc_elev=resid_cc_binned_elev,resid_i_elev=resid_i_binned_elev,
                       count_ba=count_binned_ba,ba=ba_binned,
                       resid_c_ba=resid_c_binned_ba,#resid_ci_ba=resid_ci_binned_ba,
                       resid_cc_ba=resid_cc_binned_ba,resid_i_ba=resid_i_binned_ba,
                       count_ppt=count_binned_ppt,ppt=ppt_binned,
                       resid_c_ppt=resid_c_binned_ppt,#resid_ci_ppt=resid_ci_binned_ppt,
                       resid_cc_ppt=resid_cc_binned_ppt,resid_i_ppt=resid_i_binned_ppt,
                       count_t=count_binned_t,t=t_binned,
                       resid_c_t=resid_c_binned_t,#resid_ci_t=resid_ci_binned_t,
                       resid_cc_t=resid_cc_binned_t,resid_i_t=resid_i_binned_t)

res2_binned<-data.frame(#count_elev=rep(count_binned_elev,4),elev=rep(elev_binned,4),
                       #resid_elev=c(resid_c2_binned_elev,resid_ci2_binned_elev,
                        #            resid_cc2_binned_elev,resid_i2_binned_elev),
                       count_ba=rep(count_binned_ba,3),ba=rep(ba_binned,3),
                       resid_ba=c(resid_c2_binned_ba,#resid_ci2_binned_ba,
                                  resid_cc2_binned_ba,resid_i2_binned_ba),
                       count_ppt=rep(count_binned_ppt,3),ppt=rep(ppt_binned,3),
                       resid_ppt=c(resid_c2_binned_ppt,#resid_ci2_binned_ppt,
                                   resid_cc2_binned_ppt,resid_i2_binned_ppt),
                       count_t=rep(count_binned_t,3),t=rep(t_binned,3),
                       resid_t=c(resid_c2_binned_t,#resid_ci2_binned_t,
                                 resid_cc2_binned_t,resid_i2_binned_t),
                       model=c(rep("c",ncuts),rep("cc",ncuts),rep("i",ncuts)))

# Variance-covariance
FIA_pres<-subset(FIA_lambda,PApied==1)
FIA_abs<-subset(FIA_lambda,PApied==0)

splom(FIA_pres[c(7,8,9)])
splom(FIA_abs[c(7,8,9)])

pres_cov<-cov(FIA_pres[c(7,8,9)])
abs_cov<-cov(FIA_abs[c(7,8,9)])
pres_cor<-cor(FIA_pres[c(7,8,9)])
abs_cor<-cor(FIA_abs[c(7,8,9)])

##Lambda histograms
FIA_lambda$PApied_f<-as.factor(FIA_lambda$PApied)
l_means<-FIA_lambda %>%
  group_by(PApied_f) %>%
  summarise(lambda_c=mean(lambda_c,na.rm=T),lambda_cc=mean(lambda_cc,na.rm=T),
            lambda_i=mean(lambda_i,na.rm=T))
l_meds<-FIA_lambda %>%
  group_by(PApied_f) %>%
  summarise(lambda_c=median(lambda_c,na.rm=T),lambda_cc=median(lambda_cc,na.rm=T),
            lambda_i=median(lambda_i,na.rm=T))

# Niche space plots
FIA_lambda<-na.omit(FIA_lambda)
interpdf_c <-interp2xyz(interp(x=FIA_lambda$PPT_yr, y=FIA_lambda$T_yr, z=FIA_lambda$lambda_c, duplicate="mean"), data.frame=TRUE)
interpdf_cc <-interp2xyz(interp(x=FIA_lambda$PPT_yr, y=FIA_lambda$T_yr, z=FIA_lambda$lambda_cc, duplicate="mean"), data.frame=TRUE)
interpdf_i <-interp2xyz(interp(x=FIA_lambda$PPT_yr, y=FIA_lambda$T_yr, z=FIA_lambda$lambda_i, duplicate="mean"), data.frame=TRUE)

interpdf_c<-na.omit(interpdf_c)
interpdf_cc<-na.omit(interpdf_cc)
interpdf_i<-na.omit(interpdf_i)

mean_lam_c<-mean(FIA_pied_pres$lambda_c)
mean_lam_cc<-mean(FIA_pied_pres$lambda_cc)
mean_lam_i<-mean(FIA_pied_pres$lambda_i)

interpdf_c$z_diff<-interpdf_c$z-mean_lam_c
interpdf_cc$z_diff<-interpdf_cc$z-mean_lam_cc
interpdf_i$z_diff<-(interpdf_i$z-mean_lam_i)

q_neg_c<-quantile(interpdf_c$z_diff[which(interpdf_c$z_diff<0)])
q_pos_c<-quantile(interpdf_c$z_diff[which(interpdf_c$z_diff>0)])
q_c<-c(q_neg_c[2:4],q_pos_c[2:4])

q_neg_cc<-quantile(interpdf_cc$z_diff[which(interpdf_cc$z_diff<0)])
q_pos_cc<-quantile(interpdf_cc$z_diff[which(interpdf_cc$z_diff>0)])
q_cc<-c(q_neg_cc[2:4],q_pos_cc[2:4])

q_neg_i<-quantile(interpdf_i$z_diff[which(interpdf_i$z_diff<0)])
q_pos_i<-quantile(interpdf_i$z_diff[which(interpdf_i$z_diff>0)])
q_i<-c(q_neg_i[2:4],q_pos_i[2:4])

## Save output
save(FIA_lambda, FIA_pied_pres, pres_binned, res_binned, res2_binned, l_means, l_meds, 
     interpdf_c, interpdf_cc, interpdf_i, q_c, q_cc, q_i, 
     file="./Output/residual_gam.rda")

save(FIA_pres, FIA_abs, pres_cov, pres_cor, abs_cov, abs_cor, file="./Output/covariance.rda")
