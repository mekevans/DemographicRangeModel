library(raster)
#library(rgdal)
#library(rgeos)
#library(dplyr)
library(glmmTMB)
library(mgcv)
library(tidyverse)
#library(brms)

# Load vital rate and IPM functions
source("./Code/IPM/BuildIPM.R")

PRISM.norm.path <-  "./ClimateData/PRISM/Normals/"

# load gam models and scaling --------------------------------------------------

load("./Code/IPM/GrRescaling_gam.Rdata")
load("./Code/IPM/SurvRescaling_gam.Rdata")
load("./Code/IPM/RecruitRescaling_gam.Rdata")
load("./Code/IPM/recrstats.rda")
#load("./Code/IPM/RecruitRescaling8_gam.Rdata")

# Load FIA survival, growth data
FIA <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)
# when running IPM without trees killed by fire, technically should filter out those trees
# prolly makes ~0 difference
FIA <- FIA[!(FIA$DSTRBCD1 %in% c(30, 31, 32, 80)), ]

# Set IPM parameters
min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500

# Set aggregation factor
#aggr <- 20

# Load climate layers
# Load climate layers
# these created using the script "current.R"
ppt_yr_raster <- raster(paste0(PRISM.norm.path, "PPT_year.tif"))
t_yr_raster <- raster(paste0(PRISM.norm.path, "T_year.tif"))
# stand-level basal area raster
ba_raster <- raster("./BA/balive_RF.tif")

ppt_yr_raster <- resample(ppt_yr_raster, ba_raster)
t_yr_raster <- resample(t_yr_raster, ba_raster)

ppt_yr_raster <- aggregate(ppt_yr_raster, 4)
t_yr_raster <- aggregate(t_yr_raster, 4)
ba_raster <- aggregate(ba_raster, 4)

#Load rasters of unperturbed lambda values
lambda_c<-raster("./Output/tifs/PIED.clim_lambda_gam.tif")
#lambda_ci<-raster("./Output/tifs/PIED.climint_lambda_gam.tif")
lambda_cc<-raster("./Output/tifs/PIED.climcomp_lambda_gam.tif")
lambda_i<-raster("./Output/tifs/PIED.int_lambda_gam.tif")

#LTRE

#Calculate mean climate conditions where PIED is present
FIA2 <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)

mean_ppt_pied <- mean(subset(FIA2,PApied==1)$PPT_yr_norm)
mean_t_pied <- mean(subset(FIA2,PApied==1)$T_yr_norm)
mean_ba_pied <- mean(subset(FIA2,PApied==1)$BALIVE)

save(mean_ppt_pied,mean_t_pied,mean_ba_pied, file="./Output/clim_means.rda")

# Calculate lambda under mean climate conditions
pred_data<-data.frame(PPT_yr=mean_ppt_pied,T_yr=mean_t_pied,BALIVE=mean_ba_pied)
K_mean <- ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
          rmodel=rmodel.int.gam, gSD=growSD.int.gam,
          data=pred_data,
          s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
lambda_mean <- Re(eigen(K_mean)$values[1])

# Calculate differences between lambda in each pixel and mean lambda
lambda_diff_i<-lambda_i-lambda_mean

lambda_ltre_g <- ppt_yr_raster
lambda_ltre_g <- setValues(lambda_ltre_g, 0)
lambda_ltre_s <- lambda_ltre_g 
lambda_ltre_r <- lambda_ltre_g 

for (i in 1:nrow(ppt_yr_raster)) { #
  for (j in 1:ncol(ppt_yr_raster)){ #
    
  # Extract climate for cell
  #  pred_data <- data.frame(PPT_yr=mean(as.numeric(ppt_yr_raster[i,j]),mean_ppt_pied),
  #                          T_yr=mean(as.numeric(t_yr_raster[i,j]),mean_t_pied),
  #                          BALIVE=mean(as.numeric(ba_raster[i,j]),mean_ba_pied))
  
  #Number of perturbation steps
  inc=1
  
  # Climate at mean
  pred_data <- data.frame(PPT_yr=mean_ppt_pied,
                          T_yr=mean_t_pied,
                          BALIVE=mean_ba_pied)
  
  # Climate at pixel
  pred_data_new <- data.frame(PPT_yr=as.numeric(ppt_yr_raster[i,j]),
                              T_yr=as.numeric(t_yr_raster[i,j]),
                              BALIVE=as.numeric(ba_raster[i,j]))
  
  if (is.na(pred_data_new$PPT_yr) | is.na(pred_data_new$T_yr) | is.na(pred_data_new$BALIVE)) {
    lambda_ltre_g[i,j] <- NA
    lambda_ltre_s[i,j] <- NA
    lambda_ltre_r[i,j] <- NA
    print("Missing predictor... Skipping!")
    next
  }
  
  # Calculate change in predictors
  d_pred_data<-pred_data_new-pred_data
  
  #Divide change in each parameter into amount to perturb each time
  d_pred_data_inc<-d_pred_data/inc
  #perturb<-c(0,0,0)
  
  lambda_old <- lambda_mean
  pred_data_inc<-list(pred_data,pred_data,pred_data)
  for (p in 1:inc){# Calculate lambda
    dlambda<-rep(NA,3)
    for (v in 1:3) {
      pred_data_inc[[v]]<-pred_data_inc[[v]]+d_pred_data_inc
      K<-ipm_ltre(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                 rmodel=rmodel.int.gam, gSD=growSD.int.gam, data_g=pred_data_inc[[1]],data_s=pred_data_inc[[2]],data_r=pred_data_inc[[3]], 
                 s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
      lambda<-Re(eigen(K)$values[1])
      dlambda[v]<-lambda-lambda_old
      lambda_old<-lambda
    }
        
    lambda_ltre_g[i,j] <- lambda_ltre_g[i,j]+dlambda[1] 
    lambda_ltre_s[i,j] <- lambda_ltre_s[i,j]+dlambda[2] 
    lambda_ltre_r[i,j] <- lambda_ltre_r[i,j]+dlambda[3] 
  }
  print(j)
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(lambda_ltre_g, main = "Lambda"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

lambda_ltre_vital<-sum(lambda_ltre_g,lambda_ltre_s,lambda_ltre_r)

writeRaster(lambda_ltre_vital, "./Output/tifs/PIED.int_ltre_vital.tif", overwrite = T)
writeRaster(lambda_ltre_g, "./Output/tifs/PIED.int_ltre_g.tif", overwrite = T)
writeRaster(lambda_ltre_s, "./Output/tifs/PIED.int_ltre_s.tif", overwrite = T)
writeRaster(lambda_ltre_r, "./Output/tifs/PIED.int_ltre_r.tif", overwrite = T)

# LTRE map environmental predictors
lambda_ltre_p <- ppt_yr_raster
lambda_ltre_p <- setValues(lambda_ltre_p, 0)
lambda_ltre_t <- lambda_ltre_p
lambda_ltre_b <- lambda_ltre_p 

for (i in 1:nrow(ppt_yr_raster)) { #(i in 6:10){ #
  for (j in 1:ncol(ppt_yr_raster)){ #(j in 1:1){ #
    
    # Extract climate for cell
    #  pred_data <- data.frame(PPT_yr=mean(as.numeric(ppt_yr_raster[i,j]),mean_ppt_pied),
    #                          T_yr=mean(as.numeric(t_yr_raster[i,j]),mean_t_pied),
    #                          BALIVE=mean(as.numeric(ba_raster[i,j]),mean_ba_pied))
    
    #Number of perturbation steps
    inc=1
    
    # Climate at mean
    pred_data <- data.frame(PPT_yr=mean_ppt_pied,
                            T_yr=mean_t_pied,
                            BALIVE=mean_ba_pied)
    
    # Climate at pixel
    pred_data_new <- data.frame(PPT_yr=as.numeric(ppt_yr_raster[i,j]),
                                T_yr=as.numeric(t_yr_raster[i,j]),
                                BALIVE=as.numeric(ba_raster[i,j]))
    
    if (is.na(pred_data_new$PPT_yr) | is.na(pred_data_new$T_yr) | is.na(pred_data_new$BALIVE) | is.nan(pred_data_new$BALIVE)) {
      lambda_ltre_p[i,j] <- NA
      lambda_ltre_t[i,j] <- NA
      lambda_ltre_b[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    
    # Calculate change in predictors and divide change in each parameter into amount to perturb each time
    dclim <- c((pred_data_new$PPT_yr-pred_data$PPT_yr),(pred_data_new$T_yr-pred_data$T_yr),
               (pred_data_new$BALIVE-pred_data$BALIVE))
    perturb_inc<-dclim/inc
    perturb<-c(0,0,0)
    
    lambda_old <- lambda_mean
    
    for (p in 1:inc){# Calculate lambda
      dlambda<-rep(NA,3)
      for (v in 1:3) {
        perturb[v]<-perturb[v]+perturb_inc[v]
        pred_data<-pred_data+perturb
        K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                   rmodel=rmodel.int.gam, gSD=growSD.int.gam, data=pred_data, 
                   s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
        lambda<-Re(eigen(K)$values[1])
        dlambda[v]<-lambda-lambda_old
        lambda_old<-lambda
        perturb<-c(0,0,0)
      }
      
      lambda_ltre_p[i,j] <- lambda_ltre_p[i,j]+dlambda[1] 
      lambda_ltre_t[i,j] <- lambda_ltre_t[i,j]+dlambda[2] 
      lambda_ltre_b[i,j] <- lambda_ltre_b[i,j]+dlambda[3] 
    }
    print(j)
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(lambda_ltre_p, main = "Lambda"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

lambda_ltre_env<-sum(lambda_ltre_p,lambda_ltre_t,lambda_ltre_b)

lambda_ltre_p <- ppt_yr_raster
lambda_ltre_p <- setValues(lambda_ltre_p, 0)
lambda_ltre_t <- lambda_ltre_p
lambda_ltre_b <- lambda_ltre_p 

perturb=c(1,1,1)

for (i in 60:nrow(ppt_yr_raster)) {
  for (j in 2:ncol(ppt_yr_raster)){
    
    # Climate at mean
    pred_data <- data.frame(PPT_yr=mean_ppt_pied,
                            T_yr=mean_t_pied,
                            BALIVE=mean_ba_pied)
    
    # Climate at pixel
    pred_data_new <- data.frame(PPT_yr=as.numeric(ppt_yr_raster[i,j]),
                                T_yr=as.numeric(t_yr_raster[i,j]),
                                BALIVE=as.numeric(ba_raster[i,j]))
    
    if (is.na(pred_data_new$PPT_yr) | is.na(pred_data_new$T_yr) | is.na(pred_data_new$BALIVE)) {
      lambda_ltre_p[i,j] <- NA
      lambda_ltre_t[i,j] <- NA
      lambda_ltre_b[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    
    pred_data_mean <- data.frame(PPT_yr=mean(as.numeric(ppt_yr_raster[i,j]),mean_ppt_pied),
                            T_yr=mean(as.numeric(t_yr_raster[i,j]),mean_t_pied),
                            BALIVE=mean(as.numeric(ba_raster[i,j]),mean_ba_pied))
    
    K_old<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam, data=pred_data_mean, 
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    lambda_old<-Re(eigen(K_old)$values[1])
    
    # Calculate change in predictors and divide change in each parameter into amount to perturb each time
    dclim <- c((pred_data_new$PPT_yr-pred_data$PPT_yr),(pred_data_new$T_yr-pred_data$T_yr),
               (pred_data_new$BALIVE-pred_data$BALIVE))
  
    dlambda<-rep(NA,3)
    
    for(p in 1:3){
      pred_data_pert<-pred_data_mean
      pred_data_pert[1,p]<-pred_data_mean[1,p]+perturb[p]
      K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                 rmodel=rmodel.int.gam, gSD=growSD.int.gam, data=pred_data_pert, 
                 s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
      lambda<-Re(eigen(K)$values[1])
      dlambda[p]<-(lambda-lambda_old)/perturb[p]
    }

    lambda_ltre_p[i,j] <- dlambda[1]*dclim[1] 
    lambda_ltre_t[i,j] <- dlambda[2]*dclim[2] 
    lambda_ltre_b[i,j] <- dlambda[3]*dclim[3] 
      
    print(j)
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(lambda_ltre_p, main = "Lambda"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

lambda_ltre_env<-sum(lambda_ltre_p,lambda_ltre_t,lambda_ltre_b)



### Elevation - lambda "LTRE"
mytheme2<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.text=element_text(size=12),legend.title=element_text(size=12),
                legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
                axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
                axis.line.x = element_line(color="black", size = 0.3),
                axis.line.y = element_line(color="black", size = 0.3))

env<-cbind(FIA2$BALIVE,FIA2$PPT_yr_norm,FIA2$T_yr_norm)

FIA2$PLT_CN_factor<-as.factor(FIA2$plot)
k=5

elev_ba <-  bam(BALIVE ~ s(elev,k=k), data = FIA2)
elev_ppt <-  bam(PPT_yr_norm ~ s(elev,k=k), data = FIA2)
elev_t <-  bam(T_yr_norm ~ s(elev,k=k), data = FIA2)

elev_seq <- seq(min(FIA2$elev,na.rm=T),max(FIA2$elev,na.rm=T),by=100)

elev_env<-data.frame(Elevation = elev_seq,
                     BA = predict(elev_ba,newdata = data.frame(elev=elev_seq)),
                     MAP = predict(elev_ppt,newdata = data.frame(elev=elev_seq)),
                     MAT = predict(elev_t,newdata = data.frame(elev=elev_seq)))

save(FIA2,elev_ba,elev_ppt,elev_t,file="./Output/elev_models.rda")

min_elev_pied <- min(subset(FIA2,PApied==1)$elev)
max_elev_pied <- max(subset(FIA2,PApied==1)$elev)
save(min_elev_pied,max_elev_pied, file="./Output/elev_limits.rda")

min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500

lambda_elev<-numeric(0)
ssd_elev<-matrix(NA,n_dim,length(elev_seq))
  
for (i in 1:length(elev_seq)) {
    # Extract climate for cell
    pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"

    # Calculate lambda
    K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.clim.int.gam, smodel=smodel.clim.int.gam, 
               rmodel=rmodel.clim.int.gam, gSD=growSD.clim.int.gam,
               data=pred_data,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    lambda_val <- Re(eigen(K)$values[1])
    ssd<-Re(eigen(K)$vectors[,1])/sum(Re(eigen(K)$vectors[,1]))
    print(lambda_val)
    lambda_elev[i] <- lambda_val
    ssd_elev[,i] <- ssd
}

lambda_elev_ltre<-matrix(NA,length(elev_seq),3)
perturb<-1
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  for(j in 1:3){
    pred_data_new<-pred_data
    pred_data_new[1,j] <- pred_data_new[1,j]+perturb
  # Calculate lambda
  K<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
             rmodel=rmodel.int.gam, gSD=growSD.int.gam,
             data=pred_data_new,
             s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  lambda_val <- Re(eigen(K)$values[1])
  lambda_elev_ltre[i,j] <- lambda_val
  }
  print(lambda_val)
}

d_lambda<-lambda_elev_ltre-lambda_elev

denv_delev<-matrix(NA,length(elev_seq),3)
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)" 
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  pred_data_new<-data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=(elev_seq[i]+1)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            T_yr=predict(elev_t, newdata = data.frame(elev=(elev_seq[i]+1)), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            BALIVE=predict(elev_ba, newdata = data.frame(elev=(elev_seq[i]+1)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
    for(j in 1:3){
    d <- pred_data_new[1,j]-pred_data[1,j]
    # Calculate lambda
    denv_delev[i,j] <- d
    }
  print(d)
}

elev_ltre<-d_lambda*denv_delev

test<-rowSums(elev_ltre)
test2<-(lambda_elev[2:length(elev_seq)]-lambda_elev[1:(length(elev_seq)-1)])/100
elev_seq2<-(elev_seq[1:(length(elev_seq)-1)]+elev_seq[2:length(elev_seq)])/2

plot(elev_seq2,test2)
lines(elev_seq,test)

lam_elev_data_c<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_c$Model="c"
lam_elev_data_ci<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_ci$Model="ci"
lam_elev_data_cc<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_cc$Model="cc"
lam_elev_data_i<-data.frame(Elevation=elev_seq,Lambda=lambda_elev)
lam_elev_data_i$Model="i"

lam_elev_data=rbind(lam_elev_data_c,lam_elev_data_cc,lam_elev_data_i)

write.csv(lam_elev_data,"./Output/lam_elev_data.csv",row.names=F)

ltre_data<-data.frame(Elev=c(elev_seq,elev_seq,elev_seq,elev_seq),
                      Env_c=c(elev_ltre[,1],elev_ltre[,2],elev_ltre[,3],test),
                      dlam=c(d_lambda[,1],d_lambda[,2],d_lambda[,3],rowSums(d_lambda)),
                      denv=c(denv_delev[,1],denv_delev[,2],denv_delev[,3],rowSums(denv_delev)),
                      Predictor=c(rep("MAP",length(elev_seq)),rep("MAT",length(elev_seq)),rep("BA",length(elev_seq)),rep("Total",length(elev_seq))))
ltre_data$Env_cc<-c(elev_ltre[,1],elev_ltre[,2],elev_ltre[,3],test)
ltre_data$Env_i<-c(elev_ltre[,1],elev_ltre[,2],elev_ltre[,3],test)

write.csv(ltre_data,"./Output/ltre_data_env.csv",row.names=F)

# Contribution of vital rates
lambda_elev_ltre_vital<-matrix(NA,length(elev_seq),3)
perturb=0.01
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
    # Calculate lambda
    K_g<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, gperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    K_s<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                 rmodel=rmodel.int.gam, gSD=growSD.int.gam,
                 data=pred_data, sperturb=perturb,
                 s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    K_r<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
                 rmodel=rmodel.int.gam, gSD=growSD.int.gam,
                 data=pred_data, rperturb=perturb,
                 s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
    lambda_elev_ltre_vital[i,1] <- Re(eigen(K_g)$values[1])
    lambda_elev_ltre_vital[i,2] <- Re(eigen(K_s)$values[1])
    lambda_elev_ltre_vital[i,3] <- Re(eigen(K_r)$values[1])
    
  print(i)
}

d_lambda_vital<-(lambda_elev-lambda_elev_ltre_vital)/perturb

dvital_delev<-matrix(NA,length(elev_seq),3)
perturb=10
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)" 
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  # Perturb elevation and calculate new climate
  pred_data_new<-data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=(elev_seq[i]+perturb)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            T_yr=predict(elev_t, newdata = data.frame(elev=(elev_seq[i]+perturb)), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                            BALIVE=predict(elev_ba, newdata = data.frame(elev=(elev_seq[i]+perturb)), #,PLT_CN_factor=14546600020004
                                           type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  # Calculate change in vital rate
  g<-g.mean(model=gmodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data)
  g_perturb<-g.mean(model=gmodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data_new)
  s<-s.x(model=smodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data)
  s_perturb<-s.x(model=smodel.int.gam, size.x = mean(FIA$PREVDIA), data=pred_data_new)
  r<-f.mean(model=rmodel.int.gam, data=pred_data)
  r_perturb<-f.mean(model=rmodel.int.gam, data=pred_data_new)

  dvital_delev[i,1] <- g-g_perturb
  dvital_delev[i,2] <- s-s_perturb
  dvital_delev[i,3] <- r-r_perturb
  
  print(i)
}

elev_ltre_vital<-d_lambda_vital*(dvital_delev)/10

test<-rowSums(elev_ltre_vital)
test2<-(lambda_elev[2:length(elev_seq)]-lambda_elev[1:(length(elev_seq)-1)])/100
elev_seq2<-(elev_seq[1:(length(elev_seq)-1)]+elev_seq[2:length(elev_seq)])/2

plot(elev_seq2,test2)
lines(elev_seq,test)

ltre_data_vital<-data.frame(Elev=c(elev_seq,elev_seq,elev_seq,elev_seq),
                      Env_c=c(elev_ltre_vital[,1],elev_ltre_vital[,2],elev_ltre_vital[,3],test),
                      dlam=c(d_lambda_vital[,1],d_lambda_vital[,2],d_lambda_vital[,3],rowSums(d_lambda_vital)),
                      denv=c((dvital_delev[,1])/10,(dvital_delev[,2])/10,(dvital_delev[,3])/10,rowSums(dvital_delev/10)),
                      Rate=c(rep("Growth",length(elev_seq)),rep("Survival",length(elev_seq)),rep("Recruitment",length(elev_seq)),rep("Total",length(elev_seq))))
ltre_data_vital$Env_cc<-c(elev_ltre_vital[,1],elev_ltre_vital[,2],elev_ltre_vital[,3],test)
ltre_data_vital$Env_i<-c(elev_ltre_vital[,1],elev_ltre_vital[,2],elev_ltre_vital[,3],test)

write.csv(ltre_data_vital,"./Output/ltre_data_vital.csv",row.names=F)

# Elasticity
elev_elast<-matrix(NA,length(elev_seq),3)
perturb=0.01
for (i in 1:length(elev_seq)) {
  # Extract climate for cell
  pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          T_yr=predict(elev_t, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                          BALIVE=predict(elev_ba, newdata = data.frame(elev=elev_seq[i]), #,PLT_CN_factor=14546600020004
                                         type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"
  
  # Calculate lambda
  K_g<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, gperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_s<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, sperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_r<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, rperturb=perturb,
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  elev_elast[i,1] <- Re(eigen(K_g)$values[1])
  elev_elast[i,2] <- Re(eigen(K_s)$values[1])
  elev_elast[i,3] <- Re(eigen(K_r)$values[1])
  
  print(i)
}

#lam_elev_data<-read.csv("./Output/lam_elev_data.csv")

#lambda_elev<-lam_elev_data_i$Lambda #[which(lam_elev_data$Model=="cc")]
elast_vital<-(elev_elast-lambda_elev)/(lambda_elev*perturb)

test<-rowSums(elast_vital)-elast_vital[,1]

elast_data_vital<-data.frame(Elev=c(elev_seq,elev_seq,elev_seq),
                            Elast_c=c(elast_vital[,1],elast_vital[,2],elast_vital[,3]),
                            Rate=c(rep("Growth",length(elev_seq)),rep("Survival",length(elev_seq)),rep("Recruitment",length(elev_seq))))

elast_data_vital$Elast_cc<-c(elast_vital[,1],elast_vital[,2],elast_vital[,3])
elast_data_vital$Elast_i<-c(elast_vital[,1],elast_vital[,2],elast_vital[,3])

write.csv(elast_data_vital,"./Output/elast_vital.csv",row.names=F)


# Perturb vitals
perturb_seq=seq(0,-1,-0.01)
elev_perturb<-matrix(NA,length(perturb_seq),3)
pred_data <- data.frame(PPT_yr=predict(elev_ppt, newdata = data.frame(elev=max_elev_pied), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                        T_yr=predict(elev_t, newdata = data.frame(elev=max_elev_pied), #,PLT_CN_factor=14546600020004
                                     type = "response", re.form = NA), #, exclude = "s(PLT_CN_factor)"
                        BALIVE=predict(elev_ba, newdata = data.frame(elev=max_elev_pied), #,PLT_CN_factor=14546600020004
                                       type = "response", re.form = NA)) #, exclude = "s(PLT_CN_factor)"

for (i in 1:length(perturb_seq)) {

  # Calculate lambda
  K_g<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, gperturb=perturb_seq[i],
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_s<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, sperturb=perturb_seq[i],
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  K_r<-ipm_fun(min=min.size, max=max.size, n=n_dim, gmodel=gmodel.int.gam, smodel=smodel.int.gam, 
               rmodel=rmodel.int.gam, gSD=growSD.int.gam,
               data=pred_data, elas=T, rperturb=perturb_seq[i],
               s.t.clamp=F, g.t.clamp=F, g.ba.clamp=F,r.ba.clamp=F)
  elev_perturb[i,1] <- Re(eigen(K_g)$values[1])
  elev_perturb[i,2] <- Re(eigen(K_s)$values[1])
  elev_perturb[i,3] <- Re(eigen(K_r)$values[1])
  
  print(i)
}


perturb_data_vital<-data.frame(Perturb=c(perturb_seq,perturb_seq,perturb_seq),
                             Elast_c=c(elev_perturb[,1],elev_perturb[,2],elev_perturb[,3]),
                             Rate=c(rep("Growth",length(perturb_seq)),rep("Survival",length(perturb_seq)),rep("Recruitment",length(perturb_seq))))

perturb_data_vital$Elast_cc<-c(elev_perturb[,1],elev_perturb[,2],elev_perturb[,3])
perturb_data_vital$Elast_i<-c(elev_perturb[,1],elev_perturb[,2],elev_perturb[,3])

perturb_data_vital2<-gather(perturb_data_vital,'Elast_c','Elast_cc','Elast_i',key="Model",value="Lambda")

write.csv(perturb_data_vital2,"./Output/perturb_vital.csv",row.names=F)

perturb_plot_growth<-ggplot(data=subset(perturb_data_vital2,Rate=="Growth"),aes(x=(-1*Perturb),y=Lambda,col=Model))+
  geom_abline(intercept=1,slope=0)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Elast_c","Elast_cc","Elast_i"), 
                      values=c("Elast_c"="#1b9e77","Elast_cc"="#d95f02","Elast_i"="#7570b3"),
                      labels=c("Clim","ClimComp","ClimCompInt"),
                      name="Model")+
  labs(x = "Percent perturbation", y ="Lambda")+
  mytheme2
ggsave(file="perturb_plot_growth.png", plot=perturb_plot_growth,width=4,height=3,units="in",dpi=600)

x_c_s<-(-1)*max(subset(perturb_data_vital2,Rate=="Survival" & Model=="Elast_c" & Lambda<1)$Perturb)
x_cc_s<-(-1)*max(subset(perturb_data_vital2,Rate=="Survival" & Model=="Elast_cc" & Lambda<1)$Perturb)
x_i_s<-(-1)*max(subset(perturb_data_vital2,Rate=="Survival" & Model=="Elast_i" & Lambda<1)$Perturb)

perturb_plot_survival<-ggplot(data=subset(perturb_data_vital2,Rate=="Survival"),aes(x=(-1*Perturb),y=Lambda,col=Model))+
  coord_cartesian(xlim = c(0,0.1))+
  geom_abline(intercept=1,slope=0)+
  geom_segment(x=x_c_s,xend=x_c_s,y=1,yend=0,col="#1b9e77",linetype="dashed")+
  geom_segment(x=x_cc_s,xend=x_cc_s,y=1,yend=0,col="#d95f02",linetype="dashed")+
  geom_segment(x=x_i_s,xend=x_i_s,y=1,yend=0,col="#7570b3",linetype="dashed")+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Elast_c","Elast_cc","Elast_i"), 
                      values=c("Elast_c"="#1b9e77","Elast_cc"="#d95f02","Elast_i"="#7570b3"),
                      labels=c("Clim","ClimComp","ClimCompInt"),
                      name="Model")+
  labs(x = "Percent perturbation", y ="Lambda")+
  mytheme2
ggsave(file="perturb_plot_survival.png", plot=perturb_plot_survival,width=4,height=3,units="in",dpi=600)

x_c_r<-(-1)*max(subset(perturb_data_vital2,Rate=="Recruitment" & Model=="Elast_c" & Lambda<1)$Perturb)
x_cc_r<-(-1)*max(subset(perturb_data_vital2,Rate=="Recruitment" & Model=="Elast_cc" & Lambda<1)$Perturb)
x_i_r<-(-1)*max(subset(perturb_data_vital2,Rate=="Recruitment" & Model=="Elast_i" & Lambda<1)$Perturb)

perturb_plot_recruit<-ggplot(data=subset(perturb_data_vital2,Rate=="Recruitment"),aes(x=(-1*Perturb),y=Lambda,col=Model))+
  geom_abline(intercept=1,slope=0)+
  geom_segment(x=x_c_r,xend=x_c_r,y=1,yend=0,col="#1b9e77",linetype="dashed")+
  geom_segment(x=x_cc_r,xend=x_cc_r,y=1,yend=0,col="#d95f02",linetype="dashed")+
  geom_segment(x=x_i_r,xend=x_i_r,y=1,yend=0,col="#7570b3",linetype="dashed")+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Elast_c","Elast_cc","Elast_i"), 
                      values=c("Elast_c"="#1b9e77","Elast_cc"="#d95f02","Elast_i"="#7570b3"),
                      labels=c("Clim","ClimComp","ClimCompInt"),
                      name="Model")+
  labs(x = "Percent perturbation", y ="Lambda")+
  mytheme2
ggsave(file="perturb_plot_recruit.png", plot=perturb_plot_recruit,width=4,height=3,units="in",dpi=600)
