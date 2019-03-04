library(ggplot2)
library(cowplot)
library(glmmTMB)

# load models and scaling --------------------------------------------------

# growth model + scaling
# from modelSelection_Growth.R
load("./Code/IPM/GrRescaling.Rdata")

# survival model + scaling
# from modelSelection_Survival.R
# load(paste0(path, "Code/IPM/SurvRescaling.Rdata"))
load("./Code/IPM/SurvRescalingNoFire.Rdata")
#load(paste0(path, "Code/IPM/SurvRescalingBA.Rdata"))

# recruitment model + scaling
# from modelSelection_Recruit.R
load("./Code/IPM/RecruitRescaling.Rdata")

# information on the size distribution of recruits (ingrowth)
# from dataPrepRecruitment.R
load("./Code/IPM/recrstats.rda")

# Load FIA survival, growth data
FIA <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)
# when running IPM without trees killed by fire, technically should filter out those trees
# prolly makes ~0 difference
FIA <- FIA[!(FIA$DSTRBCD1 %in% c(30, 31, 32, 80)), ]

# Survival ---------------------------------------
#s.x <- function(size.x, Tann, PPTann, interval = 1) {
s.x <- function(model, size.x, data, interval = 1) {
  
  # raw (unscaled) data
  # on the LH side is the name of the variable in the cloglog regression (model object)
  # on the RH side is the name of the variable in this function (s.x)
  sdata <- data.frame(PREVDIA = size.x, 
                      BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
                      # T_yr_norm = Tann, 
                      #T_yr_norm = ifelse(Tann < 12, Tann, 12), # clamp response where Tann > 12
                      #PPT_yr_norm = PPTann
                      T_yr = ifelse(data$T_yr < 12, data$T_yr, 12), # clamp response where Tann > 12
                      PPT_yr = data$PPT_yr)  
  # rescaled data
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), names(surv.scaling$scale))], # match assures that variables are matched up
                                  center = surv.scaling$center[match(names(sdata), names(surv.scaling$center))])) # regardless of the order they are entered
  
  # add offset to the scaled data
  scaled.sdata = cbind(scaled.sdata,
                       CENSUS_INTERVAL = interval) 
  #scaled.sdata$CENSUS_INTERVAL <- interval # would this be equivalent? consistent with fec below
  
  # apply the fitted model to the scaled input data
  spred <- predict(model, newdata = scaled.sdata, type = "response", re.form = NA) # no plot rnd effects used by predict here
  
  return(1-spred)
}


g.yx <- function(model, growSD, size.y, size.x, data, interval = 1) { # used to be xp, x
  # find scaled value of PREVDIA = 4.5 (for clamping)
  # scale(26, scale=gr.scaling$scale[1], center=gr.scaling$center[1]) # clamp at DRC = 26 inches
  # find scaled value of BALIVE = 1.5 (for clamping)
  # scale(190, scale=gr.scaling$scale[4], center=gr.scaling$center[4]) # clamp at BALIVE = 190
  
  # raw (unscaled) data
  gdata <- data.frame(#PREVDIA = size.x,
    PREVDIA = ifelse(size.x < 26, size.x, 26),
    #BALIVE = balive,
    BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
    PPT_yr = data$PPT_yr,
    T_yr = data$T_yr)
  # rescaled data
  scaled.gdata = data.frame(scale(gdata, 
                                  scale = gr.scaling$scale[match(names(gdata), names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), names(gr.scaling$center))]))  
  
  # apply the fitted model to the scaled input data
  gpred <- predict(model, newdata = scaled.gdata, re.form = NA) * interval # plot rnd effects not included in prediction
  return(dnorm(size.y - size.x, gpred, growSD)) # returns probability density for quantiles ranging from size.y[1]-size.x to size.y[n]-size.x
  # with mean = gpred and sd = growSD
  # growSD is the SD of the residuals from the growth model
  # i.e., probability of transitioning to a new size at time t+1
  
}

# g.mean is used for making a map (raster) of predicted growth
g.mean <- function(model, size.x, data, interval = 1) {
  # raw (unscaled) data
  gdata <- data.frame(#PREVDIA = size.x,
    PREVDIA = ifelse(size.x < 26, size.x, 26),
    #BALIVE = balive,
    BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
    PPT_yr = data$PPT_yr,
    T_yr = data$T_yr)
  # rescaled data
  scaled.gdata = data.frame(scale(gdata, 
                                  scale = gr.scaling$scale[match(names(gdata), names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), names(gr.scaling$center))]))  
  # apply the fitted model to the scaled input data
  gpred <- predict(model, newdata = scaled.gdata, re.form = NA) * interval
  return(gpred)
}

# Fecundity -----------------------------------------------
# function fec is used for building the IPM
fec <- function(model, size.y, size.x, data, interval = 1) {
  
  # raw (unscaled) data
  rdata <- data.frame(BALIVE = data$BALIVE, T_yr_norm = data$T_yr,
                      PPT_yr_norm = data$PPT_yr, T_wd_norm = data$T_wd_norm,
                      T_c_norm = data$T_c_norm, T_m_norm = data$T_m_norm)
  # rescaled data
  scaled.rdata = data.frame(scale(rdata, 
                                  scale = r.scaling$scale[match(names(rdata), names(r.scaling$scale))], 
                                  center = r.scaling$center[match(names(rdata), names(r.scaling$center))]))  
  # add census interval offset to the scaled data
  scaled.rdata$CENSUS_INTERVAL <- interval
  
  # size threshold above which trees are treated as reproductively mature
  # threshold defined here is 1 inch DRC, which is not correct biologically
  # but <is> consistent with how the ZIP model of recruitment was built
  # scaled.rdata$PIEDadults1 <- ifelse(size.x >= 1, 1, 0)
  scaled.rdata$PIEDadults1 <- 1 # remove the size threshold
  
  # apply the fitted model to the scaled version of the input data
  rpred <- predict(model, newdata = scaled.rdata, type = "response")
  #rpred <- predict(rmodel_zip3seas, newdata = scaled.rdata, type = "response") # an alternative model worth considering
  #rpred <- predict(rmodel_zip3seasTint, newdata = scaled.rdata, type = "response") # an alternative model worth considering
  return(dnorm(log(size.y), r.sizemean, r.sizesd) * rpred) # returns the probability density for quantiles ranging from size.y[1] to size.y[n]
  # with mean = sizemean and sd = recruitSD (from the object recruitstats.rda)
  # i.e., the size distribution of recruits
  # multiplied by the predicted number (count) of recruits from the ZIP model, given covariate data
}

# function f.mean is used for making a map (raster) of predicted fecundity (recruitment)
f.mean <- function(model, data, interval = 1) { 
  # raw (unscaled) data
  rdata <- data.frame(BALIVE = data$BALIVE, T_yr_norm = data$T_yr,
                      PPT_yr_norm = data$PPT_yr, T_wd_norm = data$T_wd_norm,
                      T_c_norm = data$T_c_norm, T_m_norm = data$T_m_norm)
  # rescaled data
  scaled.rdata = data.frame(scale(rdata, 
                                  scale = r.scaling$scale[match(names(rdata), names(r.scaling$scale))], 
                                  center = r.scaling$center[match(names(rdata), names(r.scaling$center))]))  
  # add census interval offset to the scaled data
  scaled.rdata$CENSUS_INTERVAL <- interval
  
  # size threshold above which trees are treated as reproductively mature
  # threshold defined here is 1 inch DRC, which is not correct biologically
  # but <is> consistent with how the ZIP model of recruitment was built
  #scaled.rdata$PIEDadults1 <- ifelse(size.x >= 1, 1, 0)
  scaled.rdata$PIEDadults1 <- 1 # removed the size threshold
  
  rpred <- predict(model, newdata = scaled.rdata, type = "response")
  return(rpred)
}

# Set IPM parameters
ipm_fun<-function(min, max, n=500, gmodel, smodel, rmodel, gSD,
                  data){
  b <- min.size+c(0:n)*(max.size-min.size)/n # these are the n+1 "edges" of the size bins
  y <- 0.5*(b[1:n]+b[2:(n+1)]) # these are the mid-points of the n size classes
  # y <- b[1:n]+0.5*((max.size-min.size)/n)
  h <- y[2]-y[1] # bin width 
  
  # Growth and survival
  G <- h*outer(y, y, g.yx, model = gmodel, growSD = gSD, data = data) # G is an n*n matrix, currently 500*500 = 250,000
  # each column is one of n sizes at time t (the mid-points), each row is the transition probability to a new size at time t+dt
  #plot(y, G[,10], type = "l", ylab = "density", xlim = c(0,55)) # try with larger and larger values for the G column, the PD marches down/across
  # S <- s.x(y, PPTann = ppt_yr_val, Tann = t_yr_val, interval = 1)
  S <- s.x(model = smodel, y, data = data, interval = 1)
  P <- G
  for (k in 1:n) P[,k] <- G[,k]*S #survival*growth subkernel
  
  # Recruitment
  R <- h*outer(y, y, fec, model = rmodel, data = data)
  # Entire kernel
  K <- P + R
  return(K)
}

min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
n_dim <- 500

# Prepare prediction data frame
noPoints <- 50
predictorNames <- c("BALIVE", "PPT_yr_norm", "T_yr_norm") #, "T_wd_norm", "T_c_norm", "T_m_norm")
predictors <- FIA[, predictorNames]
predictorDFs <- list()
for (i in 1:ncol(predictors)) {
  currentPredictor <- predictorNames[i]
  print(currentPredictor)
  predictorVals <- data.frame(BALIVE = rep(median(FIA$BALIVE), noPoints), 
                              PPT_yr = rep(median(FIA$PPT_yr_norm), noPoints), 
                              T_yr = rep(median(FIA$T_yr_norm), noPoints), 
                              T_wd_norm = rep(median(FIA$T_wd_norm), noPoints), 
                              T_c_norm = rep(median(FIA$T_c_norm), noPoints), 
                              T_m_norm = rep(median(FIA$T_m_norm), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(0.5*min(predictors[, i]), 1.5*max(predictors[, i]), length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- Re(eigen(ipm_fun(min=min.size, max=max.size, n=n_dim, 
                                                   gmodel=gmodel_int_q, 
                                                   smodel=sbase_int_q, 
                                                   rmodel=r_int_q, 
                                                   gSD=growSD_int_q,
                                          data=predictorVals[j,]))$values[1])
  }
  predictorDFs[[i]] <- predictorVals
}

A <- ggplot(data = predictorDFs[[1]], aes(x = BALIVE, y = lambda)) + geom_point() + geom_line() +
  xlab("Plot basal area") + ylab("Lambda")
B <- ggplot(data = predictorDFs[[2]], aes(x = PPT_yr, y = lambda)) + geom_point() + geom_line() +
  xlab("MAP") + ylab("Lambda")
C <- ggplot(data = predictorDFs[[3]], aes(x = T_yr, y = lambda)) + geom_point() + geom_line() +
  xlab("MAT") + ylab("Lambda")
D <- ggplot(data = predictorDFs[[4]], aes(x = T_wd_norm, y = lambda)) + geom_point() + geom_line() +
  xlab("Warm/dry-season temp") + ylab("Lambda")
E <- ggplot(data = predictorDFs[[5]], aes(x = T_c_norm, y = lambda)) + geom_point() + geom_line() +
  xlab("Cool-season temp") + ylab("Lambda")
F <- ggplot(data = predictorDFs[[6]], aes(x = T_m_norm, y = lambda)) + geom_point() + geom_line() +
  xlab("Monsoon temp") + ylab("Lambda")

all <- plot_grid(A, B, C, labels = c("a", "b", "c"), align = "hv")
save_plot("./Output/Lambda_plots_int_q.png", all, base_aspect_ratio = 2)

# Plot
pdf("D:/EvansLab/Final/Output/BC/Diagnostics.pdf")
par(mfrow = c(3,2))
for (i in 1:ncol(predictors)) {
  plot(predictorDFs[[i]][, "lambda"] ~ predictorDFs[[i]][, i], main = predictorNames[i], ylab = "lambda", xlab = predictorNames[i])
}
dev.off()