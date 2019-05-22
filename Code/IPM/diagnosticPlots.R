library(ggplot2)
library(cowplot)
library(glmmTMB)

# Load IPM functions
source("./Code/IPM/BuildIPM.R")

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
                              #T_wd_norm = rep(median(FIA$T_wd_norm), noPoints), 
                              #T_c_norm = rep(median(FIA$T_c_norm), noPoints), 
                              #T_m_norm = rep(median(FIA$T_m_norm), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(0.5*min(predictors[, i]), 1.5*max(predictors[, i]), length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- Re(eigen(ipm_fun(min=min.size, max=max.size, n=n_dim, 
                                                   gmodel=gmodel.clim, 
                                                   smodel=smodel.clim, 
                                                   rmodel=rmodel.clim, 
                                                   gSD=growSD.clim,
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