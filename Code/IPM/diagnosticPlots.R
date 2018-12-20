library(ggplot2)
library(cowplot)

# Load models and parameters
load("D:/EvansLab/Final/Models/BC/grow.rda")
load("D:/EvansLab/Final/Models/BC/surv.rda")
load("D:/EvansLab/Final/Models/BC/recr.rda")
load("D:/EvansLab/Final/Models/BC/recrstats.rda")

# Load FIA data
FIA <- read.csv("D:/EvansLab/Final/Data/Processed/Growth/GrowthData.csv")

# Survival
s.x <- function(x, elev, ba, ppt_c, ppt_w, vpd_c, vpd_w) {
  sdata <- data.frame(PREV_DRYBIO_AG = x,
                      elev = elev,
                      baLive = ba,
                      PPT_c = ppt_c,
                      PPT_w = ppt_w,
                      VPD_c = vpd_c, 
                      VPD_w = vpd_w)
  spred <- predict(smodel, sdata, type = "response")
  return(1 - spred)
}

# Growth
g.yx <- function(xp, x, elev, ba, ppt_c, ppt_w, vpd_c, vpd_w) {
  gdata <- data.frame(DRYBIO_AG = x,
                      elev = elev,
                      baLive = ba,
                      PPT_c = ppt_c,
                      PPT_w = ppt_w,
                      VPD_c = vpd_c, 
                      VPD_w = vpd_w)
  gpred <- predict(gmodel, gdata, type = "response")
  return(dnorm(xp - x, gpred, growSD))
}

g.mean <- function(x, elev, ba, ppt_c, ppt_w, vpd_c, vpd_w) {
  gdata <- data.frame(DRYBIO_AG = x,
                      elev = elev,
                      baLive = ba,
                      PPT_c = ppt_c,
                      PPT_w = ppt_w,
                      VPD_c = vpd_c, 
                      VPD_w = vpd_w)
  gpred <- predict(gmodel, gdata, type = "response")
  return(gpred)
}

# Fecundity
fec <- function(xp, x, elev, ba, ppt_c, ppt_w, vpd_c, vpd_w) {
  rdata <- data.frame(AGB_intra = x,
                      measInterval = 1,
                      elev = elev,
                      baLive = ba,
                      PPT_c_window_25 = ppt_c,
                      PPT_w_window_25 = ppt_w,
                      VPD_c_window_25 = vpd_c,
                      VPD_w_window_25 = vpd_w)
  rpred <- predict(rmodel, rdata, type = "response")
  return(dnorm(log(xp), sizemean, sizesd) * rpred)
}

f.mean <- function(x, elev, ba, ppt_c, ppt_w, vpd_c, vpd_w) {
  rdata <- data.frame(AGB_intra = x,
                      measInterval = 1,
                      elev = elev,
                      baLive = ba,
                      PPT_c_window_25 = ppt_c,
                      PPT_w_window_25 = ppt_w,
                      VPD_c_window_25 = vpd_c,
                      VPD_w_window_25 = vpd_w)
  rpred <- predict(rmodel, rdata, type = "response")
  return(rpred)
}
# Set IPM parameters
min.size <- 0.01*min(FIA$DRYBIO_AG, na.rm = T)
max.size <- 1.5*max(FIA$DRYBIO_AG, na.rm = T)
n <- 500
b <- min.size+c(0:n)*(max.size-min.size)/n 
y <- 0.5*(b[1:n]+b[2:(n+1)])
h <- y[2]-y[1]

# Create IPM function
lambdaCalc <- function(elev_val, ba_val, ppt_c_val, ppt_w_val, vpd_c_val, vpd_w_val) {
  # Growth and survival
  G <- h*outer(y, y, g.yx, elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
  S <- s.x(x = y, elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
  P <- G
  for (k in 1:n) P[,k] <- G[,k]*S
  # Recruitment
  R <- h*outer(y, y, fec, elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
  # Entire kernel
  K <- P + R
  # Calculate lambda
  lambda_val <- Re(eigen(K)$values[1])
  print(lambda_val)
  return(lambda_val)
}

# Prepare prediction data frame
noPoints <- 50
predictorNames <- c("elev", "baLive", "PPT_c", "PPT_w", "VPD_c", "VPD_w")
predictors <- FIA[, predictorNames]
predictorDFs <- list()
for (i in 1:ncol(predictors)) {
  currentPredictor <- predictorNames[i]
  print(currentPredictor)
  predictorVals <- data.frame(elev = rep(median(FIA$elev), noPoints), 
                              baLive = rep(median(FIA$baLive), noPoints), 
                              PPT_c = rep(median(FIA$PPT_c), noPoints), 
                              PPT_w = rep(median(FIA$PPT_w), noPoints), 
                              VPD_c = rep(median(FIA$VPD_c), noPoints), 
                              VPD_w = rep(median(FIA$VPD_w), noPoints), 
                              lambda = rep(0, noPoints))
  predictorVals[, i] <- seq(0.5*min(predictors[, i]), 1.5*max(predictors[, i]), length.out = noPoints)
  for (j in 1:nrow(predictorVals)) {
    predictorVals[j, "lambda"] <- lambdaCalc(predictorVals[j, 1], predictorVals[j, 2], predictorVals[j, 3], predictorVals[j, 4], predictorVals[j, 5], predictorVals[j, 6])
  }
  predictorDFs[[i]] <- predictorVals
}

A <- ggplot(data = predictorDFs[[1]], aes(x = elev, y = lambda)) + geom_point() + geom_line() +
  xlab("Elevation") + ylab("Lambda")
B <- ggplot(data = predictorDFs[[2]], aes(x = baLive, y = lambda)) + geom_point() + geom_line() +
  xlab("Plot basal area") + ylab("Lambda")
C <- ggplot(data = predictorDFs[[3]], aes(x = PPT_c, y = lambda)) + geom_point() + geom_line() +
  xlab("Cool-season precipitation") + ylab("Lambda")
D <- ggplot(data = predictorDFs[[4]], aes(x = PPT_w, y = lambda)) + geom_point() + geom_line() +
  xlab("Warm-season precipitation") + ylab("Lambda")
E <- ggplot(data = predictorDFs[[5]], aes(x = VPD_c, y = lambda)) + geom_point() + geom_line() +
  xlab("Cool-season VPD") + ylab("Lambda")
F <- ggplot(data = predictorDFs[[6]], aes(x = VPD_w, y = lambda)) + geom_point() + geom_line() +
  xlab("Warm-season VPD") + ylab("Lambda")
all <- plot_grid(A, B, C, D, F, labels = c("A", "B", "C", "D", "E"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/Fig2.png", all, base_aspect_ratio = 2)

# Plot
pdf("D:/EvansLab/Final/Output/BC/Diagnostics.pdf")
par(mfrow = c(3,2))
for (i in 1:ncol(predictors)) {
  plot(predictorDFs[[i]][, "lambda"] ~ predictorDFs[[i]][, i], main = predictorNames[i], ylab = "lambda", xlab = predictorNames[i])
}
dev.off()