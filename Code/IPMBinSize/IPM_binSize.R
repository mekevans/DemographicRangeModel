library(ggplot2)

# Load models and parameters
load("D:/EvansLab/Final/Models/B/grow.rda")
load("D:/EvansLab/Final/Models/B/surv.rda")
load("D:/EvansLab/Final/Models/B/recr.rda")
load("D:/EvansLab/Final/Models/B/recrstats.rda")

# Set seed
set.seed(2017)

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

# Calculate mean 'climate'
climate <- c(mean(FIA$PPT_c), mean(FIA$PPT_w), mean(FIA$VPD_c), mean(FIA$VPD_w), mean(FIA$baLive), mean(FIA$elev))

# Decide on bin size
convergence <- data.frame(bins = NA, lambdas = NA)
converged <- FALSE
min.size <- 0.01*min(FIA$DRYBIO_AG, na.rm = T)
max.size <- 1.5*max(FIA$DRYBIO_AG, na.rm = T)
n <- 100
b <- min.size+c(0:n)*(max.size-min.size)/n
y <- 0.5*(b[1:n]+b[2:(n+1)])
h <- y[2]-y[1] 
G <- h*outer(y, y, g.yx, elev = mean(FIA$elev), ba = mean(FIA$baLive), ppt_c = mean(FIA$PPT_c), ppt_w = mean(FIA$PPT_w), vpd_c = mean(FIA$VPD_c), vpd_w = mean(FIA$VPD_w))
S <- s.x(x = y, elev = mean(FIA$elev), ba = mean(FIA$baLive), ppt_c = mean(FIA$PPT_c), ppt_w = mean(FIA$PPT_w), vpd_c = mean(FIA$VPD_c), vpd_w = mean(FIA$VPD_w))
P <- G
for (k in 1:n) P[,k] <- G[,k]*S
R <- h*outer(y, y, fec, elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
K <- P + R
lambda <- Re(eigen(K)$values[1])
print(lambda)
convergence[1,"bins"] <- n
convergence[1,"lambdas"] <- lambda
loopCount <- 2
while (converged == FALSE) {
  n <- n + 5
  b <- min.size+c(0:n)*(max.size-min.size)/n
  y <- 0.5*(b[1:n]+b[2:(n+1)])
  h <- y[2]-y[1] 
  G <- h*outer(y, y, g.yx, elev = mean(FIA$elev), ba = mean(FIA$baLive), ppt_c = mean(FIA$PPT_c), ppt_w = mean(FIA$PPT_w), vpd_c = mean(FIA$VPD_c), vpd_w = mean(FIA$VPD_w))
  S <- s.x(x = y, elev = mean(FIA$elev), ba = mean(FIA$baLive), ppt_c = mean(FIA$PPT_c), ppt_w = mean(FIA$PPT_w), vpd_c = mean(FIA$VPD_c), vpd_w = mean(FIA$VPD_w))
  P <- G
  for (k in 1:n) P[,k] <- G[,k]*S
  R <- h*outer(y, y, fec, elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
  K <- P + R
  lambda_new <- Re(eigen(K)$values[1])
  print(lambda_new)
  if (abs(lambda_new - lambda) < 0.001) {
    converged <- TRUE
    print(n)
  }
  lambda <- lambda_new
  convergence[loopCount, "bins"] <- n
  convergence[loopCount, "lambdas"] <- lambda
  loopCount <- loopCount + 1
}

# Plot lambda prediction convergence
ggplot(data = convergence, aes(x = bins, y = lambdas)) +
  geom_point() +
  geom_line() +
  xlab("Number of rows/columns") +
  ylab("Lambda")