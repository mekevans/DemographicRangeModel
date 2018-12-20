library(raster)
library(rgdal)
library(rgeos)

# Load models and parameters
load("D:/EvansLab/IPM/Code/Elasticities/grow.rda")
load("D:/EvansLab/IPM/Code/Elasticities/surv.rda")
load("D:/EvansLab/IPM/Code/Elasticities/recr.rda")
load("D:/EvansLab/IPM/Code/Elasticities/recrstats.rda")

# Load FIA data
FIA <- read.csv("D:/EvansLab/IPM/Data/Processed/Growth/GrowthData.csv")

# Perturb coefficient
coefN <- 5
oldCoef <- as.numeric(smodel$coefficients[coefN])
smodel$coefficients[coefN] <- smodel$coefficients[coefN] * 1.01

# Survival
s.x <- function(x, vpd_w, elev, ba) {
  sdata <- data.frame(PREV_DRYBIO_AG = x,
                      VPD_w = vpd_w,
                      elev = elev,
                      baLive = ba)
  spred <- predict(smodel, sdata, type = "response")
  return(1 - spred)
}

# Growth
g.yx <- function(xp, x, ba, ppt_w, elev) {
  gdata <- data.frame(DRYBIO_AG = x,
           baLive = ba,
           PPT_w = ppt_w,
           elev = elev)
  gpred <- predict(gmodel, gdata, type = "response")
  return(dnorm(xp - x, gpred, growSD))
}

g.mean <- function(x, ba, ppt_w, elev) {
  gdata <- data.frame(DRYBIO_AG = x,
                      baLive = ba,
                      PPT_w = ppt_w,
                      elev = elev)
  gpred <- predict(gmodel, gdata, type = "response")
  return(gpred)
}

# Fecundity
fec <- function(xp, x, ppt_w, ppt_c, ba, elev) {
  rdata <- data.frame(AGB_intra = x,
                      measInterval = 1,
                      PPT_w_window_25 = ppt_w,
                      PPT_c_window_25 = ppt_c,
                      baLive = ba,
                      elev = elev)
  rpred <- predict(rmodel, rdata, type = "response")
  return(dnorm(log(xp), sizemean, sizesd) * rpred)
}

f.mean <- function(x, ppt_w, ppt_c, ba, elev) {
  rdata <- data.frame(AGB_intra = x,
                      measInterval = 1,
                      PPT_w_window_25 = ppt_w,
                      PPT_c_window_25 = ppt_c,
                      baLive = ba,
                      elev = elev)
  rpred <- predict(rmodel, rdata, type = "response")
  return(rpred)
}

# Set IPM parameters
min.size <- 0.01*min(FIA$DRYBIO_AG, na.rm = T)
max.size <- 1.5*max(FIA$DRYBIO_AG, na.rm = T)
n <- 725
b <- min.size+c(0:n)*(max.size-min.size)/n 
y <- 0.5*(b[1:n]+b[2:(n+1)])
h <- y[2]-y[1]

# Set aggregation factor
aggr <- 20

# Load climate layers
ppt_c_raster <- raster("D:/EvansLab/IPM/Data/Processed/PRISM/PPT_cool.tif")
ppt_w_raster <- raster("D:/EvansLab/IPM/Data/Processed/PRISM/PPT_warm.tif")
vpd_c_raster <- raster("D:/EvansLab/IPM/Data/Processed/PRISM/VPD_cool.tif")
vpd_w_raster <- raster("D:/EvansLab/IPM/Data/Processed/PRISM/VPD_warm.tif")
ba_raster <- raster("D:/EvansLab/IPM/Data/BA/BA_aligned_resampled.tif")
elev_raster <- raster("D:/EvansLab/IPM/Data/SRTM/SRTM_all_resampled.tif")

# Aggregate
ppt_c_raster <- aggregate(ppt_c_raster, aggr)
ppt_w_raster <- aggregate(ppt_w_raster, aggr)
vpd_c_raster <- aggregate(vpd_c_raster, aggr)
vpd_w_raster <- aggregate(vpd_w_raster, aggr)

# Crop BA raster
fourState <- readOGR("D:/EvansLab/IPM/Data/BA", "4state")
fourState <- gUnaryUnion(fourState)
ba_raster <- mask(ba_raster, fourState)

# Prepare empty raster
lambda <- raster("D:/EvansLab/IPM/Data/Processed/PRISM/PPT_cool.tif")
lambda <- setValues(lambda, 0)
lambda <- aggregate(lambda, aggr)

# Build IPMs and calculate lambda
for (i in 1:nrow(ppt_w_raster)) {
  for (j in 1:ncol(ppt_w_raster)) {
    # Extract climate for cell
    ppt_c_val <- as.numeric(ppt_c_raster[i,j])
    ppt_w_val <- as.numeric(ppt_w_raster[i,j])
    vpd_c_val <- as.numeric(vpd_c_raster[i,j])
    vpd_w_val <- as.numeric(vpd_w_raster[i,j])
    ba_val <- as.numeric(ba_raster[i,j])
    elev_val <- as.numeric(elev_raster[i, j])
    # Check for missing value
    if (is.na(ppt_c_val) | is.na(ppt_w_val) | is.na(vpd_c_val) | is.na(vpd_w_val) | is.na(ba_val) | is.na(elev_val)) {
      lambda[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    # Growth and survival
    G <- h*outer(y, y, g.yx, ba = ba_val, ppt_w = ppt_w_val, elev = elev_val)
    S <- s.x(x = y, ba = ba_val, elev = elev_val, vpd_w = vpd_w_val)
    P <- G
    for (k in 1:n) P[,k] <- G[,k]*S
    # Recruitment
    R <- h*outer(y, y, fec, ppt_w = ppt_w_val, ppt_c = ppt_c_val, elev = elev_val, ba = ba_val)
    # Entire kernel
    K <- P + R
    # Calculate lambda
    lambda_val <- Re(eigen(K)$values[1])
    print(lambda_val)
    lambda[i,j] <- lambda_val
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_w_raster)))
  plot(lambda)
}

# Export maps
unperturbed <- raster("D:/EvansLab/IPM/Code/IPM/lambda.tif")
elasticity <- (lambda-unperturbed)/(unperturbed*0.01)
writeRaster(elasticity, "D:/EvansLab/IPM/Code/Elasticities/elasticity_survival_VPD_w.tif", overwrite = T)