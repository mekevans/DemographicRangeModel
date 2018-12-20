library(raster)
library(rgdal)
library(rgeos)

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

# Load climate layers
ppt_c_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
ppt_w_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_warm.tif")
vpd_c_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/VPD_cool.tif")
vpd_w_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/VPD_warm.tif")
ba_raster <- raster("D:/EvansLab/Final/Data/BA/BA_resampled.tif")
elev_raster <- raster("D:/EvansLab/Final/Data/SRTM/SRTM_resampled.tif")

# Aggregate
aggr <- 20
ppt_c_raster <- aggregate(ppt_c_raster, aggr)
ppt_w_raster <- aggregate(ppt_w_raster, aggr)
vpd_c_raster <- aggregate(vpd_c_raster, aggr)
vpd_w_raster <- aggregate(vpd_w_raster, aggr)
ba_raster <- aggregate(ba_raster, aggr)
elev_raster <- aggregate(elev_raster, aggr)

# Prepare empty rasters
lambda <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
lambda <- setValues(lambda, 0)
lambda <- aggregate(lambda, aggr)
growth <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
growth <- setValues(growth, 0)
growth <- aggregate(growth, aggr)
survival <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
survival <- setValues(survival, 0)
survival <- aggregate(survival, aggr)
reproduction <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
reproduction <- setValues(reproduction, 0)
reproduction <- aggregate(reproduction, aggr)

# Build IPMs and calculate lambda
for (i in 1:nrow(ppt_c_raster)) {
  for (j in 1:ncol(ppt_c_raster)) {
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
    lambda[i,j] <- lambda_val
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_c_raster)))
  plot(lambda)
}

# Fill growth, survival, reproduction rasters
for (i in 1:nrow(ppt_c_raster)) {
  for (j in 1:ncol(ppt_c_raster)) {
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
    G <- g.mean(x = median(FIA$DRYBIO_AG), elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
    S <- s.x(x = median(FIA$DRYBIO_AG), elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
    R <- f.mean(x = median(FIA$DRYBIO_AG), elev = elev_val, ba = ba_val, ppt_c = ppt_c_val, ppt_w = ppt_w_val, vpd_c = vpd_c_val, vpd_w = vpd_w_val)
    growth[i,j] <- G
    survival[i,j] <- S
    reproduction[i,j] <- R
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_c_raster)))
  plot(growth)
}

# Crop rasters
fourState <- readOGR("D:/EvansLab/Final/Data/Mask", "4state")
fourState <- gUnaryUnion(fourState)
lambda <- mask(lambda, fourState)
growth <- mask(growth, fourState)
survival <- mask(survival, fourState)
reproduction <- mask(reproduction, fourState)

# Export
writeRaster(lambda, "D:/EvansLab/Final/Output/BC/lambda.tif", overwrite = T)
writeRaster(growth, "D:/EvansLab/Final/Output/BC/growth.tif", overwrite = T)
writeRaster(survival, "D:/EvansLab/Final/Output/BC/survival.tif", overwrite = T)
writeRaster(reproduction, "D:/EvansLab/Final/Output/BC/reproduction.tif", overwrite = T)
pdf("D:/EvansLab/Final/Output/BC/Maps.pdf")
plot(lambda, main = "Lambda")
plot(growth, main = "Growth")
plot(survival, main = "Survival")
plot(reproduction, main = "Reproduction")
dev.off()