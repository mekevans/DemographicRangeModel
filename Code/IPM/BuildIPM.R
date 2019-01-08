library(raster)
library(rgdal)
library(rgeos)

# define path
#path = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/"
path = "C:/Users/mekevans/Documents/Cdrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PIED_IPM/MEKEvans/"

#PRISM.norm.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/Normals/"
PRISM.norm.path <-  "F:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/Normals/"

# load models and scaling --------------------------------------------------

# growth model + scaling
# from modelSelection_Growth.R
load(paste0(path, "Code/IPM/GrRescaling.Rdata"))

# survival model + scaling
# from modelSelection_Survival.R
load(paste0(path, "Code/IPM/SurvRescaling.Rdata"))

# recruitment model + scaling
# from modelSelection_Recruit.R
load(paste0(path, "Code/IPM/RecruitRescaling.Rdata"))

# information on the size distribution of recruits (ingrowth)
# from dataPrepRecruitment.R
load(paste0(path, "Code/IPM/recrstats.rda"))

# Load FIA survival, growth data
FIA <- read.csv(paste0(path, "Processed/Survival/SurvivalData.csv"), header = T, stringsAsFactors = F)
#FIA <- read.csv("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

xy = FIA[, c("LON", "LAT")]
spFIA = sp::SpatialPoints(xy)
spFIA = SpatialPointsDataFrame(spFIA, FIA)

# Survival ---------------------------------------
s.x <- function(size.x, Tann, PPTann, interval = 1) {
  
  # raw (unscaled) data
  # on the LH side is the name of the variable in the cloglog regression (model object)
  # on the RH side is the name of the variable in this function (s.x)
  sdata <- data.frame(PREVDIA = size.x, 
                      #T_yr_norm = Tann, 
                      T_yr_norm = ifelse(Tann < 12, Tann, 12), # clamp response where Tann > 12
                      PPT_yr_norm = PPTann)  
  # rescaled data
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), names(surv.scaling$scale))], # match assures that variables are matched up
                                  center = surv.scaling$center[match(names(sdata), names(surv.scaling$center))])) # regardless of the order they are entered
  
  # add offset to the scaled data
  scaled.sdata = cbind(scaled.sdata,
                       CENSUS_INTERVAL = interval) 
  #scaled.sdata$CENSUS_INTERVAL <- interval # would this be equivalent? consistent with fec below
  
  # apply the fitted model to the scaled input data
  spred <- predict(smodel4.q, newdata = scaled.sdata, type = "response", re.form = NA) # no plot rnd effects used by predict here
  
  return(1-spred)
}

# try out s.x on some realistic combinations of values
#s.x(size.x = c(8.4, 8.4, 8.4), Tann = c(2.0, 8.9, 15.0), PPTann = c(395, 395, 395), interval = 10)

# Growth ------------------------------

# g.yx is used for building the IPM
# it returns the probability density of size at time t+1 as a function size at time t and covariate data
g.yx <- function(size.y, size.x, balive, PPTann, Tann, interval = 1) { # used to be xp, x
  # find scaled value of PREVDIA = 4.5 (for clamping)
  # scale(26, scale=gr.scaling$scale[1], center=gr.scaling$center[1]) # clamp at DRC = 26 inches
  # find scaled value of BALIVE = 1.5 (for clamping)
  # scale(190, scale=gr.scaling$scale[4], center=gr.scaling$center[4]) # clamp at BALIVE = 190
  
  # raw (unscaled) data
  gdata <- data.frame(#PREVDIA = size.x,
                      PREVDIA = ifelse(size.x < 26, size.x, 26),
                      #BALIVE = balive,
                      BALIVE = ifelse(balive < 190, balive, 190),
                      PPT_yr_norm = PPTann,
                      T_yr_norm = Tann)

  # rescaled data
  scaled.gdata = data.frame(scale(gdata, 
                                  scale = gr.scaling$scale[match(names(gdata), names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), names(gr.scaling$center))]))  

  # apply the fitted model to the scaled input data
  gpred <- predict(gmodel.7, newdata = scaled.gdata, re.form = NA) * interval # plot rnd effects not included in prediction
  return(dnorm(size.y - size.x, gpred, growSD)) # returns probability density for quantiles ranging from size.y[1]-size.x to size.y[n]-size.x
  # with mean = gpred and sd = growSD
  # growSD is the SD of the residuals from the growth model
  # i.e., probability of transitioning to a new size at time t+1
  
  }

# try out g.yx
#d_growth <- g.yx(y, 10, balive = 110, PPTann = 348, Tann = 17) # time step = 1 year
d_growth <- g.yx(y, 39, balive = 109, PPTann = 395, Tann = 8.9) # time step = 1 year
plot(y, d_growth, type = "l", ylab = "density", xlim = c(30,40))
# sum of d_growth should = 1.0 (but isn't)

# g.mean is used for making a map (raster) of predicted growth
g.mean <- function(size.x, balive, PPTann, Tann, interval = 1) {
  # raw (unscaled) data
    gdata <- data.frame(#PREVDIA = size.x,
                      PREVDIA = ifelse(size.x < 26, size.x, 26),
                      #BALIVE = balive,
                      BALIVE = ifelse(balive < 190, balive, 190),
                      PPT_yr_norm = PPTann,
                      T_yr_norm = Tann)
  # rescaled data
  scaled.gdata = data.frame(scale(gdata, 
                                  scale = gr.scaling$scale[match(names(gdata), names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), names(gr.scaling$center))]))  
  # apply the fitted model to the scaled input data
  gpred <- predict(gmodel.7, newdata = scaled.gdata, re.form = NA) * interval
  return(gpred)
}

# Fecundity -----------------------------------------------
# function fec is used for building the IPM
fec <- function(size.y, size.x, balive, Tann, PPTann, interval = 1) {

  # raw (unscaled) data
  rdata <- data.frame(BALIVE = balive,
                      T_yr_norm = Tann,
                      PPT_yr_norm = PPTann)
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
  rpred <- predict(rmodel_zipoiss, newdata = scaled.rdata, type = "response")
  #rpred <- predict(rmodel_zip3seas, newdata = scaled.rdata, type = "response") # an alternative model worth considering
  #rpred <- predict(rmodel_zip3seasTint, newdata = scaled.rdata, type = "response") # an alternative model worth considering
  return(dnorm(log(size.y), r.sizemean, r.sizesd) * rpred) # returns the probability density for quantiles ranging from size.y[1] to size.y[n]
  # with mean = sizemean and sd = recruitSD (from the object recruitstats.rda)
  # i.e., the size distribution of recruits
  # multiplied by the predicted number (count) of recruits from the ZIP model, given covariate data
}

# try out fec (y's are the 500 mid-point sizes)
# d_recruit <- fec(y, 10, balive = 110, PPTann = 348, Tann = 17) # time step = 1 year
# d_recruit <- fec(y, 1.0581, balive = 88.6, PPTann = 244, Tann = 7.5) # smallest tree
# d_recruit <- fec(y, 35, balive = 88.6, PPTann = 244, Tann = 7.5) # largest tree
# the probability density looks the same, no matter what the size of the "parent" tree
# because our model of recruitment doesn't account for the influence of tree size
# plot(y, d_recruit, type = "l", ylab = "density")


# function f.mean is used for making a map (raster) of predicted fecundity (recruitment)
f.mean <- function(balive, Tann, PPTann, interval = 1) { 
  # raw (unscaled) data
  rdata <- data.frame(BALIVE = balive,
                      T_yr_norm = Tann,
                      PPT_yr_norm = PPTann)
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
  
  rpred <- predict(rmodel_zipoiss, newdata = scaled.rdata, type = "response")
  return(rpred)
}

# Set IPM parameters
min.size <- 1*min(FIA$PREVDIA, na.rm = T) # minimum size is 1.0 inches diameter at root collar
max.size <- 1.5*max(FIA$PREVDIA, na.rm = T) # maximum size is 39.4 inches DRC
# I think max.size should be smaller...59.1 inches DRC is totally unrealistic
n <- 500 # number of bins (size classes)
b <- min.size+c(0:n)*(max.size-min.size)/n # these are the n+1 "edges" of the size bins
y <- 0.5*(b[1:n]+b[2:(n+1)]) # these are the mid-points of the n size classes
# y <- b[1:n]+0.5*((max.size-min.size)/n)
h <- y[2]-y[1] # bin width

# Load climate layers
# these created using the script "current.R"
ppt_yr_raster <- raster(paste0(PRISM.norm.path, "PPT_year.tif"))
t_yr_raster <- raster(paste0(PRISM.norm.path, "T_year.tif"))
# stand-level basal area raster
ba_raster <- raster(paste0(path, "BA/BA.tif"))
# ba_raster <- raster("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/BA/BA.tif")

ppt_yr_raster <- resample(ppt_yr_raster, ba_raster)
t_yr_raster <- resample(t_yr_raster, ba_raster)

# aggregate for now, just to make it faster
#ppt_yr_raster <- aggregate(ppt_yr_raster, 2)
#t_yr_raster <- aggregate(t_yr_raster, 2)
#ba_raster <- aggregate(ba_raster, 2)


# Prepare empty rasters
lambda <- ppt_yr_raster
lambda <- setValues(lambda, NA)
growth <- ppt_yr_raster
growth <- setValues(growth, 0)
survival <- ppt_yr_raster
survival <- setValues(survival, 0)
reproduction <- ppt_yr_raster
reproduction <- setValues(reproduction, 0)

# Build IPMs and calculate lambda
for (i in 1:nrow(ppt_yr_raster)) {
  for (j in 1:ncol(ppt_yr_raster)) {
    # Extract climate for cell
    ppt_yr_val <- as.numeric(ppt_yr_raster[i,j])
    t_yr_val <- as.numeric(t_yr_raster[i,j])
    ba_val <- as.numeric(ba_raster[i,j])
    # Check for missing value
    if (is.na(ppt_yr_val) | is.na(t_yr_val) | is.na(ba_val)) {
      lambda[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    # Growth and survival
    G <- h*outer(y, y, g.yx, balive = ba_val, PPTann = ppt_yr_val, Tann = t_yr_val) # G is an n*n matrix, currently 500*500 = 250,000
    # each column is one of n sizes at time t (the mid-points), each row is the transition probability to a new size at time t+dt
    #plot(y, G[,10], type = "l", ylab = "density", xlim = c(0,55)) # try with larger and larger values for the G column, the PD marches down/across
    S <- s.x(y, PPTann = ppt_yr_val, Tann = t_yr_val, interval = 1)
    P <- G
    for (k in 1:n) P[,k] <- G[,k]*S #survival*growth subkernel
    # Recruitment
    R <- h*outer(y, y, fec, balive = ba_val, PPTann = ppt_yr_val, Tann = t_yr_val)
    # Entire kernel
    K <- P + R
    # Calculate lambda
    lambda_val <- Re(eigen(K)$values[1])
    print(lambda_val)
    lambda[i,j] <- lambda_val
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(lambda, main = "Lambda"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

# Fill growth, survival, reproduction rasters
for (i in 1:nrow(ppt_yr_raster)) {
  for (j in 1:ncol(ppt_yr_raster)) {
    # Extract climate for cell
    ppt_yr_val <- as.numeric(ppt_yr_raster[i,j])
    t_yr_val <- as.numeric(t_yr_raster[i,j])
    ba_val <- as.numeric(ba_raster[i,j])
    # Check for missing value
    if (is.na(ppt_yr_val) | is.na(t_yr_val) | is.na(ba_val)) {
      lambda[i,j] <- NA
      print("Missing predictor... Skipping!")
      next
    }
    G <- g.mean(size.x = median(FIA$DIA, na.rm = T), balive = ba_val, PPTann = ppt_yr_val, Tann = t_yr_val) # interval = 1 in g.mean function
    S <- s.x(size.x = median(FIA$DIA, na.rm = T), PPTann = ppt_yr_val, Tann = t_yr_val)
    R <- f.mean(balive = ba_val, PPTann = ppt_yr_val, Tann = t_yr_val)
    growth[i,j] <- G
    survival[i,j] <- S
    reproduction[i,j] <- R
  }
  print(paste0("Finished row ", i, " of ", nrow(ppt_yr_raster)))
  plot(growth); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
}

pdf("C:/Users/mekevans/Documents/Cdrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PIED_IPM/MEKEvans/Output/vitalRates.pdf")
plot(growth, main = "Growth"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(survival, main = "Survival", zlim = c(0.95, 1)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(reproduction, main = "Reproduction"); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()


# Crop rasters
# fourState <- readOGR("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans", "4state")
# fourState <- gUnaryUnion(fourState)
# lambda <- mask(lambda, fourState)
# growth <- mask(growth, fourState)
# survival <- mask(survival, fourState)
# reproduction <- mask(reproduction, fourState)

# Export
# writeRaster(lambda, "D:/EvansLab/Final/Output/BC/lambda.tif", overwrite = T)
# writeRaster(growth, "D:/EvansLab/Final/Output/BC/growth.tif", overwrite = T)
# writeRaster(survival, "D:/EvansLab/Final/Output/BC/survival.tif", overwrite = T)
# writeRaster(reproduction, "D:/EvansLab/Final/Output/BC/reproduction.tif", overwrite = T)
