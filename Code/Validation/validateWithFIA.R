library(raster)
library(pROC)
library(rgeos)
library(rgdal)
library(lattice)
library(ggplot2)
library(cowplot)

# # Set aggregation factor
# aggFactor <- 1

# Read presence-absence FIA layer
pa <- raster("D:/EvansLab/Final/Data/Processed/Validation/presenceAbsenceRaster.tif")

# Read in lambda rasters
lambda <- raster("D:/EvansLab/Final/Output/C/lambda.tif")
plot(lambda)

# # Aggregate
# if (aggFactor != 1) {
#   pa <- aggregate(pa, aggFactor, fun = max)
#   lambda <- aggregate(lambda, aggFactor)
# }

# AUC
# valData <- data.frame(pa = integer(ncell(pa)), lambda = integer(ncell(lambda)))
# valData[, "pa"] <- getValues(pa)
# valData[, "lambda"] <- getValues(lambda)
# valData <- valData[complete.cases(valData),]
# AUC <- roc(response = valData$pa, predictor = valData$lambda)

# # RMSE (simple)
# lambdaMin <- cellStats(lambda, min)
# lambdaMax <- cellStats(lambda, max)
# lambdaThresholds <- seq(lambdaMin, lambdaMax, 0.001)
# RMSEs <- numeric(length(lambdaThresholds))
# NRMSEs <- numeric(length(lambdaThresholds))
# 
# for (i in 1:length(lambdaThresholds)) {
#   threshold <- lambdaThresholds[i]
#   valData$thresholded <- ifelse(valData$lambda < threshold, 0, 1)
#   valData$error <- valData$thresholded - valData$pa
#   rmse <- sqrt(sum(valData$error^2) / nrow(valData))
#   nrmse <- rmse / mean(valData$pa)
#   RMSEs[i] <- rmse
#   NRMSEs[i] <- nrmse
#   print(i)
# }
# 
# pdf("D:/EvansLab/Final/Output/C/RMSE_simple.pdf")
# plot(NRMSEs ~ lambdaThresholds)
# dev.off()

# RMSE (complex)
lambdaMin <- cellStats(lambda, min)
lambdaMax <- cellStats(lambda, max)
lambdaThresholds <- seq(lambdaMin, lambdaMax, 0.005)
RMSEs <- numeric(length(lambdaThresholds))
NRMSEs <- numeric(length(lambdaThresholds))
valOutput <- data.frame(lowThreshold = numeric(), highThreshold = numeric(), RMSE = numeric(), NRMSE = numeric())

counter <- 1
for (i in 1:length(lambdaThresholds)) {
  for (j in 1:length(lambdaThresholds)) {
    lowThreshold <- lambdaThresholds[i]
    highThreshold <- lambdaThresholds[j]
    if (j <= i) next
    valOutput[counter, "lowThreshold"] <- lowThreshold
    valOutput[counter, "highThreshold"] <- highThreshold
    valData$thresholded <- ifelse(valData$lambda < highThreshold & valData$lambda > lowThreshold, 1, 0)
    valData$error <- valData$thresholded - valData$pa
    rmse <- sqrt(sum(valData$error^2) / nrow(valData))
    nrmse <- rmse / mean(valData$pa)
    valOutput[counter, "RMSE"] <- rmse
    valOutput[counter, "NRMSE"] <- nrmse
    counter <- counter + 1
  }
  print(i)
}

# pdf("D:/EvansLab/Final/Output/C/RMSE_complex.pdf")
# wireframe(NRMSE ~ lowThreshold*highThreshold, data = valOutput, main = "NRMSE landscape",
#           drape = T,
#           colorkey = T)
# ggplot(aes(x = lowThreshold, y = lowThreshold), data = valOutput) +
#   geom_point(aes(x = lowThreshold, y = highThreshold, col = NRMSE)) +
#   xlab("Lower threshold") + ylab("Higher threshold") + scale_colour_gradientn(colors = terrain.colors(1000))
# dev.off()

# AUC (quadrant)
# nwObserved <- getValuesBlock(pa, row = 1, nrows = nrow(pa) / 2 + 0.5, col = 1, ncols = ncol(pa) / 2)
# nwPredicted <- getValuesBlock(lambda, row = 1, nrows = nrow(lambda) / 2 + 0.5, col = 1, ncols = ncol(lambda) / 2)
# nwAUC <- roc(response = nwObserved, predictor = nwPredicted)
# 
# neObserved <- getValuesBlock(pa, row = 1, nrows = nrow(pa) / 2 + 0.5, col = (ncol(pa) / 2) + 1, ncols = (ncol(pa) / 2))
# nePredicted <- getValuesBlock(lambda, row = 1, nrows = nrow(lambda) / 2 + 0.5, col = (ncol(lambda) / 2) + 1, ncols = (ncol(lambda) / 2))
# neAUC <- roc(response = neObserved, predictor = nePredicted)
# 
# swObserved <- getValuesBlock(pa, row = (nrow(pa) / 2) + 1.5, nrows = nrow(pa), col = 1, ncols = ncol(pa) / 2)
# swPredicted <- getValuesBlock(lambda, row = (nrow(lambda) / 2) + 1, nrows = nrow(lambda), col = 1, ncols = ncol(lambda) / 2)
# swAUC <- roc(response = swObserved, predictor = swPredicted)
# 
# seObserved <- getValuesBlock(pa, row = (nrow(pa) / 2) + 1.5, nrows = nrow(pa), col = (ncol(pa) / 2) + 1, ncols = (ncol(pa) / 2))
# sePredicted <- getValuesBlock(lambda, row = (nrow(lambda) / 2) + 1, nrows = nrow(lambda), col = (ncol(lambda) / 2) + 1, ncols = (ncol(lambda) / 2))
# seAUC <- roc(response = seObserved, predictor = sePredicted)

# Print AUC
# print(AUC)
# print(nwAUC)
# print(neAUC)
# print(swAUC)
# print(seAUC)