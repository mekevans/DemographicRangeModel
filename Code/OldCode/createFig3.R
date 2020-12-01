library(raster)
library(pROC)
library(rgeos)
library(rgdal)
library(lattice)
library(ggplot2)
library(cowplot)

# Read presence-absence FIA layer
pa <- raster("D:/EvansLab/Final/Data/Processed/Validation/presenceAbsenceRaster.tif")

# Read in lambda raster
lambda <- raster("D:/EvansLab/Final/Output/BC/lambda.tif")

# Model B
valData <- data.frame(pa = integer(ncell(pa)), lambda = integer(ncell(lambda)))
valData[, "pa"] <- getValues(pa)
valData[, "lambda"] <- getValues(lambda)
valData <- valData[complete.cases(valData),]
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

plotBC <- ggplot(aes(x = lowThreshold, y = lowThreshold), data = valOutput) +
  geom_point(aes(x = lowThreshold, y = highThreshold, col = NRMSE)) +
  xlab("Lower threshold") + ylab("Higher threshold") + scale_colour_gradientn(colors = terrain.colors(1000))

all <- plot_grid(plotB, plotC, plotBC, labels = c("A", "B", "C"), align = "hv")
save_plot("D:/EvansLab/Final/Manuscript/Fig3.png", all, base_aspect_ratio = 2)
