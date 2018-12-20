library(raster)

lambdaB <- raster("D:/EvansLab/Final/Output/B/lambda.tif")
lambdaC <- raster("D:/EvansLab/Final/Output/C/lambda.tif")
lambdaBC <- raster("D:/EvansLab/Final/Output/BC/lambda.tif")

diff <- lambdaB - (lambdaC + lambdaBC) / 2

plot(diff)
writeRaster(diff, "D:/EvansLab/Final/Output/rawS12.tif")
