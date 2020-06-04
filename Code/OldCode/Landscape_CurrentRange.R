library(raster)

# Read in rasters
lambda <- raster("D:/EvansLab/Final/Output/C/lambda.tif")
presAbs <- raster("D:/EvansLab/Final/Data/Processed/Validation/presenceAbsenceRaster.tif")

# Create data frame
vals <- data.frame(lambda = getValues(lambda), presAbs = getValues(presAbs))
vals <- vals[complete.cases(vals),]

# Create histogram
pdf("D:/EvansLab/Final/Output/C/CurrentVsBackground.pdf")
hist(subset(vals, presAbs == 1)$lambda, col = rgb(0,1,0,0.5), freq = F, main = "Lambda in current (green) vs. background (red) range", xlab = "Lambda")
abline(v = mean(subset(vals, presAbs == 1)$lambda), col = "green", lwd = 10)
hist(vals$lambda, col = rgb(1,0,0,0.5), add = T, freq = F)
abline(v = mean(vals$lambda), col = "red", lwd = 10)
dev.off()