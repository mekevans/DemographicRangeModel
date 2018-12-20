library(raster)
library(dismo)
library(rgdal)

# Load rasters
ppt_c_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
ppt_w_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_warm.tif")
vpd_c_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/VPD_cool.tif")
vpd_w_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/VPD_warm.tif")
ba_raster <- raster("D:/EvansLab/Final/Data/BA/BA_resampled.tif")
elev_raster <- raster("D:/EvansLab/Final/Data/SRTM/SRTM_resampled.tif")

# Stack rasters
pred <- stack(ppt_c_raster, ppt_w_raster, vpd_c_raster, vpd_w_raster, ba_raster, elev_raster)

# Read plots where PIED occurs
occ <- readOGR("D:/EvansLab/Final/Data/Processed", "plotPresences")

# Create MaxEnt model
model <- maxent(x = pred, p = occ)

# Predict for current conditions
prediction <- predict(model, pred, ext = extent(pred), progress = "text")

# Export predictions
writeRaster(prediction, paste0("D:/EvansLab/Final/Output/MaxEnt/MaxEnt.tif"), overwrite = T)

# Plot predictions
plot(prediction)
