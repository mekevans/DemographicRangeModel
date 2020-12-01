library(raster)
library(rgdal)
library(rgeos)

# Load layers
template_raster <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
ba_raster <- raster("D:/EvansLab/Final/Data/BA/BA.tif")
elev_raster <- raster("D:/EvansLab/Final/Data/SRTM/SRTM_all.tif")

# Resample
ba_raster <- resample(ba_raster, template_raster)
elev_raster <- resample(elev_raster, template_raster)

# Export
writeRaster(ba_raster, "D:/EvansLab/Final/Data/BA/BA_resampled.tif")
writeRaster(elev_raster, "D:/EvansLab/Final/Data/SRTM/SRTM_resampled.tif")