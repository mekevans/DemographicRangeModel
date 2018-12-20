library(raster)
library(rgdal)

# Read data
PPT <- list.files(path = "D:/EvansLab/Final/Data/WorldClim/", pattern = "tif$", full.names = T)
PPT <- stack(PPT)

# Create seasonal normals
PPT_cool <- PPT$he85pr7010 + PPT$he85pr7011 + PPT$he85pr7012 + PPT$he85pr701 + PPT$he85pr702 + PPT$he85pr703
PPT_warm <- PPT$he85pr705 + PPT$he85pr706 + PPT$he85pr707 

# Crop normals to extent of BA (for better visualization)
BA <- raster("D:/EvansLab/Final/Data/BA/BA.tif")
ext <- extent(BA)
PPT_cool <- crop(PPT_cool, ext)
PPT_warm <- crop(PPT_warm, ext)

# Resample future rasters
historic <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif")
PPT_cool <- resample(PPT_cool, historic)

historic <- raster("D:/EvansLab/Final/Data/PRISM/Normals/PPT_warm.tif")
PPT_warm <- resample(PPT_warm, historic)

# Export seasonal normals
writeRaster(PPT_cool, "D:/EvansLab/Final/Data/WorldClim/PPT_cool.tif", overwrite = T)
writeRaster(PPT_warm, "D:/EvansLab/Final/Data/WorldClim/PPT_warm.tif", overwrite = T)