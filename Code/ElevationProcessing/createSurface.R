library(raster)

# Find tiles
tiles <- list.files(path = "D:/EvansLab/Final/Data/SRTM/", pattern = "tif$", full.names = T)

# Read tiles
rasters <- list()
for (i in 1:length(tiles)) {
  rasters[[i]] <- raster(tiles[i])
}

# Create mosaic
rasters$fun <- mean
elev <- do.call(mosaic, rasters)

# Write to file
writeRaster(elev, "D:/EvansLab/Final/Data/SRTM/SRTM_all.tif")