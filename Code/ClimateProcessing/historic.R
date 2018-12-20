library(raster)

# Search for PRISM files
climateDir <- "D:/EvansLab/Final/Data/PRISM/Historic/"
pptFiles <- list.files(path = climateDir, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
tmpFiles <- list.files(path = climateDir, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
vpdminFiles <- list.files(path = climateDir, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
vpdmaxFiles <- list.files(path = climateDir, pattern = glob2rx("*vpdmax*.bil"), full.names = TRUE)

# Stack monthly data
pptStack <- stack()
for (i in pptFiles) {
  print(i)
  pptStack <- stack(pptStack, raster(i))
}
tmpStack <- stack()
for (i in tmpFiles) {
  print(i)
  tmpStack <- stack(tmpStack, raster(i))
}
vpdStack <- stack()
for (i in 1:length(vpdminFiles)) {
  print(i)
  vpdStack <- stack(vpdStack, (raster(vpdminFiles[i]) + raster(vpdmaxFiles[i])) / 2)
}

# Crop to extent of FIA Pinus edulis occurences
cropExtent <- extent(raster("D:/EvansLab/Final/Data/BA/BA.tif"))
pptStackCropped <- crop(pptStack, cropExtent)
tmpStackCropped <- crop(tmpStack, cropExtent)
vpdStackCropped <- crop(vpdStack, cropExtent)

# Export rasters
writeRaster(pptStackCropped, paste0(climateDir, "pptStack2.tif"), overwrite = T)
writeRaster(tmpStackCropped, paste0(climateDir, "tmpStack2.tif"), overwrite = T)
writeRaster(vpdStackCropped, paste0(climateDir, "vpdStack2.tif"), overwrite = T)