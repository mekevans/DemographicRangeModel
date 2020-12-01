#### Use PRISM climate data to create rasters of monthly climate variables in PIED study region

library(raster)

### PRISM download January 22, 2019
### January 1981 through June 2018
### (37*12) + 6 = 450 files

# Search for PRISM files
PRISM.path <-  "./ClimateData/PRISM/"
pptFiles <- list.files(path = PRISM.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
tmpFiles <- list.files(path = PRISM.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
#vpdminFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmin*.bil"), full.names = TRUE)
vpdmaxFiles <- list.files(path = PRISM.path, pattern = glob2rx("*vpdmax*.bil"), full.names = TRUE)

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
for (i in 1:length(vpdmaxFiles)) {
  print(i)
  rast<-raster(vpdmaxFiles[i])
  if(i == 431 | i == 432){crs(rast)<-crs(raster(vpdmaxFiles[1]))}
  vpdStack <- stack(vpdStack, rast)
}

# Crop to extent of FIA Pinus edulis occurrences
cropExtent <- extent(raster("./BA/BA.tif"))
pptStackCropped <- crop(pptStack, cropExtent)
tmpStackCropped <- crop(tmpStack, cropExtent)
vpdStackCropped <- crop(vpdStack, cropExtent)

# Export rasters
clim.path <-  "./ClimateData/"
writeRaster(pptStackCropped, paste0(clim.path, "pptStack.tif"), overwrite = T)
writeRaster(tmpStackCropped, paste0(clim.path, "tmpStack.tif"), overwrite = T)
writeRaster(vpdStackCropped, paste0(clim.path, "vpdStack.tif"), overwrite = T)
