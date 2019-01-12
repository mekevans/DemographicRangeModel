library(raster)
library(rgdal)

# Read PRISM normals
PRISM.norm.path <-  "./ClimateData/PRISM/Normals"
ppt.norm.files <- list.files(path = PRISM.norm.path, pattern = glob2rx("*ppt*.bil"), full.names = TRUE)
tmp.norm.files <- list.files(path = PRISM.norm.path, pattern = glob2rx("*tmean*.bil"), full.names = TRUE)
vpd.norm.files <- list.files(path = PRISM.norm.path, pattern = glob2rx("*vpdmax*.bil"), full.names = TRUE)
ppt.normals <- stack(ppt.norm.files)
tmp.normals <- stack(tmp.norm.files)
vpd.normals <- stack(vpd.norm.files)

# Crop normals to extent of BA (for better visualization)
cropExtent <- extent(raster("./BA/BA.tif"))
PPT.norm <- crop(ppt.normals, cropExtent)
TMP.norm <- crop(tmp.normals, cropExtent)
VPD.norm <- crop(vpd.normals, cropExtent)

# Export monthly normals
clim.path <-  "./ClimateData/"
writeRaster(PPT.norm, paste0(clim.path, "pptNormals.tif"), overwrite = T)
writeRaster(TMP.norm, paste0(clim.path, "tmpNormals.tif"), overwrite = T)
writeRaster(VPD.norm, paste0(clim.path, "vpdNormals.tif"), overwrite = T)
