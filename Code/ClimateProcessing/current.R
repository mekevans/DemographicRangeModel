library(raster)
library(rgdal)

# Read PRISM normals
prism_files <- list.files(path = "D:/EvansLab/Final/Data/PRISM/Normals/", pattern = "bil$", full.names = T)
prism <- stack(prism_files)

# Create seasonal normals
PPT_cool <- 
  prism$PRISM_ppt_30yr_normal_800mM2_10_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_11_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_12_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_01_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_02_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_03_bil
VPDmin_cool <- (
  prism$PRISM_vpdmin_30yr_normal_800mM2_10_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_11_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_12_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_01_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_02_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_03_bil) / 6
VPDmax_cool <- (
  prism$PRISM_vpdmax_30yr_normal_800mM2_10_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_11_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_12_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_01_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_02_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_03_bil) / 6
VPD_cool <- (VPDmin_cool + VPDmax_cool) / 2

PPT_warm <- 
  prism$PRISM_ppt_30yr_normal_800mM2_05_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_06_bil + 
  prism$PRISM_ppt_30yr_normal_800mM2_07_bil
VPDmin_warm <- (
  prism$PRISM_vpdmin_30yr_normal_800mM2_05_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_06_bil + 
  prism$PRISM_vpdmin_30yr_normal_800mM2_07_bil) / 3
VPDmax_warm <- (
  prism$PRISM_vpdmax_30yr_normal_800mM2_05_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_06_bil + 
  prism$PRISM_vpdmax_30yr_normal_800mM2_07_bil) / 3
VPD_warm <- (VPDmax_warm + VPDmin_warm) / 2

# Crop normals to extent of BA (for better visualization)
BA <- raster("D:/EvansLab/Final/Data/BA/BA.tif")
ext <- extent(BA)
PPT_cool <- crop(PPT_cool, ext)
PPT_warm <- crop(PPT_warm, ext)
VPD_cool <- crop(VPD_cool, ext)
VPD_warm <- crop(VPD_warm, ext)

# Export seasonal normals
writeRaster(PPT_cool, "D:/EvansLab/Final/Data/PRISM/Normals/PPT_cool.tif", overwrite = T)
writeRaster(PPT_warm, "D:/EvansLab/Final/Data/PRISM/Normals/PPT_warm.tif", overwrite = T)
writeRaster(VPD_cool, "D:/EvansLab/Final/Data/PRISM/Normals/VPD_cool.tif", overwrite = T)
writeRaster(VPD_warm, "D:/EvansLab/Final/Data/PRISM/Normals/VPD_warm.tif", overwrite = T)