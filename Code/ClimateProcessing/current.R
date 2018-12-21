library(raster)
library(rgdal)

# Read PRISM normals
PRISM.norm.path <-  "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/PRISM/Normals/"
prism_files <- list.files(path = PRISM.norm.path, pattern = "bil$", full.names = T)
prism <- stack(prism_files)

# Create water year normals
PPT_year <- # cumulative precipitation
  prism$PRISM_ppt_30yr_normal_4kmM2_09_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_10_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_11_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_12_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_01_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_02_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_03_bil +
  prism$PRISM_ppt_30yr_normal_4kmM2_04_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_05_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_06_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_07_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_08_bil
T_year <- ( # average T
  prism$PRISM_tmean_30yr_normal_4kmM2_09_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_10_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_11_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_12_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_01_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_02_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_03_bil +
  prism$PRISM_tmean_30yr_normal_4kmM2_04_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_05_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_06_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_07_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_08_bil) / 12

# Crop normals to extent of BA (for better visualization)
BA <- raster("C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/BA/BA.tif")
ext <- extent(BA)
PPT_year <- crop(PPT_year, ext)
T_year <- crop(T_year, ext)

# Export seasonal normals
writeRaster(PPT_year, paste0(PRISM.norm.path, "PPT_year.tif"), overwrite = T)
writeRaster(T_year, paste0(PRISM.norm.path, "T_year.tif"), overwrite = T)
