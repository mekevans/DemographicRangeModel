library(raster)
library(rgdal)

# Read PRISM normals
PRISM.norm.path <-  "./ClimateData/PRISM/Normals/"
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
PPT_c <- # cumulative precipitation
  prism$PRISM_ppt_30yr_normal_4kmM2_11_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_12_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_01_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_02_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_03_bil
PPT_pf <- # cumulative precipitation
  prism$PRISM_ppt_30yr_normal_4kmM2_09_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_10_bil
PPT_fs <- # cumulative precipitation
  prism$PRISM_ppt_30yr_normal_4kmM2_04_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_05_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_06_bil
PPT_wd <- # cumulative precipitation
  prism$PRISM_ppt_30yr_normal_4kmM2_09_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_10_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_04_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_05_bil + 
  prism$PRISM_ppt_30yr_normal_4kmM2_06_bil
PPT_m <- # cumulative precipitation
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
T_c <- ( # average T
  prism$PRISM_tmean_30yr_normal_4kmM2_11_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_12_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_01_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_02_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_03_bil) / 5
T_pf <- ( # average T
  prism$PRISM_tmean_30yr_normal_4kmM2_09_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_10_bil) / 2
T_fs <- ( # average T
  prism$PRISM_tmean_30yr_normal_4kmM2_04_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_05_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_06_bil) / 3
T_wd <- ( # average T
  prism$PRISM_tmean_30yr_normal_4kmM2_09_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_10_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_04_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_05_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_06_bil) / 5
T_m <- ( # average T
  prism$PRISM_tmean_30yr_normal_4kmM2_07_bil + 
  prism$PRISM_tmean_30yr_normal_4kmM2_08_bil) / 2

# Crop normals to extent of BA (for better visualization)
BA <- raster("./BA/BA.tif")
ext <- extent(BA)
PPT_year <- crop(PPT_year, ext)
PPT_c <- crop(PPT_c, ext)
PPT_pf <- crop(PPT_pf, ext)
PPT_fs <- crop(PPT_fs, ext)
PPT_wd <- crop(PPT_wd, ext)
PPT_m <- crop(PPT_m, ext)
T_year <- crop(T_year, ext)
T_c <- crop(T_c, ext)
T_pf <- crop(T_pf, ext)
T_fs <- crop(T_fs, ext)
T_wd <- crop(T_wd, ext)
T_m <- crop(T_m, ext)

# Export seasonal normals
writeRaster(PPT_year, paste0(PRISM.norm.path, "PPT_year.tif"), overwrite = T)
writeRaster(PPT_c, paste0(PRISM.norm.path, "PPT_c.tif"), overwrite = T)
writeRaster(PPT_pf, paste0(PRISM.norm.path, "PPT_pf.tif"), overwrite = T)
writeRaster(PPT_fs, paste0(PRISM.norm.path, "PPT_fs.tif"), overwrite = T)
writeRaster(PPT_wd, paste0(PRISM.norm.path, "PPT_wd.tif"), overwrite = T)
writeRaster(PPT_m, paste0(PRISM.norm.path, "PPT_m.tif"), overwrite = T)
writeRaster(T_year, paste0(PRISM.norm.path, "T_year.tif"), overwrite = T)
writeRaster(T_c, paste0(PRISM.norm.path, "T_c.tif"), overwrite = T)
writeRaster(T_pf, paste0(PRISM.norm.path, "T_pf.tif"), overwrite = T)
writeRaster(T_fs, paste0(PRISM.norm.path, "T_fs.tif"), overwrite = T)
writeRaster(T_wd, paste0(PRISM.norm.path, "T_wd.tif"), overwrite = T)
writeRaster(T_m, paste0(PRISM.norm.path, "T_m.tif"), overwrite = T)
