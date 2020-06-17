##### Create shapefile of PIED presences
### Must first run dataPrepRecruitment.R to create RecruitmentData.csv
### Last modified: 17 Jun 2020

library(raster)
library(rgdal)

data <- read.csv("./Processed/Recruitment/RecruitmentData.csv")
data <- subset(data, PIED == 1)
proj <- CRS("+proj=longlat +datum=NAD83")
points <- SpatialPointsDataFrame(data.frame(lon = data$lon, lat = data$lat), 
                                 data.frame(PIED = data$PIED),
                                 proj4string = proj)
writeOGR(points, "./Processed", "plotPresences", driver = "ESRI Shapefile")