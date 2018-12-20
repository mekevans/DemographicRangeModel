library(raster)
library(rgdal)

data <- read.csv("D:/EvansLab/Final/Data/Processed/Recruitment/RecruitmentData.csv")
data <- subset(data, PIED == 1)
proj <- CRS("+proj=longlat +datum=NAD83")
points <- SpatialPointsDataFrame(data.frame(lon = data$lon, lat = data$lat), 
                                 data.frame(PIED = data$PIED),
                                 proj4string = proj)
writeOGR(points, "D:/EvansLab/Final/Data/Processed", "plotPresences", driver = "ESRI Shapefile")