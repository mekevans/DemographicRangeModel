
# Load packages -----------------------------------------------------------

library(mlr)
library(missForest)
library(raster)
library(ranger)

set.seed(123)


# Path structure ---------------------------------------------------------

# path = "C:/Users/mekevans/Documents/old_user/Documents/CDrive/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/"

# climate.path <- "E:/Bayes/DemogRangeMod/ProofOfConcept/FIA-data/westernData/NewData/IWStates/PiedIPM/MEKEvans/ClimateData/"
climate.path = "./ClimateData/"




# Read data ---------------------------------------------------------------


######## basal area
BAdata <- read.csv("./BA/BALIVEdata2.csv", header = T, stringsAsFactors = F)
BAdataold <- read.csv("./BA/BALIVEdata.csv", header = T, stringsAsFactors = F)

######## climate layers
# these created using the script "current.R"
ppt_yr_raster <- raster(paste0(climate.path, "/PRISM/Normals/PPT_year.tif"))
t_yr_raster <- raster(paste0(climate.path, "/PRISM/Normals/T_year.tif"))



# Prepare random forest settings ------------------------------------------


data = na.omit(BAdata[, c("BALIVE", "PPT_yr_norm", "T_yr_norm", "LAT", "LON")])

task = makeRegrTask(data = data, target = "BALIVE")


######## setup learner = model
getParamSet("regr.ranger") # see possible parameter
learner = makeLearner("regr.ranger", importance = "impurity")


######## parameter tuning
params = makeParamSet(makeIntegerParam("mtry", 1,3),
                      makeIntegerParam("min.node.size", 1, 100))

tune_Control = makeTuneControlRandom(maxit = 10) 




# RF tuning ---------------------------------------------------------------

######## tune RF parameters
tune_Result = tuneParams(learner, task, measures = rmse, par.set = params, control = tune_Control,resampling = cv3)
learner_updated = setHyperPars(learner, par.vals = tune_Result$x)




# RF training -------------------------------------------------------------


######## training: 80% train, 20% test
sel = sample(1:nrow(data), ceiling(0.8*nrow(data)))
model = train(learner = learner_updated, task = task, subset = sel)
pred = predict(model, newdata = task$env$data[-sel,])

getFeatureImportance(model)$res/1e6
performance(pred, rmse)
performance(pred, rsq)
# 50.4 % Rsquared in independent sample
#ELS update: 53.1% Rsquared


######## training: 100% train for predicting with rasters
model = train(learner = learner_updated, task = task)
pred = predict(model, task = task) 

getFeatureImportance(model)$res/1e6
performance(pred, rmse)
performance(pred, rsq)
# note: these performance measures are based on the training data and therefore the model fits much better!



# Generate BA map ---------------------------------------------------------

rasterData = data.frame(PPT_yr_norm = raster::getValues(ppt_yr_raster),
                        T_yr_norm = raster::getValues(t_yr_raster),
                        LAT = yFromCell(ppt_yr_raster, 1:ncell(ppt_yr_raster)),
                        LON = xFromCell(ppt_yr_raster, 1:ncell(ppt_yr_raster)))
subs = complete.cases(rasterData)

pred = predict(model, newdata = rasterData[subs, ]) 
rasterPred = rep(NA, nrow(rasterData))
rasterPred[subs] = pred$data$response

balive <- ppt_yr_raster
balive <- setValues(balive, as.vector(rasterPred))

plot(balive)



# Save map ----------------------------------------------------------------

pdf("./BA/BAmap_RF.pdf")
plot(balive, main = "BALIVE")
dev.off()

# as raster (for use in building IPM)
writeRaster(balive, "./BA/balive_RF.tif", overwrite = T)






