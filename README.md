# DemographicRangeModel
Below are descriptions of the scripts used to make a range-wide integral projection model for Pinus edulis based on USFS forest inventory data

Preprocessing
1)	Code/BasalArea/createSurface.R: creates a points shapefile of plot basal area to be used as input for creating a spatial surface of basal area; uses FIA data; creating basal area map is done using interpolation in ArcGIS
2)	Code/CensusProcessing/censusProcessing.R: takes all FIA plots, subsets AZ, CO, NM, UT plots; combines all tree and condition records for these states as well
3)	Code/CensusProcessing/createOccurrences.R: takes the processed recruitment data, and uses it to create a points shapefile of pinon occurrences   
4)	Code/ClimateProcessing/current.R: creates annual and seasonal normals from PRISM data for current time
5)	Code/ClimateProcessing/normals.R: creates monthly normals from PRISM data for current time
6)	Code/ClimateProcessing/future.R: creates seasonal normals from WorldClim data for future time
7)	Code/ClimateProcessing/historic.R: creates raster stacks from PRISM data for historic time
8)	Code/ElevationProcessing/createSurface.R: creates elevation raster from SRTM

Vital rate models

Growth
1)	Code/Growth/dataPrepGrowth.R: creates data frame for input variables to create growth model
2)	Code/Growth/modelSelection.R: selects growth models and creates figures (e.g. marginal effect plots)

Survival
1)	Code/Survival/dataPrepSurvival.R: creates data frame for input variables to create survival model
2)	Code/Survival/modelSelection.R: selects survival models and creates figures (e.g. marginal effect plots)

Recruitment
1)	Code/Recruitment/dataPrepRecruitment.R: creates data frame for input variables to create recruitment model
2)	Code/Recruitment/modelSelection.R: selects recruitment models and creates figures (e.g. coefficient plot, observed versus predicted recruitment)
3)	Code/Recruitment/recruitSize.R: calculates mean size and SD of size of recruits (needed for kernel)

GAMS
1)	Code/GAMs/vital_gams.R: runs GAM models for growth, survival, and recruitment. Saves output to be used in IPMs.

IPM
1)	Code/IPMBinSize/IPM_binSize.R: selects a bin size for the IPM kernel
2)	Code/IPM/BuildIPM.R: creates the vital rate functions and the function for the IPM kernal
3)	Code/IPM/EvalIPM.R: uses the IPM function from BuildIPM.R and predictor rasters to do spatial projections of lambda
4)	Code/IPM/diagnosticPlots.R: creates plots of lambda and vital rates as a function of input variables
5)	Code/IPM/LTRE.R: uses the IPM function from BuildIPM.R and predictor rasters to do life table response experiment and elasticity analysis
6)	Code/IPM/Bootstrapping.R: uses the IPM function from BuildIPM.R and predictor rasters to calculate uncertainty around lambda estimates, using bootstrapping

Validation
1)	Code/Validation/prepFIAPresenceAbsence.R: creates a raster of cells coding pinon presence/absence
2)	Code/Validation/ResidualAnalysis.R: uses output of EvalIPM.R and prepFIAPresenceAbsence to validate models using occurrence data and residual analysis

Data visualization
1) Code/Plots/PIED_plots_resid.R: creates plots of predictors vs vital rates, using vital rate models (linear models, with quadratic terms)
2) Code/Plots/PIED_plots_lin_resid.R: creates plots of predictors vs vital rates, using vital rate models (linear models, without quadratic terms)
3) Code/Plots/PIED_plots_gam.R: generates data frames with output of vital rate models (GAMS) and predictors to use for plots of predictors vs vital rates; these dataframes get passed to PIED_plots_all.R
4) Code/Plots/PIED_plots_all.R: uses data from IPM scripts, validation scripts, and PIED_plots_all.R to create figures used in the manuscript Schultz et al.

Fire
1) Code/Fire/Fire_prelim.R: uses FIA data to calculate mortality rates of Pinus edulis in surface and crown fires used in Discussion in Schultz et al.
