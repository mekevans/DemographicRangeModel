### Create vital rate and IPM functions for use in analysis

# Survival ---------------------------------------
#s.x <- function(size.x, Tann, PPTann, interval = 1) {
s.x <- function(model, size.x, data, interval = 1, perturb = 0) {
  
  # raw (unscaled) data
  # on the LH side is the name of the variable in the cloglog regression (model object)
  # on the RH side is the name of the variable in this function (s.x)
  sdata <- data.frame(PREVDIA = size.x, 
                      BALIVE = data$BALIVE,
                      # T_yr_norm = Tann, 
                      #T_yr_norm = ifelse(Tann < 12, Tann, 12), # clamp response where Tann > 12
                      #PPT_yr_norm = PPTann
                      T_yr = ifelse(data$T_yr < 12, data$T_yr, 12), # clamp response where Tann > 12
                      PPT_yr = data$PPT_yr)  
  # rescaled data
  scaled.sdata = data.frame(scale(sdata, 
                                  scale = surv.scaling$scale[match(names(sdata), names(surv.scaling$scale))], # match assures that variables are matched up
                                  center = surv.scaling$center[match(names(sdata), names(surv.scaling$center))])) # regardless of the order they are entered
  
  # add offset to the scaled data
  scaled.sdata = cbind(scaled.sdata,
                       CENSUS_INTERVAL = interval) 
  #scaled.sdata$CENSUS_INTERVAL <- interval # would this be equivalent? consistent with fec below
  
  # apply the fitted model to the scaled input data
  spred <- predict(model, newdata = scaled.sdata, type = "response", re.form = NA) # no plot rnd effects used by predict here
  
  return((1-spred) + perturb)
}

# Growth ------------------------------

# g.yx is used for building the IPM
# it returns the probability density of size at time t+1 as a function size at time t and covariate data
g.yx <- function(model, growSD, size.y, size.x, data, interval = 1, perturb = 0) { # used to be xp, x
  # find scaled value of PREVDIA = 4.5 (for clamping)
  # scale(26, scale=gr.scaling$scale[1], center=gr.scaling$center[1]) # clamp at DRC = 26 inches
  # find scaled value of BALIVE = 1.5 (for clamping)
  # scale(190, scale=gr.scaling$scale[4], center=gr.scaling$center[4]) # clamp at BALIVE = 190
  
  # raw (unscaled) data
  gdata <- data.frame(#PREVDIA = size.x,
    PREVDIA = ifelse(size.x < 26, size.x, 26),
    #BALIVE = balive,
    BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
    PPT_yr = data$PPT_yr,
    #T_yr = data$T_yr
    T_yr = ifelse(data$T_yr < 12, data$T_yr, 12))
  # rescaled data
  scaled.gdata = data.frame(scale(gdata, 
                                  scale = gr.scaling$scale[match(names(gdata), names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), names(gr.scaling$center))]))  
  
  # apply the fitted model to the scaled input data
  gpred <- (predict(model, newdata = scaled.gdata, re.form = NA) * interval) + perturb # plot rnd effects not included in prediction
  return(dnorm(size.y - size.x, gpred, growSD)) # returns probability density for quantiles ranging from size.y[1]-size.x to size.y[n]-size.x
  # with mean = gpred and sd = growSD
  # growSD is the SD of the residuals from the growth model
  # i.e., probability of transitioning to a new size at time t+1
  
}

# g.mean is used for making a map (raster) of predicted growth
g.mean <- function(model, size.x, data, interval = 1, perturb = 0) {
  # raw (unscaled) data
  gdata <- data.frame(#PREVDIA = size.x,
    PREVDIA = ifelse(size.x < 26, size.x, 26),
    #BALIVE = balive,
    BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
    PPT_yr = data$PPT_yr,
    #T_yr = data$T_yr
    T_yr = ifelse(data$T_yr < 12, data$T_yr, 12))
  # rescaled data
  scaled.gdata = data.frame(scale(gdata, 
                                  scale = gr.scaling$scale[match(names(gdata), names(gr.scaling$scale))], 
                                  center = gr.scaling$center[match(names(gdata), names(gr.scaling$center))]))  
  # apply the fitted model to the scaled input data
  gpred <- predict(model, newdata = scaled.gdata, re.form = NA) * interval
  return(gpred + perturb)
}

# Fecundity -----------------------------------------------
# function fec is used for building the IPM
fec <- function(model, size.y, size.x, data, interval = 1, perturb = 0) {
  
  # raw (unscaled) data
  rdata <- data.frame(BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
                      T_yr_norm = data$T_yr,
                      PPT_yr_norm = data$PPT_yr, T_wd_norm = data$T_wd_norm,
                      T_c_norm = data$T_c_norm, T_m_norm = data$T_m_norm)
  # rescaled data
  scaled.rdata = data.frame(scale(rdata, 
                                  scale = r.scaling$scale[match(names(rdata), names(r.scaling$scale))], 
                                  center = r.scaling$center[match(names(rdata), names(r.scaling$center))]))  
  # add census interval offset to the scaled data
  scaled.rdata$CENSUS_INTERVAL <- interval
  
  # size threshold above which trees are treated as reproductively mature
  # threshold defined here is 1 inch DRC, which is not correct biologically
  # but <is> consistent with how the ZIP model of recruitment was built
  # scaled.rdata$PIEDadults1 <- ifelse(size.x >= 1, 1, 0)
  scaled.rdata$PIEDadults1 <- 1 # remove the size threshold
  
  # apply the fitted model to the scaled version of the input data
  rpred <- predict(model, newdata = scaled.rdata, type = "response") + perturb
  #rpred <- predict(rmodel_zip3seas, newdata = scaled.rdata, type = "response") # an alternative model worth considering
  #rpred <- predict(rmodel_zip3seasTint, newdata = scaled.rdata, type = "response") # an alternative model worth considering
  return(dnorm(log(size.y), r.sizemean, r.sizesd) * rpred) # returns the probability density for quantiles ranging from size.y[1] to size.y[n]
  # with mean = sizemean and sd = recruitSD (from the object recruitstats.rda)
  # i.e., the size distribution of recruits
  # multiplied by the predicted number (count) of recruits from the ZIP model, given covariate data
}

 # function f.mean is used for making a map (raster) of predicted fecundity (recruitment)
 f.mean <- function(model, data, interval = 1, perturb = 0) { 
   # raw (unscaled) data
   rdata <- data.frame(BALIVE = ifelse(data$BALIVE < 190, data$BALIVE, 190),
                       T_yr_norm = data$T_yr,
                       PPT_yr_norm = data$PPT_yr, T_wd_norm = data$T_wd_norm,
                       T_c_norm = data$T_c_norm, T_m_norm = data$T_m_norm)
   # rescaled data
   scaled.rdata = data.frame(scale(rdata, 
                                   scale = r.scaling$scale[match(names(rdata), names(r.scaling$scale))], 
                                   center = r.scaling$center[match(names(rdata), names(r.scaling$center))]))  
   # add census interval offset to the scaled data
   scaled.rdata$CENSUS_INTERVAL <- interval
   
   # size threshold above which trees are treated as reproductively mature
   # threshold defined here is 1 inch DRC, which is not correct biologically
   # but <is> consistent with how the ZIP model of recruitment was built
   #scaled.rdata$PIEDadults1 <- ifelse(size.x >= 1, 1, 0)
   scaled.rdata$PIEDadults1 <- 1 # removed the size threshold
   
   rpred <- predict(model, newdata = scaled.rdata, type = "response") + perturb
   return(rpred)
 }
 
 
 ipm_fun<-function(min, max, n=500, gmodel, smodel, rmodel, gSD,
                         data, gperturb=0, sperturb=0, rperturb=0){
   b <- min+c(0:n)*(max-min)/n # these are the n+1 "edges" of the size bins
   y <- 0.5*(b[1:n]+b[2:(n+1)]) # these are the mid-points of the n size classes
   # y <- b[1:n]+0.5*((max.size-min.size)/n)
   h <- y[2]-y[1] # bin width 
   
   # Growth and survival
   G <- h*outer(y, y, g.yx, model = gmodel, growSD = gSD, data = data, perturb = gperturb) # G is an n*n matrix, currently 500*500 = 250,000
   # each column is one of n sizes at time t (the mid-points), each row is the transition probability to a new size at time t+dt
   #plot(y, G[,10], type = "l", ylab = "density", xlim = c(0,55)) # try with larger and larger values for the G column, the PD marches down/across
   # S <- s.x(y, PPTann = ppt_yr_val, Tann = t_yr_val, interval = 1)
   S <- s.x(model = smodel, y, data = data, interval = 1, perturb = sperturb)
   P <- G
   for (k in 1:n) P[,k] <- G[,k]*S #survival*growth subkernel
   
   # Recruitment
   R <- h*outer(y, y, fec, model = rmodel, data = data, perturb = rperturb)
   # Entire kernel
   K <- P + R
   return(K)
 }
 