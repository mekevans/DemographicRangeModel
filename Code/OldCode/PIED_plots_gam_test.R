#Recruitment
rdata$recruits1yr<-rdata$recruits1/rdata$CENSUS_INTERVAL

r_fun<-function(ba,ppt,t,pa,ci=1,model,clampba=F){
  rdata=data.frame(BALIVE = ifelse(clampba==T & ba >204, 204, 
                                   ifelse(clampba==T & ba <93, 93, ba)),
                   PPT_yr_norm = ppt, T_yr_norm = t)
  #PPT_yr_norm = ppt, T_yr_norm = t)
  scaled.rdata = data.frame(scale(rdata, 
                                  scale = r.scaling$scale[match(names(rdata), 
                                                                names(r.scaling$scale))], 
                                  center = r.scaling$center[match(names(rdata), 
                                                                  names(r.scaling$center))])) 
  scaled.rdata$CENSUS_INTERVAL = ci
  scaled.rdata$PIEDadults1 = pa
  scaled.rdata$off = log(scaled.rdata$CENSUS_INTERVAL) + log(scaled.rdata$PIEDadults1)
  return(predict(model, newdata = scaled.rdata, type = "response"))
}

means<-c(mean(rdata$BALIVE,na.rm=T),
         mean(rdata$PPT_yr_norm,na.rm=T),mean(rdata$T_yr_norm,na.rm=T),
         mean(rdata$PIEDadults1),mean(rdata$CENSUS_INTERVAL))
seq<-data.frame(ba=seq(cellStats(ba_raster,stat="min",na.rm=T),
                       max(rdata$BALIVE,na.rm=T),length=50),
                ppt=seq(cellStats(ppt_yr_raster,stat="min",na.rm=T),
                        cellStats(ppt_yr_raster,stat="max",na.rm=T),length=50),
                t=seq(cellStats(t_yr_raster,stat="min",na.rm=T),
                      cellStats(t_yr_raster,stat="max",na.rm=T),length=50))

#Calculate residuals
rdata$resid_c<-rdata$recruits1yr-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                     rdata$PIEDadults1,1,rmodel.clim.gam)
rdata$resid_ci<-rdata$recruits1yr-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                      rdata$PIEDadults1,1,rmodel.clim.int.gam)
rdata$resid_cc<-rdata$recruits1yr-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                      rdata$PIEDadults1,1,rmodel.clim.comp.gam)
rdata$resid_i<-rdata$recruits1yr-r_fun(rdata$BALIVE,rdata$PPT_yr_norm,rdata$T_yr_norm,
                                     rdata$PIEDadults1,1,rmodel.int.gam)

ncuts=30
chopsize_ba<-cut(rdata$BALIVE,ncuts)
chopsize_PPT<-cut(rdata$PPT_yr_norm,ncuts)
chopsize_T<-cut(rdata$T_yr_norm,ncuts)

ba_binned<-as.vector(sapply(split(rdata$BALIVE,chopsize_ba),mean,na.rm=T))
count_binned_ba<-as.vector(sapply(split(rdata$recruits1,chopsize_ba),length))
recr_binned_ba_c<-as.vector(sapply(split(rdata$resid_c,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],1,rmodel.clim.gam))
recr_binned_ba_ci<-as.vector(sapply(split(rdata$resid_ci,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],1,rmodel.clim.int.gam))
recr_binned_ba_cc<-as.vector(sapply(split(rdata$resid_cc,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],1,rmodel.clim.comp.gam))
recr_binned_ba_i<-as.vector(sapply(split(rdata$resid_i,chopsize_ba),mean,na.rm=T))+
  as.vector(r_fun(ba_binned,means[2],means[3],means[4],1,rmodel.int.gam))

PPT_binned<-as.vector(sapply(split(rdata$PPT_yr_norm,chopsize_PPT),mean,na.rm=T))
count_binned_PPT<-as.vector(sapply(split(rdata$recruits1,chopsize_PPT),length))
recr_binned_PPT_c<-as.vector(sapply(split(rdata$resid_c,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],1,rmodel.clim.gam))
recr_binned_PPT_ci<-as.vector(sapply(split(rdata$resid_ci,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],1,rmodel.clim.int.gam))
recr_binned_PPT_cc<-as.vector(sapply(split(rdata$resid_cc,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],1,rmodel.clim.comp.gam))
recr_binned_PPT_i<-as.vector(sapply(split(rdata$resid_i,chopsize_PPT),mean,na.rm=T))+
  as.vector(r_fun(means[1],PPT_binned,means[3],means[4],1,rmodel.int.gam))

T_binned<-as.vector(sapply(split(rdata$T_yr_norm,chopsize_T),mean,na.rm=T))
count_binned_T<-as.vector(sapply(split(rdata$recruits1,chopsize_T),length))
recr_binned_T_c<-as.vector(sapply(split(rdata$resid_c,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],1,rmodel.clim.gam))
recr_binned_T_ci<-as.vector(sapply(split(rdata$resid_ci,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],1,rmodel.clim.int.gam))
recr_binned_T_cc<-as.vector(sapply(split(rdata$resid_cc,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],1,rmodel.clim.comp.gam))
recr_binned_T_i<-as.vector(sapply(split(rdata$resid_i,chopsize_T),mean,na.rm=T))+
  as.vector(r_fun(means[1],means[2],T_binned,means[4],1,rmodel.int.gam))

r_binned<-as.data.frame(cbind(recr_binned_ba_c,recr_binned_PPT_c,recr_binned_T_c,
                              recr_binned_ba_ci,recr_binned_PPT_ci,recr_binned_T_ci,
                              recr_binned_ba_cc,recr_binned_PPT_cc,recr_binned_T_cc,
                              recr_binned_ba_i,recr_binned_PPT_i,recr_binned_T_i,
                              ba_binned,PPT_binned,T_binned,
                              count_binned_ba,count_binned_PPT,count_binned_T))
names(r_binned)<-c("recr_ba_c","recr_PPT_c","recr_T_c",
                   "recr_ba_ci","recr_PPT_ci","recr_T_ci",
                   "recr_ba_cc","recr_PPT_cc","recr_T_cc",
                   "recr_ba_i","recr_PPT_i","recr_T_i",
                   "BALIVE","PPT","T","count_ba","count_PPT","count_T")

rplot_data_clim<-cbind(data.frame(ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],1,
                                                 rmodel.clim.gam),
                                  t_pred=r_fun(means[1],means[2],seq$t,means[4],1,rmodel.clim.gam),
                                  ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],1,
                                                   rmodel.clim.gam,clampba=T),
                                  t_pred_c=r_fun(means[1],means[2],seq$t,means[4],1,rmodel.clim.gam,
                                                 clampba=T)),seq)
rplot_data_climint<-cbind(data.frame(ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],1,
                                                    rmodel.clim.int.gam),
                                     t_pred=r_fun(means[1],means[2],seq$t,means[4],1,rmodel.clim.int.gam),
                                     ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],1,
                                                      rmodel.clim.int.gam,clampba=T),
                                     t_pred_c=r_fun(means[1],means[2],seq$t,means[4],1,rmodel.clim.int.gam,
                                                    clampba=T)),seq)
rplot_data_climcomp<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],1,
                                                    rmodel.clim.comp.gam),
                                      ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],1,
                                                     rmodel.clim.comp.gam),
                                      t_pred=r_fun(means[1],means[2],seq$t,means[4],1,
                                                   rmodel.clim.comp.gam),
                                      ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],1,
                                                      rmodel.clim.comp.gam,clampba=T),
                                      ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],1,
                                                       rmodel.clim.comp.gam,clampba=T),
                                      t_pred_c=r_fun(means[1],means[2],seq$t,means[4],1,
                                                     rmodel.clim.comp.gam,clampba=T)),seq)
rplot_data_int<-cbind(data.frame(ba_pred=r_fun(seq$ba,means[2],means[3],means[4],1,
                                               rmodel.int.gam),
                                 ppt_pred=r_fun(means[1],seq$ppt,means[3],means[4],1,rmodel.int.gam),
                                 t_pred=r_fun(means[1],means[2],seq$t,means[4],1,rmodel.int.gam),
                                 ba_pred_c=r_fun(seq$ba,means[2],means[3],means[4],1,
                                                 rmodel.int.gam,clampba=T),
                                 ppt_pred_c=r_fun(means[1],seq$ppt,means[3],means[4],1,rmodel.int.gam,
                                                  clampba=T),
                                 t_pred_c=r_fun(means[1],means[2],seq$t,means[4],1,rmodel.int.gam,
                                                clampba=T)),seq)
