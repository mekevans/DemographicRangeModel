#### PIED Plots ####

## Code to create figures for PIED manuscript "Demographic range modeling reveals that climate and competition are insufficient to explain a species’ distribution"
## Created by: Emily Schultz
## Created on: 18 Dec 2019
## Last modified: 18 Ded 2019

### Packages
library(tidyverse)
library(plotrix)
library(cowplot)
library(RColorBrewer)
library(scales)
library(ggalt)
library(raster)

### Set up theme for ggplot figures
mytheme<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.text=element_text(size=11),legend.title=element_text(size=12),
                legend.key = element_rect(fill = "white"),axis.text=element_text(size=12),
                axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),
                axis.line.x = element_line(color="black", size = 0.3),
                axis.line.y = element_line(color="black", size = 0.3))

legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

### Predictor covariance ----

## Data
load("./Output/covariance.rda")

## PIED present
ba_ppt_p<-ggplot(data=FIA_pres,aes(x=BALIVE,y=PPT_yr_norm))+
  geom_point()+
  mytheme+labs(x="Basal area",y="MAP")+
  annotate("label", x = 200, y = max(FIA_pres$PPT_yr_norm), 
           label = paste("Covariance =",round(pres_cov[1,3],2),"\nCorrelation =",round(pres_cor[1,3],2)),
           hjust = 0, vjust = 1, colour="#d95f02",size=5)
ba_t_p<-ggplot(data=FIA_pres,aes(x=BALIVE,y=T_yr_norm))+
  geom_point()+
  mytheme+labs(x="Basal area",y="MAT")+
  annotate("label", x = 200, y = max(FIA_pres$T_yr_norm), 
           label = paste("Covariance=",round(pres_cov[2,3],2),"\nCorrelation =",round(pres_cor[2,3],2)),
           hjust = 0, vjust = 1, colour="#d95f02",size=5)
ppt_t_p<-ggplot(data=FIA_pres,aes(x=PPT_yr_norm,y=T_yr_norm))+
  geom_point()+
  mytheme+labs(x="MAP",y="MAT")+
  annotate("label", x = 525, y = max(FIA_pres$T_yr_norm), 
           label = paste("Covariance=",round(pres_cov[1,2],2),"\nCorrelation =",round(pres_cor[1,2],2)),
           hjust = 0, vjust = 1, colour="#d95f02",size=5)

## PIED absent
ba_ppt_a<-ggplot(data=FIA_abs,aes(x=BALIVE,y=PPT_yr_norm))+
  geom_point()+
  mytheme+labs(x="Basal area",y="MAP")+
  annotate("label", x = 375, y = max(FIA_abs$PPT_yr_norm), 
           label = paste("Covariance=",round(abs_cov[1,3],2),"\nCorrelation =",round(abs_cor[1,3],2)),
           hjust = 0, vjust = 1, colour="#1b9e77",size=5)
ba_t_a<-ggplot(data=FIA_abs,aes(x=BALIVE,y=T_yr_norm))+
  geom_point()+
  mytheme+labs(x="Basal area",y="MAT")+
  annotate("label", x = 375, y = max(FIA_abs$T_yr_norm), 
           label = paste("Covariance=",round(abs_cov[2,3],2),"\nCorrelation =",round(abs_cor[2,3],2)),
           hjust = 0, vjust = 1, colour="#d95f02",size=5)
ppt_t_a<-ggplot(data=FIA_abs,aes(x=PPT_yr_norm,y=T_yr_norm))+
  geom_point()+
  mytheme+labs(x="MAP",y="MAT")+
  annotate("label", x = 1000, y = max(FIA_abs$T_yr_norm), 
           label = paste("Covariance=",round(abs_cov[1,2],2),"\nCorrelation =",round(abs_cor[1,2],2)),
           hjust = 0, vjust = 1, colour="#d95f02",size=5)

cov_plot<-plot_grid(ba_ppt_p,ba_t_p,ppt_t_p,ba_ppt_a,ba_t_a,ppt_t_a, nrow=2,ncol=3,
                    labels = c('A', 'B','C','D','E','F'), label_size = 12)

# Save plot
save_plot(file="pred_cov.png",cov_plot,base_height = 8,base_asp = 2)


### Elevation-predictor covariance ----

## Data
load("./Output/elev_models.rda")

## Plots
elev_ba_plot<-ggplot(data=FIA2,aes(x=elev,y=BALIVE))+
  geom_point(alpha=0.5)+
  geom_line(data=elev_env,aes(x=Elevation,y=BA),col="#1b9e77",size=1.25)+
  labs(x="Elevation (ft)", y="Plot basal area", title="a") + mytheme

elev_ppt_plot<-ggplot(data=FIA2,aes(x=elev,y=PPT_yr_norm))+
  geom_point(alpha=0.5)+
  geom_line(data=elev_env,aes(x=Elevation,y=MAP),col="#1b9e77",size=1.25)+
  labs(x="Elevation (ft)", y="MAP", title="b") + mytheme

elev_t_plot<-ggplot(data=FIA2,aes(x=elev,y=T_yr_norm))+
  geom_point(alpha=0.5)+
  geom_line(data=elev_env,aes(x=Elevation,y=MAT),col="#1b9e77",size=1.25)+
  labs(x="Elevation (ft)", y="MAT", title="c") + mytheme

## Save plots
ggsave(file="elev_b.png", plot=elev_ba_plot,width=4,height=3,units="in",dpi=600)
ggsave(file="elev_p.png", plot=elev_ppt_plot,width=4,height=3,units="in",dpi=600)
ggsave(file="elev_t.png", plot=elev_t_plot,width=4,height=3,units="in",dpi=600)

### Vital rate figures ----

## Data
load("./Output/vital_effects_gam.rda")

## Climate-only, no interactions
# Growth
grow_c_leg_plot_h<-ggplot(data=g_binned,aes(x=PPT,y=grow_PPT_c,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
grow_c_leg_plot_v<-ggplot(data=g_binned,aes(x=PPT,y=grow_PPT_c,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

g_c_leg_h<-legend(grow_c_leg_plot_h)
g_c_leg_v<-legend(grow_c_leg_plot_v)

grow_c_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_clim$dia),xmax=min(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$dia),xmin=max(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_c,size=count_dia),alpha=0.7)+
  geom_line(data=grplot_data_clim,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Diameter increment", title="a")+
  theme(legend.position="none")+mytheme

grow_c_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=0.015,ymax=max(grplot_data_clim$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=0.005,ymax=max(grplot_data_clim$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_c,size=count_PPT),alpha=0.7)+
  geom_line(data=grplot_data_clim,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment", title="a")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none")+mytheme
  #theme(plot.title = element_text(vjust = -8))

grow_c_t<-ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_clim$t_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_clim$t_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_c,size=count_T),alpha=0.7)+
  geom_line(data=grplot_data_clim,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment", title="b")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none")+mytheme

# Survival
surv_c_leg_plot_h<-ggplot(data=s_binned,aes(x=PPT,y=mort_PPT_c,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
surv_c_leg_plot_v<-ggplot(data=s_binned,aes(x=PPT,y=mort_PPT_c,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

s_c_leg_h<-legend(surv_c_leg_plot_h)
s_c_leg_v<-legend(surv_c_leg_plot_v)

surv_c_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_clim$dia),xmax=min(survData$PREVDIA),ymin=0.5,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$dia),xmin=max(survData$PREVDIA),ymin=0.5,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_c,size=count_dia),alpha=0.7)+ #,col="#1b9e77"
  geom_line(data=splot_data_clim,aes(x=dia,y=dia_pred),size=1.25,col="#1b9e77")+ 
  labs(x="Previous diameter", y="Survival", title="b")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none")+mytheme

surv_c_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_clim$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_c,size=count_PPT),alpha=0.7)+ #,col="#1b9e77"
  geom_line(data=splot_data_clim,aes(x=ppt,y=ppt_pred),size=1.25,col="#1b9e77")+
  labs(x="30-year precipitation norm", y="Survival", title="c")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none")+mytheme

surv_c_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_clim$t),xmax=min(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$t),xmin=max(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_c,size=count_T),alpha=0.7)+ #,col="#1b9e77"
  geom_line(data=splot_data_clim,aes(x=t,y=t_pred),size=1.25,col="#1b9e77")+
  labs(x="30-year temperature norm", y="Survival", title="d")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none")+mytheme

# Recruit
recr_c_leg_plot_h<-ggplot(data=r_binned,aes(x=PPT,y=recr_PPT_c,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
recr_c_leg_plot_v<-ggplot(data=r_binned,aes(x=PPT,y=recr_PPT_c,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

r_c_leg_h<-legend(recr_c_leg_plot_h)
r_c_leg_v<-legend(recr_c_leg_plot_v)

recr_c_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_clim$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_clim$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_clim$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_clim$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_c,size=count_PPT),alpha=0.7)+
  geom_line(data=rplot_data_clim,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits", title="e")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

recr_c_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_clim$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=max(r_binned$recr_T_c,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_clim$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=max(r_binned$recr_T_c,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_c,size=count_T),alpha=0.7)+
  geom_line(data=rplot_data_clim,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits", title="f")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

## Climate-only, interactions
# Growth
grow_ci_leg_plot_h<-ggplot(data=g_binned,aes(x=PPT,y=grow_PPT_ci,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
grow_ci_leg_plot_v<-ggplot(data=g_binned,aes(x=PPT,y=grow_PPT_ci,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

g_ci_leg_h<-legend(grow_ci_leg_plot_h)
g_ci_leg_v<-legend(grow_ci_leg_plot_v)

grow_ci_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_clim$dia),xmax=min(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$dia),xmin=max(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_ci,size=count_dia),alpha=0.7)+
  geom_line(data=grplot_data_climint,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Diameter increment", title="c")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

grow_ci_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climint$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=0.015,ymax=max(grplot_data_climint$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climint$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=0.005,ymax=max(grplot_data_climint$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_ci,size=count_PPT),alpha=0.7)+
  geom_line(data=grplot_data_climint,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment", title="a")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

grow_ci_t<-ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_clim$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_climint$t_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_clim$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_climint$t_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_ci,size=count_T),alpha=0.7)+
  geom_line(data=grplot_data_climint,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment", title="b")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

# Survival
surv_ci_leg_plot_h<-ggplot(data=s_binned,aes(x=PPT,y=mort_PPT_ci,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
surv_ci_leg_plot_v<-ggplot(data=s_binned,aes(x=PPT,y=mort_PPT_ci,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

s_ci_leg_h<-legend(surv_ci_leg_plot_h)
s_ci_leg_v<-legend(surv_ci_leg_plot_v)

surv_ci_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_clim$dia),xmax=min(survData$PREVDIA),ymin=0.4,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$dia),xmin=max(survData$PREVDIA),ymin=0.4,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_ci,size=count_dia),alpha=0.7)+ #,col="#1b9e77"
  geom_line(data=splot_data_climint,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Survival", title="d")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

surv_ci_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_clim$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_ci,size=count_PPT),alpha=0.7)+ #,col="#1b9e77"
  geom_line(data=splot_data_climint,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Survival", title="c")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

surv_ci_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_clim$t),xmax=min(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_clim$t),xmin=max(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_ci,size=count_T),alpha=0.7)+ #,col="#1b9e77"
  geom_line(data=splot_data_climint,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Survival", title="d")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

# Recruit
recr_ci_leg_plot_h<-ggplot(data=r_binned,aes(x=PPT,y=recr_PPT_ci,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
recr_ci_leg_plot_v<-ggplot(data=r_binned,aes(x=PPT,y=recr_PPT_ci,size=count_PPT))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

r_ci_leg_h<-legend(recr_ci_leg_plot_h)
r_ci_leg_v<-legend(recr_ci_leg_plot_v)

recr_ci_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climint$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_climint$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climint$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_climint$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_ci,size=count_PPT),alpha=0.7)+
  geom_line(data=rplot_data_climint,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_climint,aes(x=ppt,y=ppt_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits", title="e")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

recr_ci_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climint$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=max(r_binned$recr_T_ci,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climint$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=max(r_binned$recr_T_ci,na.rm=T)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_ci,size=count_T),alpha=0.7)+
  geom_line(data=rplot_data_climint,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_climint,aes(x=t,y=t_pred_c),col="#1b9e77",linetype="dotted",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits", title="f")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

## Climate + competition, no interactions
# Growth
grow_cc_leg_plot_h<-ggplot(data=g_binned,aes(x=BALIVE,y=grow_ba_cc,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
grow_cc_leg_plot_v<-ggplot(data=g_binned,aes(x=BALIVE,y=grow_ba_cc,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

g_cc_leg_h<-legend(grow_cc_leg_plot_h)
g_cc_leg_v<-legend(grow_cc_leg_plot_v)

grow_cc_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$dia),xmax=min(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp$dia),xmin=max(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_cc,size=count_dia),alpha=0.7)+
  geom_line(data=grplot_data_climcomp,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Diameter increment", title="c")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

grow_cc_b <- ggplot(data=grdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$ba),xmax=min(grdata$BALIVE),
                ymin=-0.002,ymax=0.1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=BALIVE,y=grow_ba_cc,size=count_ba),alpha=0.7)+
  geom_line(data=grplot_data_climcomp,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  labs(x="Live basal area", y="Diameter increment", title="a")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

grow_cc_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=min(grplot_data_climcomp$ppt_pred),ymax=max(grplot_data_climcomp$ppt_pred)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=min(grplot_data_climcomp$ppt_pred),ymax=max(grplot_data_climcomp$ppt_pred)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_cc,size=count_PPT),alpha=0.7)+
  geom_line(data=grplot_data_climcomp,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Diameter increment", title="b")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

grow_cc_t <- ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_climcomp$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_climcomp$t_pred)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_climcomp$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_climcomp$t_pred)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_cc,size=count_T),alpha=0.7)+
  geom_line(data=grplot_data_climcomp,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Diameter increment", title="c")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

# Survival
surv_cc_leg_plot_h<-ggplot(data=s_binned,aes(x=BALIVE,y=mort_ba_cc,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
surv_cc_leg_plot_v<-ggplot(data=s_binned,aes(x=BALIVE,y=mort_ba_cc,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

s_cc_leg_h<-legend(surv_cc_leg_plot_h)
s_cc_leg_v<-legend(surv_cc_leg_plot_v)

surv_cc_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_climcomp$dia),xmax=min(survData$PREVDIA),
                ymin=0.5,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp$dia),xmin=max(survData$PREVDIA),
                ymin=0.5,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_cc,size=count_dia),alpha=0.7)+
  geom_line(data=splot_data_climcomp,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Survival", title="d")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

surv_cc_b <- ggplot(data=survData,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(splot_data_climcomp$ba),xmax=min(survData$BALIVE),ymin=0.2,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=BALIVE,y=mort_ba_cc,size=count_ba),alpha=0.7)+
  geom_line(data=splot_data_climcomp,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  labs(x="Live basal area", y="Survival", title="d")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

surv_cc_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_climcomp$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_cc,size=count_PPT),alpha=0.7)+
  geom_line(data=splot_data_climcomp,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Survival", title="e")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

surv_cc_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_climcomp$t),xmax=min(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_climcomp$t),xmin=max(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_cc,size=count_T),alpha=0.7)+
  geom_line(data=splot_data_climcomp,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Survival", title="f")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

# Recruit
recr_cc_leg_plot_h<-ggplot(data=r_binned,aes(x=BALIVE,y=recr_ba_cc,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
recr_cc_leg_plot_v<-ggplot(data=r_binned,aes(x=BALIVE,y=recr_ba_cc,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

r_cc_leg_h<-legend(recr_cc_leg_plot_h)
r_cc_leg_v<-legend(recr_cc_leg_plot_v)

recr_cc_b <- ggplot(data=rdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(rplot_data_climcomp$ba),xmax=min(rdata$BALIVE),
                ymin=-0.05,ymax=2),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=BALIVE,y=recr_ba_cc,size=count_ba),alpha=0.7)+
  geom_line(data=rplot_data_climcomp,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  labs(x="Live basal area", y="Number recruits", title="g")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

recr_cc_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climcomp$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_climcomp$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climcomp$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_climcomp$ppt_pred)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_cc,size=count_PPT),alpha=0.7)+
  geom_line(data=rplot_data_climcomp,aes(x=ppt,y=ppt_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year precipitation norm", y="Number recruits", title="h")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

recr_cc_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_climcomp$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_climcomp$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_cc,size=count_T),alpha=0.7)+
  geom_line(data=rplot_data_climcomp,aes(x=t,y=t_pred),col="#1b9e77",size=1.25)+
  labs(x="30-year temperature norm", y="Number recruits", title="i")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="none") + mytheme

## Climate + competition, interactions
# Growth
grow_i_leg_plot_h<-ggplot(data=g_binned,aes(x=BALIVE,y=grow_ba_i,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
grow_i_leg_plot_v<-ggplot(data=g_binned,aes(x=BALIVE,y=grow_ba_i,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

g_i_leg_h<-legend(grow_i_leg_plot_h)
g_i_leg_v<-legend(grow_i_leg_plot_v)

grow_i_d <- ggplot(data=grdata,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(grplot_data_int$dia),xmax=min(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int$dia),xmin=max(grdata$PREVDIA),ymin=-0.3,ymax=0.1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PREVDIA,y=grow_dia_i,size=count_dia),alpha=0.5)+
  geom_line(data=grplot_data_int,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Diameter increment", title="e")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top",legend.box.margin=margin(-10,-10,-10,-10)) + mytheme

grow_i_b <- ggplot(data=grdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(grplot_data_int$ba),xmax=min(grdata$BALIVE),
                ymin=-0.002,ymax=0.1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=BALIVE,y=grow_ba_i,size=count_ba),alpha=0.5)+
  geom_line(data=grplot_data_int,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  labs(x="Live basal area", y="Diameter increment", title="a")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top",legend.box.margin=margin(-10,-10,-10,-10)) + mytheme

grow_i_p <- ggplot(data=grdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_int$ppt),xmax=min(grdata$PPT_yr_norm),
                ymin=min(grplot_data_int$ppt_pred),ymax=max(grplot_data_int$ppt_pred_75)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int$ppt),xmin=max(grdata$PPT_yr_norm),
                ymin=min(grplot_data_int$ppt_pred),ymax=max(grplot_data_int$ppt_pred_75)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=PPT,y=grow_PPT_i,size=count_PPT),alpha=0.5)+
  geom_line(data=grplot_data_int,aes(x=ppt,y=ppt_pred,linetype="50"),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=ppt,y=ppt_pred_25,linetype="25"),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=ppt,y=ppt_pred_75,linetype="75"),col="#1b9e77",size=1.25)+
  scale_linetype_manual(breaks=c("25","50","75"),values=c("25"="dashed","50"="solid","75"="dotted"),
                        labels=c("25th","50th","75th"), name = "MAT quantile")+
  labs(x="30-year precipitation norm", y="Diameter increment", title="b")+
  guides(size=guide_legend(title="Count",order=1)) +
  theme(legend.position="top",legend.box = "vertical", legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(-10,-10,-10,-10), legend.key.width = unit(1.25,"cm")) + mytheme

grow_i_t <- ggplot(data=grdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(grplot_data_int$t),xmax=min(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_int$t_pred_75)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(grplot_data_int$t),xmin=max(grdata$T_yr_norm),
                ymin=-0.002,ymax=max(grplot_data_int$t_pred_75)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=g_binned,aes(x=T,y=grow_T_i,size=count_T),alpha=0.5)+
  geom_line(data=grplot_data_int,aes(x=t,y=t_pred,linetype="50"),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=t,y=t_pred_25,linetype="25"),col="#1b9e77",size=1.25)+
  geom_line(data=grplot_data_int,aes(x=t,y=t_pred_75,linetype="75"),col="#1b9e77",size=1.25)+
  scale_linetype_manual(breaks=c("25","50","75"),values=c("25"="dashed","50"="solid","75"="dotted"),
                        labels=c("25th","50th","75th"), name = "MAP quantile")+
  labs(x="30-year temperature norm", y="Diameter increment", title="c")+
  guides(size=guide_legend(title="Count",order=1)) +
  theme(legend.position="top",legend.box = "vertical", legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(-10,-10,-10,-10), legend.key.width = unit(1.25,"cm")) + mytheme

# Survival
surv_i_leg_plot_h<-ggplot(data=s_binned,aes(x=BALIVE,y=mort_ba_i,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
surv_i_leg_plot_v<-ggplot(data=s_binned,aes(x=BALIVE,y=mort_ba_i,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

s_i_leg_h<-legend(surv_i_leg_plot_h)
s_i_leg_v<-legend(surv_i_leg_plot_v)

surv_i_d <- ggplot(data=survData,aes(x=PREVDIA))+
  geom_rect(aes(xmin=min(splot_data_int$dia),xmax=min(survData$PREVDIA),
                ymin=0.5,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int$dia),xmin=max(survData$PREVDIA),
                ymin=0.5,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PREVDIA,y=mort_dia_i,size=count_dia),alpha=0.5)+ #,col="#1b9e77"
  geom_line(data=splot_data_int,aes(x=dia,y=dia_pred),col="#1b9e77",size=1.25)+
  labs(x="Previous diameter", y="Survival", title="f")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top",legend.box.margin=margin(-10,-10,-10,-10)) + mytheme

surv_i_b <- ggplot(data=survData,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(splot_data_int$ba),xmax=min(survData$BALIVE),ymin=0.3,ymax=1),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=BALIVE,y=mort_ba_i,size=count_ba),alpha=0.5)+ #,col="#1b9e77"
  geom_line(data=splot_data_int,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  labs(x="Live basal area", y="Survival", title="d")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top",legend.box.margin=margin(-10,-10,-10,-10)) + mytheme

surv_i_p <- ggplot(data=survData,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_int$ppt),xmax=min(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int$ppt),xmin=max(survData$PPT_yr_norm),
                ymin=0.6,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=PPT,y=mort_PPT_i,size=count_PPT),alpha=0.5)+ #,col="#1b9e77"
  geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred,linetype="50"),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred_25,linetype="25"),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=ppt,y=ppt_pred_75,linetype="75"),col="#1b9e77",size=1.25)+
  scale_linetype_manual(breaks=c("25","50","75"),values=c("25"="dashed","50"="solid","75"="dotted"),
                        labels=c("25th","50th","75th"), name = "MAT quantile")+
  labs(x="30-year precipitation norm", y="Survival", title="e")+
  guides(size=guide_legend(title="Count",order=1)) +
  theme(legend.position="top",legend.box = "vertical", legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(-10,-10,-10,-10), legend.key.width = unit(1.25,"cm")) + mytheme

surv_i_t <- ggplot(data=survData,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(splot_data_int$t),xmax=min(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(splot_data_int$t),xmin=max(survData$T_yr_norm),
                ymin=0,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=s_binned,aes(x=T,y=mort_T_i,size=count_T),alpha=0.5)+ #,col="#1b9e77"
  geom_line(data=splot_data_int,aes(x=t,y=t_pred,linetype="50"),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=t,y=t_pred_25,linetype="25"),col="#1b9e77",size=1.25)+
  geom_line(data=splot_data_int,aes(x=t,y=t_pred_75,linetype="75"),col="#1b9e77",size=1.25)+
  scale_linetype_manual(breaks=c("25","50","75"),values=c("25"="dashed","50"="solid","75"="dotted"),
                        labels=c("25th","50th","75th"), name = "MAP quantile")+
  labs(x="30-year temperature norm", y="Survival", title="f")+
  guides(size=guide_legend(title="Count",order=1)) +
  theme(legend.position="top",legend.box = "vertical", legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(-10,-10,-10,-10), legend.key.width = unit(1.25,"cm")) + mytheme

# Recruit
recr_i_leg_plot_h<-ggplot(data=r_binned,aes(x=BALIVE,y=recr_ba_i,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + theme(legend.position="top") + mytheme
recr_i_leg_plot_v<-ggplot(data=r_binned,aes(x=BALIVE,y=recr_ba_i,size=count_ba))+
  geom_point()+
  guides(size=guide_legend(title="Count")) + mytheme

r_i_leg_h<-legend(recr_i_leg_plot_h)
r_i_leg_v<-legend(recr_i_leg_plot_v)

recr_i_b <- ggplot(data=rdata,aes(x=BALIVE))+
  geom_rect(aes(xmin=min(rplot_data_int$ba),xmax=min(rdata$BALIVE),
                ymin=-0.4,ymax=3),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=BALIVE,y=recr_ba_i,size=count_ba),alpha=0.5)+
  geom_line(data=rplot_data_int,aes(x=ba,y=ba_pred),col="#1b9e77",size=1.25)+
  labs(x="Live basal area", y="Number recruits", title="g")+
  guides(size=guide_legend(title="Count")) +
  theme(legend.position="top",legend.box.margin=margin(-10,-10,-10,-10)) + mytheme

recr_i_p <- ggplot(data=rdata,aes(x=PPT_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_int$ppt),xmax=min(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_int$ppt_pred_75)),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_int$ppt),xmin=max(rdata$PPT_yr_norm),
                ymin=-0.05,ymax=max(rplot_data_int$ppt_pred_75)),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=PPT,y=recr_PPT_i,size=count_PPT),alpha=0.5)+
  geom_line(data=rplot_data_int,aes(x=ppt,y=ppt_pred,linetype="50"),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=ppt,y=ppt_pred_25,linetype="25"),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=ppt,y=ppt_pred_75,linetype="75"),col="#1b9e77",size=1.25)+
  scale_linetype_manual(breaks=c("25","50","75"),values=c("25"="dashed","50"="solid","75"="dotted"),
                        labels=c("25th","50th","75th"), name = "MAT quantile")+
  labs(x="30-year precipitation norm", y="Number recruits", title="h")+
  guides(size=guide_legend(title="Count",order=1)) +
  theme(legend.position="top",legend.box = "vertical", legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(-10,-10,-10,-10), legend.key.width = unit(1.25,"cm")) + mytheme

recr_i_t <- ggplot(data=rdata,aes(x=T_yr_norm))+
  geom_rect(aes(xmin=min(rplot_data_int$t),xmax=min(rdata$T_yr_norm),
                ymin=-0.05,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmax=max(rplot_data_int$t),xmin=max(rdata$T_yr_norm),
                ymin=-0.05,ymax=1),fill="grey80",col="grey80",alpha=0.1)+
  geom_point(data=r_binned,aes(x=T,y=recr_T_i,size=count_T),alpha=0.5)+
  geom_line(data=rplot_data_int,aes(x=t,y=t_pred,linetype="50"),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=t,y=t_pred_25,linetype="25"),col="#1b9e77",size=1.25)+
  geom_line(data=rplot_data_int,aes(x=t,y=t_pred_75,linetype="75"),col="#1b9e77",size=1.25)+
  scale_linetype_manual(breaks=c("25","50","75"),values=c("25"="dashed","50"="solid","75"="dotted"),
                        labels=c("25th","50th","75th"), name = "MAP quantile")+
  labs(x="30-year temperature norm", y="Number recruits", title="i")+
  guides(size=guide_legend(title="Count",order=1), linetype=guide_legend()) +
  theme(legend.position="top",legend.box = "vertical", legend.margin=margin(-5,-5,-5,-5),
        legend.box.margin=margin(-10,-10,-10,-10), legend.key.width = unit(1.25,"cm")) + mytheme

## Save plots

# Climate-only, no interactions
ggsave(file="gam_PIED_manuscript_grow_c_d.png", plot=grow_c_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_c_p.png", plot=grow_c_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_c_t.png", plot=grow_c_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="grow_c_leg_h.png", plot=g_c_leg_h)
ggsave(file="grow_c_leg_v.png", plot=g_c_leg_v)

ggsave(file="gam_PIED_manuscript_surv_c_d.png", plot=surv_c_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_c_p.png", plot=surv_c_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_c_t.png", plot=surv_c_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="surv_c_leg_h.png", plot=s_c_leg_h)
ggsave(file="surv_c_leg_v.png", plot=s_c_leg_v)

ggsave(file="gam_PIED_manuscript_recr_c_p.png", plot=recr_c_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_recr_c_t.png", plot=recr_c_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="recr_c_leg_h.png", plot=r_c_leg_h)
ggsave(file="recr_c_leg_v.png", plot=r_c_leg_v)

# Climate-only, interactions
ggsave(file="gam_PIED_manuscript_grow_ci_d.png", plot=grow_ci_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_ci_p.png", plot=grow_ci_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_ci_t.png", plot=grow_ci_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="grow_ci_leg_h.png", plot=g_ci_leg_h)
ggsave(file="grow_ci_leg_v.png", plot=g_ci_leg_v)

ggsave(file="gam_PIED_manuscript_surv_ci_d.png", plot=surv_ci_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_ci_p.png", plot=surv_ci_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_ci_t.png", plot=surv_ci_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="surv_ci_leg_h.png", plot=s_ci_leg_h)
ggsave(file="surv_ci_leg_v.png", plot=s_ci_leg_v)

ggsave(file="gam_PIED_manuscript_recr_ci_p.png", plot=recr_ci_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_recr_ci_t.png", plot=recr_ci_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="recr_ci_leg_h.png", plot=r_ci_leg_h)
ggsave(file="recr_ci_leg_v.png", plot=r_ci_leg_v)

# Climate + competition, no interactions
ggsave(file="gam_PIED_manuscript_grow_cc_d.png", plot=grow_cc_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_cc_b.png", plot=grow_cc_b,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_cc_p.png", plot=grow_cc_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_cc_t.png", plot=grow_cc_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="grow_cc_leg_h.png", plot=g_cc_leg_h)
ggsave(file="grow_cc_leg_v.png", plot=g_cc_leg_v)

ggsave(file="gam_PIED_manuscript_surv_cc_d.png", plot=surv_cc_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_cc_b.png", plot=surv_cc_b,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_cc_p.png", plot=surv_cc_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_cc_t.png", plot=surv_cc_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="surv_cc_leg_h.png", plot=s_cc_leg_h)
ggsave(file="surv_cc_leg_v.png", plot=s_cc_leg_v)

ggsave(file="gam_PIED_manuscript_recr_cc_b.png", plot=recr_cc_b,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_recr_cc_p.png", plot=recr_cc_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_recr_cc_t.png", plot=recr_cc_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="recr_cc_leg_h.png", plot=r_cc_leg_h)
ggsave(file="recr_cc_leg_v.png", plot=r_cc_leg_v)

# Climate + competition, interactions
ggsave(file="gam_PIED_manuscript_grow_i_d.png", plot=grow_i_d,
       width=5,height=3.75,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_i_b.png", plot=grow_i_b,
       width=5,height=3.75,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_i_p.png", plot=grow_i_p,
       width=5,height=3.75,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_grow_i_t.png", plot=grow_i_t,
       width=5,height=3.75,units="in",dpi=600)
ggsave(file="grow_i_leg_h.png", plot=g_i_leg_h)
ggsave(file="grow_i_leg_v.png", plot=g_i_leg_v)

ggsave(file="gam_PIED_manuscript_surv_i_d.png", plot=surv_i_d,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_i_b.png", plot=surv_i_b,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_i_p.png", plot=surv_i_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_surv_i_t.png", plot=surv_i_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="surv_i_leg_h.png", plot=s_i_leg_h)
ggsave(file="surv_i_leg_v.png", plot=s_i_leg_v)

ggsave(file="gam_PIED_manuscript_recr_i_b.png", plot=recr_i_b,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_recr_i_p.png", plot=recr_i_p,
       width=4,height=3,units="in",dpi=600)
ggsave(file="gam_PIED_manuscript_recr_i_t.png", plot=recr_i_t,
       width=4,height=3,units="in",dpi=600)
ggsave(file="recr_i_leg_h.png", plot=r_i_leg_h)
ggsave(file="recr_i_leg_v.png", plot=r_i_leg_v)

### Lambda figures ----

## Data
load("./Output/lambda_effects_gam.rda")

## Climate-only, no interactions
lambda_c_b <- ggplot(data = predictorDFs_clim[[1]], aes(x = BALIVE, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_clim[[1]]$BALIVE),xmax=min(FIA$BALIVE),
                ymin=min(predictorDFs_clim[[1]]$lambda),ymax=max(predictorDFs_clim[[1]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("Plot basal area") + ylab("Lambda") + mytheme

lambda_c_p <- ggplot(data = predictorDFs_clim[[2]], aes(x = PPT_yr, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_clim[[2]]$PPT_yr),xmax=min(FIA$PPT_yr),
                ymin=min(predictorDFs_clim[[2]]$lambda),ymax=max(predictorDFs_clim[[2]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max(FIA$PPT_yr),xmax=max(predictorDFs_clim[[2]]$PPT_yr),
                ymin=min(predictorDFs_clim[[2]]$lambda),ymax=max(predictorDFs_clim[[2]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("MAP") + ylab("Lambda") + mytheme

lambda_c_t <- ggplot(data = predictorDFs_clim[[3]], aes(x = T_yr, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_clim[[3]]$T_yr),xmax=min(FIA$T_yr),
                ymin=min(predictorDFs_clim[[3]]$lambda),ymax=max(predictorDFs_clim[[3]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max(FIA$T_yr),xmax=max(predictorDFs_clim[[3]]$T_yr),
                ymin=min(predictorDFs_clim[[3]]$lambda),ymax=max(predictorDFs_clim[[3]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("MAT") + ylab("Lambda") + mytheme

## Climate + competition, no interactions
lambda_cc_b <- ggplot(data = predictorDFs_climcomp[[1]], aes(x = BALIVE, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_climcomp[[1]]$BALIVE),xmax=min(FIA$BALIVE),
                ymin=min(predictorDFs_climcomp[[1]]$lambda),ymax=max(predictorDFs_climcomp[[1]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("Plot basal area") + ylab("Lambda") + mytheme

lambda__cc_p <- ggplot(data = predictorDFs_climcomp[[2]], aes(x = PPT_yr, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_climcomp[[2]]$PPT_yr),xmax=min(FIA$PPT_yr),
                ymin=min(predictorDFs_climcomp[[2]]$lambda),ymax=max(predictorDFs_climcomp[[2]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max(FIA$PPT_yr),xmax=max(predictorDFs_climcomp[[2]]$PPT_yr),
                ymin=min(predictorDFs_climcomp[[2]]$lambda),ymax=max(predictorDFs_climcomp[[2]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("MAP") + ylab("Lambda") + mytheme

lambda_cc_t <- ggplot(data = predictorDFs_climcomp[[3]], aes(x = T_yr, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_climcomp[[3]]$T_yr),xmax=min(FIA$T_yr),
                ymin=min(predictorDFs_climcomp[[3]]$lambda),ymax=max(predictorDFs_climcomp[[3]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max(FIA$T_yr),xmax=max(predictorDFs_climcomp[[3]]$T_yr),
                ymin=min(predictorDFs_climcomp[[3]]$lambda),ymax=max(predictorDFs_climcomp[[3]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("MAT") + ylab("Lambda") + mytheme

## Climate + competition, interactions
lambda_i_b <- ggplot(data = predictorDFs_int[[1]], aes(x = BALIVE, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_int[[1]]$BALIVE),xmax=min(FIA$BALIVE),
                ymin=min(predictorDFs_int[[1]]$lambda),ymax=max(predictorDFs_int[[1]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("Plot basal area") + ylab("Lambda") + mytheme

lambda_i_p <- ggplot(data = predictorDFs_int[[2]], aes(x = PPT_yr, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_int[[2]]$PPT_yr),xmax=min(FIA$PPT_yr),
                ymin=min(predictorDFs_int[[2]]$lambda),ymax=max(predictorDFs_int[[2]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max(FIA$PPT_yr),xmax=max(predictorDFs_int[[2]]$PPT_yr),
                ymin=min(predictorDFs_int[[2]]$lambda),ymax=max(predictorDFs_int[[2]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("MAP") + ylab("Lambda") + mytheme

lambda_i_t <- ggplot(data = predictorDFs_int[[3]], aes(x = T_yr, y = lambda)) + 
  geom_rect(aes(xmin=min(predictorDFs_int[[3]]$T_yr),xmax=min(FIA$T_yr),
                ymin=min(predictorDFs_int[[3]]$lambda),ymax=max(predictorDFs_int[[3]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max(FIA$T_yr),xmax=max(predictorDFs_int[[3]]$T_yr),
                ymin=min(predictorDFs_int[[3]]$lambda),ymax=max(predictorDFs_int[[3]]$lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(col="#1b9e77",size=1.25) +
  xlab("MAT") + ylab("Lambda") + mytheme

## Save plots
ggsave(file="gam_lam_c_b.png", plot=lambda_c_b,width=4,height=3,units="in",dpi=600)
ggsave(file="gam_lam_c_p.png", plot=lambda_c_p,width=4,height=3,units="in",dpi=600)
ggsave(file="gam_lam_c_t.png", plot=lambda_c_t,width=4,height=3,units="in",dpi=600)

ggsave(file="gam_lam_cc_b.png", plot=lambda_cc_b,width=4,height=3,units="in",dpi=600)
ggsave(file="gam_lam_cc_p.png", plot=lambda_cc_p,width=4,height=3,units="in",dpi=600)
ggsave(file="gam_lam_cc_t.png", plot=lambda_cc_t,width=4,height=3,units="in",dpi=600)

ggsave(file="gam_lam_i_b.png", plot=lambda_i_b,width=4,height=3,units="in",dpi=600)
ggsave(file="gam_lam_i_p.png", plot=lambda_i_p,width=4,height=3,units="in",dpi=600)
ggsave(file="gam_lam_i_t.png", plot=lambda_i_t,width=4,height=3,units="in",dpi=600)

all <- plot_grid(A, B, C, labels = c("a", "b", "c"), align = "hv")
save_plot("./Output/Lambda_plots_int_q.png", all, base_aspect_ratio = 2)

### Vital-lambda combo plot
legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

vital_lambda_plot_c<-plot_grid(grow_c_p,grow_c_t,surv_c_p,surv_c_t,recr_c_p,recr_c_t,
                             nrow=4,ncol=2,
                    labels = c('A','B','C','D','E','F','G','H'), label_size = 12)

### Maps ----

## Data
lambda_c<-raster("./Output/tifs/PIED.clim_lambda_gam.tif")
lambda_ccl<-raster("./Output/tifs/PIED.climclamp_lambda.tif")
lambda_cc<-raster("./Output/tifs/PIED.climcomp_lambda_gam.tif")
lambda_ccf<-raster("./Output/tifs/PIED.climcompfire_lambda.tif")
lambda_i<-raster("./Output/tifs/PIED.int_lambda_gam.tif")

growth_c<-raster("./Output/tifs/PIED.clim_growth_gam.tif")
growth_ccl<-raster("./Output/tifs/PIED.climclamp_growth.tif")
growth_cc<-raster("./Output/tifs/PIED.climcomp_growth_gam.tif")
growth_ccf<-raster("./Output/tifs/PIED.climcompfire_growth.tif")
growth_i<-raster("./Output/tifs/PIED.int_growth_gam.tif")

survival_c<-raster("./Output/tifs/PIED.clim_survival_gam.tif")
survival_ccl<-raster("./Output/tifs/PIED.climclamp_survival.tif")
survival_cc<-raster("./Output/tifs/PIED.climcomp_survival_gam.tif")
survival_ccf<-raster("./Output/tifs/PIED.climcompfire_survival.tif")
survival_i<-raster("./Output/tifs/PIED.int_survival_gam.tif")

reproduction_c<-raster("./Output/tifs/PIED.clim_reproduction_gam.tif")
reproduction_ccl<-raster("./Output/tifs/PIED.climclamp_reproduction.tif")
reproduction_cc<-raster("./Output/tifs/PIED.climcomp_reproduction_gam.tif")
reproduction_ccf<-raster("./Output/tifs/PIED.climcompfire_reproduction.tif")
reproduction_i<-raster("./Output/tifs/PIED.int_reproduction_gam.tif")

extrap <- raster("./Output/tifs/extrap.tif")

## Color palettes
pal_seq <- colorRampPalette(brewer.pal(n=9, name = "PuBuGn"))
pal_div <- colorRampPalette(brewer.pal(n=11, name = "BrBG"))

## Climate raster maps
pdf("./Output/Climate_maps.pdf")
plot(ba_raster, main = "Live Basal Area", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(ppt_yr_raster, main = "MAP", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
plot(t_yr_raster, main = "MAT", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

## Climate-only, no interactions
pdf("./Output/PIED_clim_gam.pdf")
plot(lambda_c, main = "Lambda", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(growth_c, main = "Growth", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(survival_c, main = "Survival",col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(reproduction_c, main = "Reproduction", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
dev.off()

png(file="./Output/PIED_clim_gam_lam.png",4,4,units="in",type="cairo",res=600)
plot(lambda_c, col=pal_seq(50), xlab="Longitude", ylab="Latitude",cex.axis=0.8,cex.lab=0.9,
     legend.args=list(text='Lambda',adj=0,cex=0.9),
     axis.args=list(cex.axis=0.8)) 
plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

cuts<-c(seq(min(na.omit(values(lambda_c))),1,length=31),
        seq(1,max(na.omit(values(lambda_c))),length=31)[2:31])

cuts_round<-round(cuts,2)
legend_text<-(c(cuts_round[1],cuts_round[11],cuts_round[21],cuts_round[31],cuts_round[41],cuts_round[51],cuts_round[61]))

png(file="./Output/PIED_clim_gam_lam_div.png",4,4,units="in",type="cairo",res=600)
plot(lambda_c, breaks=cuts, col=pal_div(61), xlab="Longitude", ylab="Latitude",cex.axis=0.8,cex.lab=0.9,
     axis.args=list(cex.axis=0.8),legend=F) 
color.legend(-100,33,-99,40.5,
             legend=legend_text,rect.col=pal_div(61),align="rb",gradient="y",cex=0.9)
plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
legend(x=-102,y=43,legend=expression(lambda),bty = "n",cex=1.5)
dev.off()

# Climate + competition, no interactions
pdf("./Output/PIED_climcomp_gam.pdf")
plot(lambda_cc, main = "Lambda", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(growth_cc, main = "Growth", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(survival_cc, main = "Survival",col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(reproduction_cc, main = "Reproduction", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
dev.off()

png(file="./Output/PIED_climcomp_gam_lam.png",4,4,units="in",type="cairo",res=600)
plot(lambda_cc, col=pal_seq(50), xlab="Longitude", ylab="Latitude",cex.axis=0.8,cex.lab=0.9,
     legend.args=list(text='Lambda',adj=0,cex=0.9),
     axis.args=list(cex.axis=0.8))
plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

cuts<-c(seq(min(na.omit(values(lambda_cc))),1,length=31),
        seq(1,max(na.omit(values(lambda_cc))),length=31)[2:31])

cuts_round<-round(cuts,2)
legend_text<-(c(cuts_round[1],cuts_round[11],cuts_round[21],cuts_round[31],cuts_round[41],cuts_round[51],cuts_round[61]))

png(file="./Output/PIED_climcomp_gam_lam_div.png",4,4,units="in",type="cairo",res=600)
plot(lambda_cc, breaks=cuts, col=pal_div(61), xlab="Longitude", ylab="Latitude",cex.axis=0.8,cex.lab=0.9,
     axis.args=list(cex.axis=0.8),legend=F)
color.legend(-100,33,-99,40.5,
             legend=legend_text,rect.col=pal_div(61),align="rb",gradient="y",cex=0.9)
plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
legend(x=-102,y=43,legend=expression(lambda),bty = "n",cex=1.5)
dev.off()

# Climate + competition, interactions
pdf("./Output/PIED_int_gam.pdf")
plot(lambda_i, main = "Lambda", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(growth_i, main = "Growth", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(survival_i, main = "Survival",col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
plot(reproduction_i, main = "Reproduction", col=pal_seq(50)); points(LAT ~ LON, FIA, pch = 19, cex = 0.05); plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
dev.off()

png(file="./Output/PIED_int_gam_lam.png",4,4,units="in",type="cairo",res=600)
plot(lambda_i, col=pal_seq(50), xlab="Longitude", ylab="Latitude",cex.axis=0.8,cex.lab=0.9,
     legend.args=list(text='Lambda',adj=0,cex=0.9),
     axis.args=list(cex.axis=0.8)) 
plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
dev.off()

cuts<-c(seq(min(na.omit(values(lambda_i))),1,length=31),
        seq(1,max(na.omit(values(lambda_i))),length=31)[2:31])

cuts_round<-round(cuts,2)
legend_text<-(c(cuts_round[1],cuts_round[11],cuts_round[21],cuts_round[31],cuts_round[41],cuts_round[51],cuts_round[61]))

png(file="./Output/PIED_int_gam_lam_div.png",4,4,units="in",type="cairo",res=600)
plot(lambda_i, breaks=cuts, col=pal(61), xlab="Longitude", ylab="Latitude",cex.axis=0.8,cex.lab=0.9,
     axis.args=list(cex.axis=0.8),legend=F) 
color.legend(-100,33,-99,40.5,
             legend=legend_text,rect.col=pal(61),align="rb",gradient="y",cex=0.9)
plot(extrap,col=grey(0.1),legend=F,alpha=0.2,add=T)
points(LAT ~ LON, FIA, pch = 19, cex = 0.05)
legend(x=-102,y=43,legend=expression(lambda),bty = "n",cex=1.5)
dev.off()


### Presence-absence ----

## Data
load("./Output/residual_gam.rda")

## Lambda vs occurrence
pres_plot_c<-ggplot(data=subset(pres_binned,model=="c"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  annotate("label", x = 0.89, y = 0.5, 
           label = paste("Deviance=",round(pa_c$deviance,2),
                         "\nAIC =",round(pa_c$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "lambda", y = "Probability of occurrence")

pres_plot_cc<-ggplot(data=subset(pres_binned,model=="cc"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  #stat_function(aes(x=seq(1,1.9,length=50)),fun=test,col="red")+
  annotate("label", x = 0.88, y = 0.55, 
           label = paste("Deviance=",round(pa_cc$deviance,2),
                         "\nAIC =",round(pa_cc$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "lambda", y = "Probability of occurrence")

pres_plot_i<-ggplot(data=subset(pres_binned,model=="i"),aes(x=lam,y=pres))+
  geom_point(aes(size=count_lam))+
  geom_line(aes(y=pred),size=1)+
  #stat_function(aes(x=seq(1,1.5,length=50)),fun=test,col="red")+
  annotate("label", x = 0.82, y = 0.55, 
           label = paste("Deviance=",round(pa_i$deviance,2),
                         "\nAIC =",round(pa_i$aic,2)),
           hjust = 0, vjust = 1, size=5)+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "lambda", y = "Probability of occurrence")

# Save plots
save_plot(file="pres_plot_c_gam.png",pres_plot_c)
save_plot(file="pres_plot_cc_gam.png",pres_plot_cc)
save_plot(file="pres_plot_i_gam.png",pres_plot_i)

## Lambda histograms
lambdaHist_c<-ggplot(data=FIA_lambda,aes(x=lambda_c,fill=PApied_f)) +
  geom_histogram(aes(y = ..density..),alpha=0.7)+
  geom_vline(data=l_means, aes(xintercept=lambda_c,  colour=PApied_f),
             linetype="dashed", size=1)+
  scale_fill_manual(values = c("#1b9e77","#d95f02"), limits = c("0","1"), breaks = c("0","1"),
                    name = "PIED presence", labels = c("Absent","Present"))+
  scale_colour_manual(values = c("#1b9e77","#d95f02"), limits = c("0","1"), breaks = c("0","1"),
                      name = "PIED presence", labels = c("Absent","Present"))+
  labs(x="Lambda", y="Density")+
  mytheme

lambdaHist_cc<-ggplot(data=FIA_lambda,aes(x=lambda_cc,fill=PApied_f)) +
  geom_histogram(aes(y = ..density..),alpha=0.7)+
  geom_vline(data=l_means, aes(xintercept=lambda_cc,  colour=PApied_f),
             linetype="dashed", size=1)+
  scale_fill_manual(values = c("#1b9e77","#d95f02"), limits = c("0","1"), breaks = c("0","1"),
                    name = "PIED presence", labels = c("Absent","Present"))+
  scale_colour_manual(values = c("#1b9e77","#d95f02"), limits = c("0","1"), breaks = c("0","1"),
                      name = "PIED presence", labels = c("Absent","Present"))+
  labs(x="Lambda", y="Density")+
  theme_classic()

lambdaHist_i<-ggplot(data=FIA_lambda,aes(x=lambda_i,fill=PApied_f)) +
  geom_histogram(aes(y = ..density..),alpha=0.7)+
  geom_vline(data=l_means, aes(xintercept=lambda_i,  colour=PApied_f),
             linetype="dashed", size=1)+
  scale_fill_manual(values = c("#1b9e77","#d95f02"), limits = c("0","1"), breaks = c("0","1"),
                    name = "PIED presence", labels = c("Absent","Present"))+
  scale_colour_manual(values = c("#1b9e77","#d95f02"), limits = c("0","1"), breaks = c("0","1"),
                      name = "PIED presence", labels = c("Absent","Present"))+
  labs(x="Lambda", y="Density")+
  mytheme

# Save plots
save_plot(file="lam_hist_c.png",lambdaHist_c,base_asp = 1.2)
save_plot(file="lam_hist_cc.png",lambdaHist_cc,base_asp = 1.2)
save_plot(file="lam_hist_i.png",lambdaHist_i,base_asp = 1.2)

## Niche space plots
# Raw lambda
niche_c<-ggplot(data=interpdf_c,aes(x=x,y=y))+
  geom_tile(aes(fill=z))+
  geom_density2d(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),color="black")+
  geom_encircle(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),col="black",expand=0.01)+
  scale_fill_distiller(palette="PuBuGn",direction=1,name="Lambda")+
  labs(x="MAP",y="MAT")

niche_cc<-ggplot(data=interpdf_cc,aes(x=x,y=y))+
  geom_tile(aes(fill=z))+
  geom_density2d(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),color="black")+
  geom_encircle(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),col="black",expand=0.01)+
  scale_fill_distiller(palette="PuBuGn",direction=1,name="Lambda")+
  labs(x="MAP",y="MAT")

niche_i<-ggplot(data=interpdf_i,aes(x=x,y=y))+
  geom_tile(aes(fill=z))+
  geom_density2d(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),color="black")+
  geom_encircle(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),col="black",expand=0.01)+
  scale_fill_distiller(palette="PuBuGn",direction=1,name="Lambda")+
  labs(x="MAP",y="MAT")

# Deviation from mean lambda
niche_c_diff<-ggplot(data=interpdf_c,aes(x=x,y=y))+
  geom_tile(aes(fill=z_diff))+
  geom_density2d(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),color="black")+
  geom_encircle(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),col="black",expand=0.01)+
  scale_fill_distiller(palette="BrBG",direction=1,values=rescale(q_c),
                       name="Deviation \nin lambda")+
  labs(x="MAP",y="MAT")+mytheme

niche_cc_diff<-ggplot(data=interpdf_cc,aes(x=x,y=y))+
  geom_tile(aes(fill=z_diff))+
  geom_density2d(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),color="black")+
  geom_encircle(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),col="black",expand=0.01)+
  scale_fill_distiller(palette="BrBG",direction=1,values=rescale(q_cc),
                       name="Deviation \nin lambda")+
  labs(x="MAP",y="MAT")+mytheme

niche_i_diff<-ggplot(data=interpdf_i,aes(x=x,y=y))+
  geom_tile(aes(fill=z_diff))+
  geom_density2d(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),color="black")+
  geom_encircle(data=FIA_pied_pres,aes(x=PPT_yr,y=T_yr),col="black",expand=0.01)+
  scale_fill_distiller(palette="BrBG",direction=1,values=rescale(q_i),
                       name="Deviation \nin lambda")+
  labs(x="MAP",y="MAT")+mytheme

# Save plots
save_plot(file="niche_c.png",niche_c,base_asp = 1.5)
save_plot(file="niche_cc.png",niche_cc,base_asp = 1.5)
save_plot(file="niche_i.png",niche_i,base_asp = 1.5)

save_plot(file="niche_c_diff.png",niche_c_diff,base_asp = 1.5)
save_plot(file="niche_cc_diff.png",niche_cc_diff,base_asp = 1.5)
save_plot(file="niche_i_diff.png",niche_i_diff,base_asp = 1.5)


### Residual analysis ----
## Data
load("./Output/residual_gam.rda") # Same data used in presence-absence plots

## Residuals vs predictors

# Scaled residuals
ba_resid_c_plot<-ggplot(data=res_binned,aes(x=ba_binned,y=resid_c_ba,size=count_ba))+
  geom_point()+mytheme+labs(x = "Basal area", y = "Scaled residuals")
ba_resid_ci_plot<-ggplot(data=res_binned,aes(x=ba_binned,y=resid_ci_ba,size=count_ba))+
  geom_point()+mytheme+labs(x = "Basal area", y = "Scaled residuals")
ba_resid_cc_plot<-ggplot(data=res_binned,aes(x=ba_binned,y=resid_cc_ba,size=count_ba))+
  geom_point()+mytheme+labs(x = "Basal area", y = "Scaled residuals")
ba_resid_i_plot<-ggplot(data=res_binned,aes(x=ba_binned,y=resid_i_ba,size=count_ba))+
  geom_point()+mytheme+labs(x = "Basal area", y = "Scaled residuals")

ppt_resid_c_plot<-ggplot(data=res_binned,aes(x=ppt_binned,y=resid_c_ppt,size=count_ppt))+
  geom_point()+mytheme+labs(x = "MAP", y = "Scaled residuals")
ppt_resid_ci_plot<-ggplot(data=res_binned,aes(x=ppt_binned,y=resid_ci_ppt,size=count_ppt))+
  geom_point()+mytheme+labs(x = "MAP", y = "Scaled residuals")
ppt_resid_cc_plot<-ggplot(data=res_binned,aes(x=ppt_binned,y=resid_cc_ppt,size=count_ppt))+
  geom_point()+mytheme+labs(x = "MAP", y = "Scaled residuals")
ppt_resid_i_plot<-ggplot(data=res_binned,aes(x=ppt_binned,y=resid_i_ppt,size=count_ppt))+
  geom_point()+mytheme+labs(x = "MAP", y = "Scaled residuals")

t_resid_c_plot<-ggplot(data=res_binned,aes(x=t_binned,y=resid_c_t,size=count_t))+
  geom_point()+mytheme+labs(x = "MAT", y = "Scaled residuals")
t_resid_ci_plot<-ggplot(data=res_binned,aes(x=t_binned,y=resid_ci_t,size=count_t))+
  geom_point()+mytheme+labs(x = "MAT", y = "Scaled residuals")
t_resid_cc_plot<-ggplot(data=res_binned,aes(x=t_binned,y=resid_cc_t,size=count_t))+
  geom_point()+mytheme+labs(x = "MAT", y = "Scaled residuals")
t_resid_i_plot<-ggplot(data=res_binned,aes(x=t_binned,y=resid_i_t,size=count_t))+
  geom_point()+mytheme+labs(x = "MAT", y = "Scaled residuals")

# Unscaled residuals
ba_resid2_plot<-ggplot(data=res2_binned,aes(x=ba,y=resid_ba,colour=model))+
  geom_abline(intercept=0,slope=0)+
  geom_point(aes(size=count_ba))+
  scale_colour_manual(breaks=c("c","ci","cc","i"), 
                      values=c("c"="#1b9e77","ci"="#d95f02","cc"="#7570b3",
                               "i"="#e7298a"),
                      labels=c("Clim","ClimInt","ClimComp","ClimCompInt"),
                      name="Model")+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "Basal area", y = "Residuals")

ppt_resid2_plot<-ggplot(data=res2_binned,aes(x=ppt,y=resid_ppt,colour=model))+
  geom_abline(intercept=0,slope=0)+
  geom_point(aes(size=count_ppt))+
  scale_colour_manual(breaks=c("c","ci","cc","i"), 
                      values=c("c"="#1b9e77","ci"="#d95f02","cc"="#7570b3",
                               "i"="#e7298a"),
                      labels=c("Clim","ClimInt","ClimComp","ClimCompInt"),
                      name="Model")+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "MAP", y = "Residuals")

t_resid2_plot<-ggplot(data=res2_binned,aes(x=t,y=resid_t,colour=model))+
  geom_abline(intercept=0,slope=0)+
  geom_point(aes(size=count_t))+
  scale_colour_manual(breaks=c("c","ci","cc","i"), 
                      values=c("c"="#1b9e77","ci"="#d95f02","cc"="#7570b3",
                               "i"="#e7298a"),
                      labels=c("Clim","ClimInt","ClimComp","ClimCompInt"),
                      name="Model")+
  guides(size=guide_legend(title="Count")) +
  mytheme+labs(x = "MAT", y = "Residuals")

# Save plots
save_plot(file="ba_resid.png",ba_resid2_plot)
save_plot(file="ppt_resid.png",ppt_resid2_plot)
save_plot(file="t_resid.png",t_resid2_plot)

## Plot spatially
# Climate-only, no interactions
BA_raster_c<-gplot(ba_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_c2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr")+
  scale_colour_distiller(palette="PRGn")

PPT_raster_c<-gplot(ppt_yr_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_c2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr")+
  scale_colour_distiller(palette="PRGn")

T_raster_c<-gplot(t_yr_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_c2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr",direction=1)+
  scale_colour_distiller(palette="PRGn")

# Climate + competition, no interactions
BA_raster_cc<-gplot(ba_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_cc2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr")+
  scale_colour_distiller(palette="PRGn")

PPT_raster_cc<-gplot(ppt_yr_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_cc2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr")+
  scale_colour_distiller(palette="PRGn")

T_raster_cc<-gplot(t_yr_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_cc2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr",direction=1)+
  scale_colour_distiller(palette="PRGn")

# Climate + competition, interactions
BA_raster_i<-gplot(ba_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_i2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr")+
  scale_colour_distiller(palette="PRGn")

PPT_raster_i<-gplot(ppt_yr_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_i2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr")+
  scale_colour_distiller(palette="PRGn")

T_raster_i<-gplot(t_yr_raster)+ geom_tile(aes(fill = value))+
  geom_point(data=FIA_lambda,aes(x=lon,y=lat,col=resid_i2),size=0.8)+
  scale_fill_distiller(palette="YlOrBr",direction=1)+
  scale_colour_distiller(palette="PRGn")

# Save plots
save_plot(file="BA_raster_c.png",BA_raster_c,base_asp = 1.5)
save_plot(file="PPT_raster_c.png",PPT_raster_c,base_asp = 1.5)
save_plot(file="T_raster_c.png",T_raster_c,base_asp = 1.5)

save_plot(file="BA_raster_cc.png",BA_raster_c,base_asp = 1.5)
save_plot(file="PPT_raster_cc.png",PPT_raster_cc,base_asp = 1.5)
save_plot(file="T_raster_cc.png",T_raster_c,base_asp = 1.5)

save_plot(file="BA_raster_i.png",BA_raster_c,base_asp = 1.5)
save_plot(file="PPT_raster_i.png",PPT_raster_i,base_asp = 1.5)
save_plot(file="T_raster_i.png",T_raster_c,base_asp = 1.5)


### Elevation-lambda plots ----

## Data
lam_elev_data<-read.csv("./Output/lam_elev_data.csv")
load("./Output/elev_limits.rda")

## Plots
lam_elev_plot_c<-ggplot(data=subset(lam_elev_data, Model=="c"),aes(x=Elevation,y=Lambda))+
  geom_rect(aes(xmin=min(lam_elev_data$Elevation),xmax=min_elev_pied,
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(lam_elev_data$Elevation),
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  labs(x = "Elevation (ft)", y = expression(lambda))+
  mytheme

lam_elev_plot_cc<-ggplot(data=subset(lam_elev_data, Model=="cc"),aes(x=Elevation,y=Lambda))+
  geom_rect(aes(xmin=min(lam_elev_data$Elevation),xmax=min_elev_pied,
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(lam_elev_data$Elevation),
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  labs(x = "Elevation (ft)", y = expression(lambda))+
  mytheme

lam_elev_plot_i<-ggplot(data=subset(lam_elev_data, Model=="i"),aes(x=Elevation,y=Lambda))+
  geom_rect(aes(xmin=min(lam_elev_data$Elevation),xmax=min_elev_pied,
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(lam_elev_data$Elevation),
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  labs(x = "Elevation (ft)", y = expression(lambda))+
  mytheme

lam_elev_plot<-ggplot(data=lam_elev_data,aes(x=Elevation,y=Lambda,colour=Model))+
  geom_rect(aes(xmin=min(lam_elev_data$Elevation),xmax=min_elev_pied,
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(lam_elev_data$Elevation),
                ymin=min(lam_elev_data$Lambda),ymax=max(lam_elev_data$Lambda)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("c","cc","i"), 
                      values=c("c"="#1b9e77","cc"="#d95f02","i"="#7570b3"),
                      labels=c("c","cc","i"),
                      name="Model")+
  labs(x = "Elevation (ft)", y = expression(lambda))+
  mytheme

## Save plots
ggsave(file="lam_elev_plot_c.png", plot=lam_elev_plot_c,width=4,height=3,units="in",dpi=600)
ggsave(file="lam_elev_plot_cc.png", plot=lam_elev_plot_cc,width=4,height=3,units="in",dpi=600)
ggsave(file="lam_elev_plot_i.png", plot=lam_elev_plot_i,width=4,height=3,units="in",dpi=600)
ggsave(file="lam_elev_plot.png", plot=lam_elev_plot,width=4,height=3,units="in",dpi=600)


### Elasticity analysis ----

## Data
elast_data_vital<-read.csv("./Output/elast_vital.csv")
load("./Output/elev_limits.rda")

## Plots
elast_plot_vital_c<-#ggplot(data=elast_data_vital,aes(x=Elev,y=Elast_c,col=Rate))+
  ggplot(data=subset(elast_data_vital,Rate!="Growth"),aes(x=Elev,y=Elast_c,col=Rate))+
  geom_rect(aes(xmin=min(elast_data_vital$Elev),xmax=min_elev_pied,
                ymin=min(elast_data_vital$Elast_c),ymax=max(elast_data_vital$Elast_c)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(elast_data_vital$Elev),
                ymin=min(elast_data_vital$Elast_c),ymax=max(elast_data_vital$Elast_c)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  #scale_colour_manual(breaks=c("Growth","Survival","Recruitment"), 
  #                    values=c("Growth"="#1b9e77","Survival"="#d95f02","Recruitment"="#7570b3"),
  #                    labels=c("Growth","Survival","Recruitment"),
  #                    name="Vital rate")+
  scale_colour_manual(breaks=c("Survival","Recruitment"), 
                      values=c("Survival"="#1b9e77","Recruitment"="#7570b3"),
                      labels=c("Grow/Surv","Recr"),
                      name="Kernel")+
  labs(x = "Elevation (ft)", y ="Elasticity")+
  mytheme

elast_plot_vital_cc<-#ggplot(data=elast_data_vital,aes(x=Elev,y=Elast_cc,col=Rate))+
  ggplot(data=subset(elast_data_vital,Rate!="Growth"),aes(x=Elev,y=Elast_cc,col=Rate))+
  geom_rect(aes(xmin=min(elast_data_vital$Elev),xmax=min_elev_pied,
                ymin=min(elast_data_vital$Elast_cc),ymax=max(elast_data_vital$Elast_cc)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(elast_data_vital$Elev),
                ymin=min(elast_data_vital$Elast_cc),ymax=max(elast_data_vital$Elast_cc)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  #scale_colour_manual(breaks=c("Growth","Survival","Recruitment"), 
  #                    values=c("Growth"="#1b9e77","Survival"="#d95f02","Recruitment"="#7570b3"),
  #                    labels=c("Growth","Survival","Recruitment"),
  #                    name="Vital rate")+
  scale_colour_manual(breaks=c("Survival","Recruitment"), 
                      values=c("Survival"="#1b9e77","Recruitment"="#7570b3"),
                      labels=c("Grow/Surv","Recr"),
                      name="Kernel")+
  labs(x = "Elevation (ft)", y = "Elasticity")+
  mytheme

elast_plot_vital_i<-#ggplot(data=elast_data_vital,aes(x=Elev,y=Elast_i,col=Rate))+
  ggplot(data=subset(elast_data_vital,Rate!="Growth"),aes(x=Elev,y=Elast_i,col=Rate))+
  geom_rect(aes(xmin=min(elast_data_vital$Elev),xmax=min_elev_pied,
                ymin=min(elast_data_vital$Elast_i),ymax=max(elast_data_vital$Elast_i)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(elast_data_vital$Elev),
                ymin=min(elast_data_vital$Elast_i),ymax=max(elast_data_vital$Elast_i)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  #scale_colour_manual(breaks=c("Growth","Survival","Recruitment"), 
  #                    values=c("Growth"="#1b9e77","Survival"="#d95f02","Recruitment"="#7570b3"),
  #                    labels=c("Growth","Survival","Recruitment"),
  #                    name="Vital rate")+
  scale_colour_manual(breaks=c("Survival","Recruitment"), 
                      values=c("Survival"="#1b9e77","Recruitment"="#7570b3"),
                      labels=c("Grow/Surv","Recr"),
                      name="Kernel")+
  labs(x = "Elevation (ft)", y = "Elasticity")+
  mytheme

## Save plots
ggsave(file="elast_vital_c.png", plot=elast_plot_vital_c,width=4,height=3,units="in",dpi=600)
ggsave(file="elast_plot_vital_cc.png", plot=elast_plot_vital_cc,width=4,height=3,units="in",dpi=600)
ggsave(file="elast_plot_vital_i.png", plot=elast_plot_vital_i,width=4,height=3,units="in",dpi=600)


### LTRE - vital rates ----

## Data
ltre_data_vital<-read.csv("./Output/ltre_data_vital.csv")
load("./Output/elev_limits.rda")

## Plots
ltre_plot_vital_c<-ggplot(data=ltre_data_vital,aes(x=Elev,y=Env_c,col=Rate))+
  geom_rect(aes(xmin=min(ltre_data_vital$Elev),xmax=min_elev_pied,
                ymin=min(ltre_data_vital$Env_c),ymax=max(ltre_data_vital$Env_c)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(ltre_data_vital$Elev),
                ymin=min(ltre_data_vital$Env_c),ymax=max(ltre_data_vital$Env_c)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Growth","Survival","Recruitment","Total"), 
                      values=c("Growth"="#1b9e77","Survival"="#d95f02","Recruitment"="#7570b3",
                               "Total"="Black"),
                      labels=c("Growth","Survival","Recruitment","Total"),
                      name="Vital rate")+
  labs(x = "Elevation (ft)", y = expression(paste("d ",lambda," / d Elevation")))+
  mytheme

ltre_plot_vital_cc<-ggplot(data=ltre_data_vital,aes(x=Elev,y=Env_cc,col=Rate))+
  geom_rect(aes(xmin=min(ltre_data_vital$Elev),xmax=min_elev_pied,
                ymin=min(ltre_data_vital$Env_cc),ymax=max(ltre_data_vital$Env_cc)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(ltre_data_vital$Elev),
                ymin=min(ltre_data_vital$Env_cc),ymax=max(ltre_data_vital$Env_cc)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Growth","Survival","Recruitment","Total"), 
                      values=c("Growth"="#1b9e77","Survival"="#d95f02","Recruitment"="#7570b3",
                               "Total"="Black"),
                      labels=c("Growth","Survival","Recruitment","Total"),
                      name="Vital rate")+
  labs(x = "Elevation (ft)", y = expression(paste("d ",lambda," / d Elevation")))+
  mytheme

ltre_plot_vital_i<-ggplot(data=ltre_data_vital,aes(x=Elev,y=Env_i,col=Rate))+
  geom_rect(aes(xmin=min(ltre_data_vital$Elev),xmax=min_elev_pied,
                ymin=min(ltre_data_vital$Env_i),ymax=max(ltre_data_vital$Env_i)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(ltre_data_vital$Elev),
                ymin=min(ltre_data_vital$Env_i),ymax=max(ltre_data_vital$Env_i)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("Growth","Survival","Recruitment","Total"), 
                      values=c("Growth"="#1b9e77","Survival"="#d95f02","Recruitment"="#7570b3",
                               "Total"="Black"),
                      labels=c("Growth","Survival","Recruitment","Total"),
                      name="Vital rate")+
  labs(x = "Elevation (ft)", y = expression(paste("d ",lambda," / d Elevation")))+
  mytheme

## Save plots
ggsave(file="ltre_plot_vital_c.png", plot=ltre_plot_vital_c,width=4,height=3,units="in",dpi=600)
ggsave(file="ltre_plot_vital_cc.png", plot=ltre_plot_vital_cc,width=4,height=3,units="in",dpi=600)
ggsave(file="ltre_plot_vital_i.png", plot=ltre_plot_vital_i,width=4,height=3,units="in",dpi=600)


### LTRE - predictors ----

## Data
ltre_data_env<-read.csv("./Output/ltre_data_env.csv")
load("./Output/elev_limits.rda")

## Plots
ltre_plot_c<-ggplot(data=ltre_data_env,aes(x=Elev,y=Env_c,col=Predictor))+
  geom_rect(aes(xmin=min(ltre_data_env$Elevation),xmax=min_elev_pied,
                ymin=min(ltre_data_env$Env_c),ymax=max(ltre_data_env$Env_c)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(ltre_data_env$Elevation),
                ymin=min(ltre_data_env$Env_c),ymax=max(ltre_data_env$Env_c)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("BA","MAP","MAT","Total"), 
                      values=c("BA"="#1b9e77","MAP"="#d95f02","MAT"="#7570b3",
                               "Total"="Black"),
                      labels=c("BA","MAP","MAT","Total"),
                      name="Predictor")+
  labs(x = "Elevation (ft)", y = expression(paste("d ",lambda," / d Elevation")))+
  mytheme

ltre_plot_cc<-ggplot(data=ltre_data_env,aes(x=Elev,y=Env_cc,col=Predictor))+
  geom_rect(aes(xmin=min(ltre_data_env$Elevation),xmax=min_elev_pied,
                ymin=min(ltre_data_env$Env_cc),ymax=max(ltre_data_env$Env_cc)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(ltre_data_env$Elevation),
                ymin=min(ltre_data_env$Env_cc),ymax=max(ltre_data_env$Env_cc)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("BA","MAP","MAT","Total"), 
                      values=c("BA"="#1b9e77","MAP"="#d95f02","MAT"="#7570b3",
                               "Total"="Black"),
                      labels=c("BA","MAP","MAT","Total"),
                      name="Predictor")+
  labs(x = "Elevation (ft)", y = expression(paste("d ",lambda," / d Elevation")))+
  mytheme

ltre_plot_i<-ggplot(data=ltre_data_env,aes(x=Elev,y=Env_i,col=Predictor))+
  geom_rect(aes(xmin=min(ltre_data_env$Elevation),xmax=min_elev_pied,
                ymin=min(ltre_data_env$Env_i),ymax=max(ltre_data_env$Env_i)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_rect(aes(xmin=max_elev_pied,xmax=max(ltre_data_env$Elevation),
                ymin=min(ltre_data_env$Env_i),ymax=max(ltre_data_env$Env_i)),
            fill="grey80",col="grey80",alpha=0.1)+
  geom_line(size=1)+
  scale_colour_manual(breaks=c("BA","MAP","MAT","Total"), 
                      values=c("BA"="#1b9e77","MAP"="#d95f02","MAT"="#7570b3",
                               "Total"="Black"),
                      labels=c("BA","MAP","MAT","Total"),
                      name="Predictor")+
  labs(x = "Elevation (ft)", y = expression(paste("d ",lambda," / d Elevation")))+
  mytheme

## Save plots
ggsave(file="ltre_plot_c.png", plot=ltre_plot_c,width=4,height=3,units="in",dpi=600)
ggsave(file="ltre_plot_cc.png", plot=ltre_plot_cc,width=4,height=3,units="in",dpi=600)
ggsave(file="ltre_plot_i.png", plot=ltre_plot_i,width=4,height=3,units="in",dpi=600)
