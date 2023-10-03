library(pROC)
library(tidyverse)
library(RColorBrewer)

##AUC
#Total
AUC_c <- roc(response = FIA_lambda$PApied, predictor = FIA_lambda$lambda_c)
AUC_ci <- roc(response = FIA_lambda$PApied, predictor = FIA_lambda$lambda_ci)
AUC_cc <- roc(response = FIA_lambda$PApied, predictor = FIA_lambda$lambda_cc)
AUC_i <- roc(response = FIA_lambda$PApied, predictor = FIA_lambda$lambda_i)

#Quartiles
FIA_lambda$q_elev<-cut(FIA_lambda$elev,quantile(FIA_lambda$elev),labels=c("1","2","3","4"))
FIA_lambda$q_ba<-cut(FIA_lambda$BALIVE,quantile(FIA_lambda$BALIVE),labels=c("1","2","3","4"))
FIA_lambda$q_ppt<-cut(FIA_lambda$PPT_yr_norm,quantile(FIA_lambda$PPT_yr_norm),labels=c("1","2","3","4"))
FIA_lambda$q_t<-cut(FIA_lambda$T_yr_norm,quantile(FIA_lambda$T_yr_norm),labels=c("1","2","3","4"))

AUC_c_elev1 <- roc(response = subset(FIA_lambda,q_elev==1)$PApied, predictor = subset(FIA_lambda,q_elev==1)$lambda_c)
AUC_c_elev2 <- roc(response = subset(FIA_lambda,q_elev==2)$PApied, predictor = subset(FIA_lambda,q_elev==2)$lambda_c)
AUC_c_elev3 <- roc(response = subset(FIA_lambda,q_elev==3)$PApied, predictor = subset(FIA_lambda,q_elev==3)$lambda_c)
AUC_c_elev4 <- roc(response = subset(FIA_lambda,q_elev==4)$PApied, predictor = subset(FIA_lambda,q_elev==4)$lambda_c)

AUC_ci_elev1 <- roc(response = subset(FIA_lambda,q_elev==1)$PApied, predictor = subset(FIA_lambda,q_elev==1)$lambda_ci)
AUC_ci_elev2 <- roc(response = subset(FIA_lambda,q_elev==2)$PApied, predictor = subset(FIA_lambda,q_elev==2)$lambda_ci)
AUC_ci_elev3 <- roc(response = subset(FIA_lambda,q_elev==3)$PApied, predictor = subset(FIA_lambda,q_elev==3)$lambda_ci)
AUC_ci_elev4 <- roc(response = subset(FIA_lambda,q_elev==4)$PApied, predictor = subset(FIA_lambda,q_elev==4)$lambda_ci)

AUC_cc_elev1 <- roc(response = subset(FIA_lambda,q_elev==1)$PApied, predictor = subset(FIA_lambda,q_elev==1)$lambda_cc)
AUC_cc_elev2 <- roc(response = subset(FIA_lambda,q_elev==2)$PApied, predictor = subset(FIA_lambda,q_elev==2)$lambda_cc)
AUC_cc_elev3 <- roc(response = subset(FIA_lambda,q_elev==3)$PApied, predictor = subset(FIA_lambda,q_elev==3)$lambda_cc)
AUC_cc_elev4 <- roc(response = subset(FIA_lambda,q_elev==4)$PApied, predictor = subset(FIA_lambda,q_elev==4)$lambda_cc)

AUC_i_elev1 <- roc(response = subset(FIA_lambda,q_elev==1)$PApied, predictor = subset(FIA_lambda,q_elev==1)$lambda_i)
AUC_i_elev2 <- roc(response = subset(FIA_lambda,q_elev==2)$PApied, predictor = subset(FIA_lambda,q_elev==2)$lambda_i)
AUC_i_elev3 <- roc(response = subset(FIA_lambda,q_elev==3)$PApied, predictor = subset(FIA_lambda,q_elev==3)$lambda_i)
AUC_i_elev4 <- roc(response = subset(FIA_lambda,q_elev==4)$PApied, predictor = subset(FIA_lambda,q_elev==4)$lambda_i)

AUC_c_ba1 <- roc(response = subset(FIA_lambda,q_ba==1)$PApied, predictor = subset(FIA_lambda,q_ba==1)$lambda_c)
AUC_c_ba2 <- roc(response = subset(FIA_lambda,q_ba==2)$PApied, predictor = subset(FIA_lambda,q_ba==2)$lambda_c)
AUC_c_ba3 <- roc(response = subset(FIA_lambda,q_ba==3)$PApied, predictor = subset(FIA_lambda,q_ba==3)$lambda_c)
AUC_c_ba4 <- roc(response = subset(FIA_lambda,q_ba==4)$PApied, predictor = subset(FIA_lambda,q_ba==4)$lambda_c)

AUC_ci_ba1 <- roc(response = subset(FIA_lambda,q_ba==1)$PApied, predictor = subset(FIA_lambda,q_ba==1)$lambda_ci)
AUC_ci_ba2 <- roc(response = subset(FIA_lambda,q_ba==2)$PApied, predictor = subset(FIA_lambda,q_ba==2)$lambda_ci)
AUC_ci_ba3 <- roc(response = subset(FIA_lambda,q_ba==3)$PApied, predictor = subset(FIA_lambda,q_ba==3)$lambda_ci)
AUC_ci_ba4 <- roc(response = subset(FIA_lambda,q_ba==4)$PApied, predictor = subset(FIA_lambda,q_ba==4)$lambda_ci)

AUC_cc_ba1 <- roc(response = subset(FIA_lambda,q_ba==1)$PApied, predictor = subset(FIA_lambda,q_ba==1)$lambda_cc)
AUC_cc_ba2 <- roc(response = subset(FIA_lambda,q_ba==2)$PApied, predictor = subset(FIA_lambda,q_ba==2)$lambda_cc)
AUC_cc_ba3 <- roc(response = subset(FIA_lambda,q_ba==3)$PApied, predictor = subset(FIA_lambda,q_ba==3)$lambda_cc)
AUC_cc_ba4 <- roc(response = subset(FIA_lambda,q_ba==4)$PApied, predictor = subset(FIA_lambda,q_ba==4)$lambda_cc)

AUC_i_ba1 <- roc(response = subset(FIA_lambda,q_ba==1)$PApied, predictor = subset(FIA_lambda,q_ba==1)$lambda_i)
AUC_i_ba2 <- roc(response = subset(FIA_lambda,q_ba==2)$PApied, predictor = subset(FIA_lambda,q_ba==2)$lambda_i)
AUC_i_ba3 <- roc(response = subset(FIA_lambda,q_ba==3)$PApied, predictor = subset(FIA_lambda,q_ba==3)$lambda_i)
AUC_i_ba4 <- roc(response = subset(FIA_lambda,q_ba==4)$PApied, predictor = subset(FIA_lambda,q_ba==4)$lambda_i)

AUC_c_ppt1 <- roc(response = subset(FIA_lambda,q_ppt==1)$PApied, predictor = subset(FIA_lambda,q_ppt==1)$lambda_c)
AUC_c_ppt2 <- roc(response = subset(FIA_lambda,q_ppt==2)$PApied, predictor = subset(FIA_lambda,q_ppt==2)$lambda_c)
AUC_c_ppt3 <- roc(response = subset(FIA_lambda,q_ppt==3)$PApied, predictor = subset(FIA_lambda,q_ppt==3)$lambda_c)
AUC_c_ppt4 <- roc(response = subset(FIA_lambda,q_ppt==4)$PApied, predictor = subset(FIA_lambda,q_ppt==4)$lambda_c)

AUC_ci_ppt1 <- roc(response = subset(FIA_lambda,q_ppt==1)$PApied, predictor = subset(FIA_lambda,q_ppt==1)$lambda_ci)
AUC_ci_ppt2 <- roc(response = subset(FIA_lambda,q_ppt==2)$PApied, predictor = subset(FIA_lambda,q_ppt==2)$lambda_ci)
AUC_ci_ppt3 <- roc(response = subset(FIA_lambda,q_ppt==3)$PApied, predictor = subset(FIA_lambda,q_ppt==3)$lambda_ci)
AUC_ci_ppt4 <- roc(response = subset(FIA_lambda,q_ppt==4)$PApied, predictor = subset(FIA_lambda,q_ppt==4)$lambda_ci)

AUC_cc_ppt1 <- roc(response = subset(FIA_lambda,q_ppt==1)$PApied, predictor = subset(FIA_lambda,q_ppt==1)$lambda_cc)
AUC_cc_ppt2 <- roc(response = subset(FIA_lambda,q_ppt==2)$PApied, predictor = subset(FIA_lambda,q_ppt==2)$lambda_cc)
AUC_cc_ppt3 <- roc(response = subset(FIA_lambda,q_ppt==3)$PApied, predictor = subset(FIA_lambda,q_ppt==3)$lambda_cc)
AUC_cc_ppt4 <- roc(response = subset(FIA_lambda,q_ppt==4)$PApied, predictor = subset(FIA_lambda,q_ppt==4)$lambda_cc)

AUC_i_ppt1 <- roc(response = subset(FIA_lambda,q_ppt==1)$PApied, predictor = subset(FIA_lambda,q_ppt==1)$lambda_i)
AUC_i_ppt2 <- roc(response = subset(FIA_lambda,q_ppt==2)$PApied, predictor = subset(FIA_lambda,q_ppt==2)$lambda_i)
AUC_i_ppt3 <- roc(response = subset(FIA_lambda,q_ppt==3)$PApied, predictor = subset(FIA_lambda,q_ppt==3)$lambda_i)
AUC_i_ppt4 <- roc(response = subset(FIA_lambda,q_ppt==4)$PApied, predictor = subset(FIA_lambda,q_ppt==4)$lambda_i)

AUC_c_t1 <- roc(response = subset(FIA_lambda,q_t==1)$PApied, predictor = subset(FIA_lambda,q_t==1)$lambda_c)
AUC_c_t2 <- roc(response = subset(FIA_lambda,q_t==2)$PApied, predictor = subset(FIA_lambda,q_t==2)$lambda_c)
AUC_c_t3 <- roc(response = subset(FIA_lambda,q_t==3)$PApied, predictor = subset(FIA_lambda,q_t==3)$lambda_c)
AUC_c_t4 <- roc(response = subset(FIA_lambda,q_t==4)$PApied, predictor = subset(FIA_lambda,q_t==4)$lambda_c)

AUC_ci_t1 <- roc(response = subset(FIA_lambda,q_t==1)$PApied, predictor = subset(FIA_lambda,q_t==1)$lambda_ci)
AUC_ci_t2 <- roc(response = subset(FIA_lambda,q_t==2)$PApied, predictor = subset(FIA_lambda,q_t==2)$lambda_ci)
AUC_ci_t3 <- roc(response = subset(FIA_lambda,q_t==3)$PApied, predictor = subset(FIA_lambda,q_t==3)$lambda_ci)
AUC_ci_t4 <- roc(response = subset(FIA_lambda,q_t==4)$PApied, predictor = subset(FIA_lambda,q_t==4)$lambda_ci)

AUC_cc_t1 <- roc(response = subset(FIA_lambda,q_t==1)$PApied, predictor = subset(FIA_lambda,q_t==1)$lambda_cc)
AUC_cc_t2 <- roc(response = subset(FIA_lambda,q_t==2)$PApied, predictor = subset(FIA_lambda,q_t==2)$lambda_cc)
AUC_cc_t3 <- roc(response = subset(FIA_lambda,q_t==3)$PApied, predictor = subset(FIA_lambda,q_t==3)$lambda_cc)
AUC_cc_t4 <- roc(response = subset(FIA_lambda,q_t==4)$PApied, predictor = subset(FIA_lambda,q_t==4)$lambda_cc)

AUC_i_t1 <- roc(response = subset(FIA_lambda,q_t==1)$PApied, predictor = subset(FIA_lambda,q_t==1)$lambda_i)
AUC_i_t2 <- roc(response = subset(FIA_lambda,q_t==2)$PApied, predictor = subset(FIA_lambda,q_t==2)$lambda_i)
AUC_i_t3 <- roc(response = subset(FIA_lambda,q_t==3)$PApied, predictor = subset(FIA_lambda,q_t==3)$lambda_i)
AUC_i_t4 <- roc(response = subset(FIA_lambda,q_t==4)$PApied, predictor = subset(FIA_lambda,q_t==4)$lambda_i)

AUC_c_plot<-data.frame(Model="c",Threshold=AUC_c$thresholds,F_pos=(1-AUC_c$specificities),
                       T_pos=AUC_c$sensitivities)
AUC_cc_plot<-data.frame(Model="cc",Threshold=AUC_cc$thresholds,F_pos=(1-AUC_cc$specificities),
                       T_pos=AUC_cc$sensitivities)
AUC_i_plot<-data.frame(Model="i",Threshold=AUC_i$thresholds,F_pos=(1-AUC_i$specificities),
                       T_pos=AUC_i$sensitivities)
AUC_plot_data<-rbind(AUC_c_plot,AUC_cc_plot,AUC_i_plot)

AUC_plot<-ggplot(data=AUC_plot_data,aes(x=F_pos,y=T_pos,col=Model))+
  geom_line(size=1)+
  geom_abline(slope=1,intercept=0)+
  scale_colour_manual(breaks=c("c","cc","i"), 
                      values=c("c"="#1b9e77","cc"="#d95f02","i"="#7570b3"),
                      labels=c("Clim","ClimInt","ClimCompInt"),
                      name="Model")+
  labs(x="False positive rate",y="True positive rate")+
  annotate("label", x = c(0.23,0.56,0.77), y = c(0.9,0.87,0.75), 
           label = c(paste("AUC =",round(AUC_c$auc,2)),
                     paste("AUC =",round(AUC_cc$auc,2)),
                     paste("AUC =",round(AUC_i$auc,2))),
           colour=c("#1b9e77","#d95f02","#7570b3"),size=3)

save_plot(file="AUC.png",AUC_plot,base_asp=1.5)

plot(AUC_c_plot$Threshold,AUC_c_plot$F_pos)
plot(AUC_c_plot$Threshold,AUC_c_plot$T_pos)

plot(AUC_cc_plot$Threshold,AUC_cc_plot$F_pos)
plot(AUC_cc_plot$Threshold,AUC_cc_plot$T_pos)

plot(AUC_i_plot$Threshold,AUC_i_plot$F_pos)
plot(AUC_i_plot$Threshold,AUC_i_plot$T_pos)
