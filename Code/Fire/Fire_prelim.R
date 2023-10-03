### Fire
rdata <- read.csv("./Processed/Recruitment/RecruitData.csv", header = T, stringsAsFactors = F)
conds <- read.csv("./FIAdata/COND_COMBINED.csv", header = T, stringsAsFactors = F)

rdata$DSTRBCD1<-conds$DSTRBCD1[match(rdata$plot, conds$PLT_CN)]

survData <- read.csv("./Processed/Survival/SurvivalData.csv", header = T, stringsAsFactors = F)

survFire<-subset(survData,DSTRBCD1==30|DSTRBCD1==31|DSTRBCD1==32)

n_Fire<-length(unique(survFire$PLT_CN))
n_Plot<-length(unique(survData$PLT_CN))
prop_Fire<-n_Fire/n_Plot

Fire_summ_total<-mean(survFire$mort)

Fire_summ_plot<-survFire %>%
  group_by(PLT_CN) %>%
  summarise(mort=mean(mort))
  
Fire_summ_type<-survFire %>%
  group_by(DSTRBCD1) %>%
  summarise(mort=mean(mort))

survData$Fire<-ifelse(
  survData$DSTRBCD1 == 30|survData$DSTRBCD1 == 31|survData$DSTRBCD1 == 32, 1, 0)

rdata$Fire<-ifelse(
  rdata$DSTRBCD1 == 30|rdata$DSTRBCD1 == 31|rdata$DSTRBCD1 == 32, 1, 0)

Fire_summ_all<-rdata %>%
  group_by(PApied) %>%
  summarise(count=length(Fire),nFire=sum(Fire,na.rm=T),freq=(sum(Fire,na.rm=T)/length(Fire)))

plot(survData$PPT_yr_norm,survData$Fire)
plot(survData$T_yr_norm,survData$Fire)

plot(rdata$PPT_yr_norm,rdata$Fire)
plot(rdata$T_yr_norm,rdata$Fire)
plot(rdata$BALIVE,rdata$Fire)
plot(rdata$BA.PIPO,rdata$Fire)
plot(rdata$BA.notPIED,rdata$Fire)


