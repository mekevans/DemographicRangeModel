library(ggeffects)
library(ggplot2)
library(coefplot)

# Read and process data
rdata <- read.csv("D:/EvansLab/Final/Data/Processed/Recruitment/RecruitmentData.csv", stringsAsFactors = F, header = T)

# Subset data
rdata <- subset(rdata, PIED == 1 & AGB_intra > 0)

# Model: elevation + BA
rmodel <- glm(data = rdata, recruits ~
                offset(log(measInterval)) + offset(log(AGB_intra + 0.00001)), 
              family = "poisson")
summary(rmodel)
save(rmodel, file = "D:/EvansLab/Final/Models/B/recr.rda")
pdf("D:/EvansLab/Final/Output/B/RecruitmentModel.pdf")
coefplot(rmodel, intercept = F)
hist(rdata$recruits, col = rgb(1,0,0,0.75), main = "Observed (red) vs. modeled (green) recruitment count", xlab = "Number of recruits")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
png("D:/EvansLab/Final/Manuscript/FigS8.png")
hist(rdata$recruits, col = rgb(1,0,0,0.75), xlab = "Number of recruits", main = "")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()

# Model: elevation + climate
rmodel <- glm(data = rdata, recruits ~ elev + PPT_c_window_25 + PPT_w_window_25 + VPD_c_window_25 + VPD_w_window_25 +
                offset(log(measInterval)) + offset(log(AGB_intra + 0.00001)), 
              family = "poisson")
rmodel <- update(rmodel, . ~ . - VPD_c_window_25)
# rmodel <- update(rmodel, . ~ . + I(VPD_w_window_25))
rmodel <- update(rmodel, . ~ . - VPD_w_window_25)
rmodel <- update(rmodel, . ~ . - elev)
summary(rmodel)
save(rmodel, file = "D:/EvansLab/Final/Models/C/recr.rda")
pdf("D:/EvansLab/Final/Output/C/RecruitmentModel.pdf")
coefplot(rmodel, intercept = F)
hist(rdata$recruits, col = rgb(1,0,0,0.75), main = "Observed (red) vs. modeled (green) recruitment count", xlab = "Number of recruits")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
coefficientplot <- coefplot(rmodel, intercept = F) + scale_y_discrete(labels = c("Cool-season\nprecipitation", "Warm-season\nprecipitation")) +
  ggtitle("")
ggsave("D:/EvansLab/Final/Manuscript/FigS7.png", coefficientplot, scale = 1.5)
# Model: elevation + BA + climate
rmodel <- glm(data = rdata, recruits ~ elev + baLive + PPT_c_window_25 + PPT_w_window_25 + VPD_c_window_25 + VPD_w_window_25 +
                offset(log(measInterval)) + offset(log(AGB_intra + 0.00001)), 
              family = "poisson")
rmodel <- update(rmodel, . ~ . - VPD_w_window_25)
rmodel <- update(rmodel, . ~ . - VPD_c_window_25)
rmodel <- update(rmodel, . ~ . - elev)
rmodel <- update(rmodel, . ~ . - baLive)
summary(rmodel)
save(rmodel, file = "D:/EvansLab/Final/Models/BC/recr.rda")
pdf("D:/EvansLab/Final/Output/BC/RecruitmentModel.pdf")
coefplot(rmodel, intercept = F)
hist(rdata$recruits, col = rgb(1,0,0,0.75), main = "Observed (red) vs. modeled (green) recruitment count", xlab = "Number of recruits")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()
png("D:/EvansLab/Final/Manuscript/FigS9.png")
hist(rdata$recruits, col = rgb(1,0,0,0.75), xlab = "Number of recruits", main = "")
hist(rmodel$fitted.values, col = rgb(0,1,0,0.75), add = T)
dev.off()