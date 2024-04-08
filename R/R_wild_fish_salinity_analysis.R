# Purpose: Analysis for predicting salinity that wild fish occupy
# Author: Brian Mahardja
# Date: 2023-12-01

library(tidyverse)
library(lubridate)
library(viridis)
library(sharpshootR)
library(zoo)
library(rpart)
library(rpart.plot)
library(caret)

# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

fish_sal_data <- read_csv(file.path(output_root,"wild_fish_salinity_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count_Sal>=20)
str(fish_sal_data)


trctrl <- trainControl(method = "cv", #cross validation
                       number = 10)   #10-fold cross validation
cp_grid <- data.frame(cp = seq(0.005, .2, .005))

tree_tc <- train(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West +WY, data = fish_sal_data, method = 'rpart',
                 trControl = trctrl,      
                 tuneGrid = cp_grid)
tree_tc
#cp 0.050 being best


ggplot(data=fish_sal_data) + geom_boxplot(aes(x=FirstFlushOccurred,y=Sal_surf_mean))
ggplot(data=fish_sal_data) + geom_point(aes(x=X2_1,y=Sal_surf_mean,color=FirstFlushOccurred))


trctrl_loo <- trainControl(
  method = 'LOOCV',                # k-fold cross validation 'cv'
  number = 1                    # number of folds
) 

tree_tc_loo <- train(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month)+Secchi_South +Secchi_West +WY, data = fish_sal_data, method = 'rpart',
                 trControl = trctrl_loo,      
                 tuneGrid = cp_grid)
tree_tc_loo
#cp 0.050 seems to be the best

fit.tree = rpart(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West +WY, data=fish_sal_data, 
                 method="anova", cp=0.050)
rpart.plot(fit.tree)

#######

tree_tc_loo <- train(FirstFlushOccurred ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + Sal_surf_mean + as.factor(Month)+Secchi_South +Secchi_West +WY, data = fish_sal_data, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo
#cp 0.050 seems to be best

fit.tree.iewpp = rpart(FirstFlushOccurred ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + as.factor(Month)+Secchi_South +Secchi_West +WY, data=fish_sal_data, 
           cp=0.025)
rpart.plot(fit.tree.iewpp)

hist(fish_sal_data$Sal_surf_mean,n=100)


### Divide data into 2
fish_sal_data$SalGrp<-as.factor(ifelse(fish_sal_data$Sal_surf_mean<=1,"1 ppt or less",">1 ppt"))

tree_tc_loo <- train(SalGrp ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + as.factor(Month)+Secchi_South +Secchi_West +WY, data = fish_sal_data, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo

fit.tree.iewpp = rpart(SalGrp ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + as.factor(Month)+Secchi_South +Secchi_West +WY, data=fish_sal_data, 
                       cp=0.050)

rpart.plot(fit.tree.iewpp)
