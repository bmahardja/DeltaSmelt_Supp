# Purpose: Analysis for predicting salinity that wild fish occupy
# Author: Brian Mahardja
# Date: 2023-12-13

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

fish_dist_data <- read_csv(file.path(output_root,"wild_fish_SWP_distance_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count>=20)

trctrl <- trainControl(method = "cv", #cross validation
                       number = 10)   #10-fold cross validation
cp_grid <- data.frame(cp = seq(0.005, .2, .005))

tree_tc <- train(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month), data = fish_dist_data, method = 'rpart',
                 trControl = trctrl,      
                 tuneGrid = cp_grid)
tree_tc
#0.035


trctrl_loo <- trainControl(
  method = 'LOOCV',                # k-fold cross validation 'cv'
  number = 1                    # number of folds
) 

tree_tc_loo <- train(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month), data = fish_dist_data, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo
#cp 0.030 or 0.040 seems to be the best


fit.tree = rpart(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) , data=fish_dist_data, 
                 method="anova", cp=0.030)
rpart.plot(fit.tree)


fit.tree2 = rpart(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) , data=fish_dist_data, 
                 method="anova", cp=0.035)
rpart.plot(fit.tree2)
summary(fit.tree2)



###########


fish_dist_data_cvp <- read_csv(file.path(output_root,"wild_fish_CVP_distance_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count>=20)

trctrl <- trainControl(method = "cv", #cross validation
                       number = 10)   #10-fold cross validation
cp_grid <- data.frame(cp = seq(0.005, .2, .005))

tree_tc <- train(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month), 
                 data = fish_dist_data_cvp, method = 'rpart',
                 trControl = trctrl,      
                 tuneGrid = cp_grid)
tree_tc
#0.020


trctrl_loo <- trainControl(
  method = 'LOOCV',                # k-fold cross validation 'cv'
  number = 1                    # number of folds
) 

tree_tc_loo <- train(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month), 
                     data = fish_dist_data_cvp, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo
#cp 0.045 to 0.070 seems to be the best


fit.tree.cvp = rpart(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) , 
                     data=fish_dist_data_cvp, 
                  method="anova", cp=0.020)
rpart.plot(fit.tree.cvp)

fit.tree.cvp2 = rpart(dist_km_facility ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) , 
                     data=fish_dist_data_cvp, 
                     method="anova", cp=0.045)
rpart.plot(fit.tree.cvp2)
