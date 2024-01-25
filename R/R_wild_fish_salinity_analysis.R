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

tree_tc <- train(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West, data = fish_sal_data, method = 'rpart',
                 trControl = trctrl,      
                 tuneGrid = cp_grid)
tree_tc
#cp 0.045 being best

fit.tree = rpart(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West, data=fish_sal_data, 
                  method="anova", cp=0.045)
rpart.plot(fit.tree)
fit.tree$variable.importance


ggplot(data=fish_sal_data) + geom_boxplot(aes(x=FirstFlushOccurred,y=Sal_surf_mean))
ggplot(data=fish_sal_data) + geom_point(aes(x=X2_1,y=Sal_surf_mean,color=FirstFlushOccurred))


trctrl_loo <- trainControl(
  method = 'LOOCV',                # k-fold cross validation 'cv'
  number = 1                    # number of folds
) 

tree_tc_loo <- train(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month), data = fish_sal_data, method = 'rpart',
                 trControl = trctrl_loo,      
                 tuneGrid = cp_grid)
tree_tc_loo
#cp 0.045 seems to be the best

fit.tree = rpart(Sal_surf_mean ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) , data=fish_sal_data, 
                 method="anova", cp=0.045)
rpart.plot(fit.tree)


############### old script

dayflow_post97<-read.csv(file=file.path(data_root,"dayflow-results-1997-2020.csv")) %>%
  mutate(Date=as.Date(Date,"%m/%d/%Y"))

dayflow_pre97<-read.csv(file=file.path(data_root,"dayflowCalculations_pre97.csv"))  %>%
  mutate(Date=as.Date(Date,"%m/%d/%Y")) %>% rename(Month=Mo)

dayflow_merged<- bind_rows(dayflow_post97 %>% dplyr::select(Year,Month,Date,SAC,EXPORTS,OUT,TOT,X2),
                           dayflow_pre97 %>% dplyr::select(Year,Month,Date,SAC,EXPORTS,OUT,TOT,X2))

monthly_avg<- dayflow_merged %>% mutate(WY=ifelse(Month>9,year(Date)+1,year(Date))) %>% filter(Month %in% c(11,12,1:3)) %>%
  group_by(WY) %>%
  summarise(mean_OUT=mean(OUT))

first_flush <- dayflow_merged  %>% mutate(WY=ifelse(Month>9,year(Date)+1,year(Date))) %>% filter(Month %in% c(11,12,1:3)) %>% left_join(monthly_avg) %>% mutate(OUT_anomaly=OUT-mean_OUT)

export_data <- dayflow_merged %>% mutate(WY=ifelse(Month>9,year(Date)+1,year(Date))) %>% group_by(WY,Month) %>%
  summarise(mean_export=mean(EXPORTS))

x2_data <- dayflow_merged  %>% mutate(WY=ifelse(Month>9,year(Date)+1,year(Date)),X2_77=ifelse(X2<77,1,0)) %>% filter(X2_77>0) %>%
  group_by(WY) %>% summarise(FirstFlush_X2=min(Date))

#Freeport flow data
FPT_daily <- CDECquery(id='FPT', sensor=20, interval='D', start=as.Date("1976-10-01"), end=as.Date("2017-12-31"))
ggplot(data=FPT_daily) + geom_point(aes(x=datetime,y=value))
str(FPT_daily)

FPT_data<- FPT_daily %>% mutate(Date=as.Date(datetime),Freeport_flow_CDEC=as.numeric(value)) %>%
  left_join((dayflow_merged %>%  dplyr::select(Date,SAC))) %>%
  mutate(Freeport_difference=SAC-Freeport_flow_CDEC) %>%
  mutate(Freeport_flow_final = ifelse(Freeport_flow_CDEC==0|Freeport_difference>=10000|Freeport_difference<=-10000,SAC,Freeport_flow_CDEC))

FPT_data_final <- FPT_data %>% dplyr::select(water_year,Date,Freeport_flow_final) %>%
  mutate(day_average_Freeport_flow = zoo::rollmean(Freeport_flow_final,3,align = "right", na.pad = TRUE)) %>%
  mutate(above_25000=ifelse(day_average_Freeport_flow>25000,1,0))

str(FPT_data_final)

FPT_data_firstflush <- FPT_data_final %>% group_by(water_year) %>% filter(month(Date) %in% c(11,12,1,2,3,4)) %>%
  filter(above_25000>0) %>% summarise(FirstFlush_firstDate=min(Date),FirstFlush_medianDate=median(Date)) %>%
  rename(WY=water_year)

fish_sal_data<-read.csv(file.path(output_root,"wild_fish_salinity_data.csv")) %>% mutate(Date_median=as.Date(Date_median)) %>%
  left_join(FPT_data_firstflush) %>%
  left_join(export_data) %>%
  left_join(dayflow_merged %>% dplyr::select(Date,X2)  %>% rename(Date_median=Date)) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(x2_data) %>%
  mutate(FirstFlushX2_Occurred=as.factor(ifelse(Date_median>=FirstFlush_X2,"Yes","No")))

str(fish_sal_data)

fish_sal_data <- fish_sal_data %>% filter(!is.na(FirstFlushOccurred))

hist(fish_sal_data$Sal_surf_mean,breaks=40)
fish_sal_data <- fish_sal_data %>% mutate(Sal_category=case_when(
  Sal_surf_mean<=1 ~ "<= 1",
  Sal_surf_mean<=4&Sal_surf_mean>1 ~ "1 to 4",
  Sal_surf_mean>4 ~ ">4",
))

str(fish_sal_data)
fish_sal_data$FirstFlushOccurred<- as.factor(fish_sal_data$FirstFlushOccurred)

fit.tree = rpart(Sal_category ~ X2 + FirstFlushOccurred + mean_export + as.factor(Month), data=fish_sal_data, method = "class", cp=0.008)
rpart.plot(fit.tree)
fit.tree
fit.tree$variable.importance

fit.tree2 = rpart(Sal_surf_mean ~ X2 + FirstFlushOccurred + mean_export + as.factor(Month) + WY, data=fish_sal_data, 
                  method="anova", cp=0.01, xval=7)


fit.tree2 = rpart(Sal_surf_mean ~ X2 + FirstFlushOccurred + mean_export + as.factor(Month) + WY +FirstFlushX2_Occurred, data=fish_sal_data, 
                  method="anova", cp=0.05, xval=0)
rpart.plot(fit.tree2)
fit.tree2$variable.importance
summary(fit.tree2)

trctrl <- trainControl(method = "cv", #cross validation
                       number = 10)   #10-fold cross validation
cp_grid <- data.frame(cp = seq(0.005, .2, .005))

tree_03 <- train(Sal_surf_mean~X2 + FirstFlushOccurred + mean_export + as.factor(Month) + WY, data = fish_sal_data, method = 'rpart',
                 trControl = trctrl,      
                 tuneGrid = cp_grid)

tree_03


fit.tree2 = rpart(Sal_surf_mean ~ X2 + FirstFlushOccurred + mean_export + as.factor(Month) + WY, data=fish_sal_data, 
                  method="anova")
rpart.plot(fit.tree2)
summary(fit.tree2)
zp <- prune(fit.tree2, cp = 0.06)
rpart.plot(zp)


fit.tree3 = rpart(Sal_surf_sd ~ FirstFlushOccurred +FirstFlushX2_Occurred + WY, data=fish_sal_data, 
                  method="anova", cp=0.05, xval=0)
rpart.plot(fit.tree3)
summary(fit.tree3)




test_mod <- lm(data=fish_sal_data,Sal_surf_mean ~ X2 * FirstFlushOccurred)
summary(test_mod)
plot(test_mod)
test_mod2 <- lm(data=fish_sal_data,Sal_surf_mean ~ X2 * FirstFlushOccurred + mean_export)
summary(test_mod2)
test_mod3 <- lm(data=fish_sal_data,Sal_surf_mean ~FirstFlushOccurred)
summary(test_mod3)
## Look into DOP tech report and years when migratory fish are more common

fish_sal_data$mod1predict<-predict(test_mod,newdata=fish_sal_data)
ggplot(data=fish_sal_data) + geom_point(aes(x=Sal_surf_mean,y=mod1predict))

fish_sal_data$residual<-fish_sal_data$Sal_surf_mean-fish_sal_data$mod1predict
ggplot(data=fish_sal_data) + geom_point(aes(x=Date_median,y=residual))
ggplot(data=fish_sal_data) + geom_point(aes(x=mod1predict,y=residual))

ggplot(data=fish_sal_data) + geom_point(aes(x=X2,y=Sal_surf_mean))
