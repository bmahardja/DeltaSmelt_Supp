# Purpose: Analysis for understanding correlation between salinity and distance
# Author: Brian Mahardja
# Date: 2024-01-24

library(tidyverse)
library(lubridate)
library(viridis)
library(sharpshootR)
library(zoo)
library(rpart)
library(rpart.plot)
library(caret)
library(MultivariateRandomForest)

# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

fish_dist_data_swp <- read_csv(file.path(output_root,"wild_fish_SWP_distance_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count>=20)

fish_dist_data_cvp <- read_csv(file.path(output_root,"wild_fish_CVP_distance_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count>=20)

fish_dist_data_cache <- read_csv(file.path(output_root,"wild_fish_Cache_distance_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count>=20)

fish_sal_data <- read_csv(file.path(output_root,"wild_fish_salinity_data_with_covariate.csv")) %>%
  # Filter out all data with <20 fish
  filter(Count_Sal>=20)

# Combine data
fish_sal_combined <- fish_dist_data_swp %>% rename(Count_SWPdist = Count,dist_km_swp=dist_km_facility,dist_km_swp_sd=dist_km_sd) %>% dplyr::select(c(1:6)) %>% 
  dplyr::select(-Date_median) %>% left_join(fish_sal_data) %>% 
  left_join((fish_dist_data_cvp %>% select(WY,Month,dist_km_facility,Count) %>% rename(Count_CVPdist=Count,dist_km_cvp=dist_km_facility))) %>%
  left_join((fish_dist_data_cache %>% select(WY,Month,dist_km_facility,Count) %>% rename(Count_Cachedist=Count,dist_km_cache=dist_km_facility))) %>%
  filter(!is.na(Count_Sal))


ggplot(data=fish_sal_combined, aes(x=Sal_surf_mean,y=dist_km_swp)) + geom_point()

summary(lm(fish_sal_combined$Sal_surf_mean~fish_sal_combined$dist_km_swp))

ggplot(data=fish_sal_combined, aes(x=Sal_surf_mean,y=dist_km_cvp)) + geom_point()

summary(lm(fish_sal_combined$Sal_surf_mean~fish_sal_combined$dist_km_cvp))

ggplot(data=fish_sal_combined, aes(x=dist_km_swp,y=dist_km_cvp)) + geom_point()

summary(lm(fish_sal_combined$dist_km_swp~fish_sal_combined$dist_km_cvp))

ggplot(data=fish_sal_combined, aes(x=dist_km_cache,y=dist_km_swp)) + geom_point()

summary(lm(fish_sal_combined$dist_km_cache~fish_sal_combined$dist_km_swp))

#####################
# K-mean clustering
fish_kmean_data <- fish_sal_combined %>% 
  select(dist_km_swp,dist_km_cvp,dist_km_cache,Sal_surf_mean,X2_1,mean_OMR,OUT_1,Secchi_South,FirstFlushOccurred) %>%
  mutate(FirstFlushOccurred=ifelse(FirstFlushOccurred=="Yes",1,0))

km.out <- kmeans(fish_kmean_data, centers = 3, nstart = 20)
km.out

# Decide how many clusters to look at
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(fish_kmean_data, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

km.out4 <- kmeans(fish_kmean_data, centers = 4, nstart = 100)
km.out4
View(km.out4$centers)

km.out3 <- kmeans(fish_kmean_data, centers = 3, nstart = 100)
km.out3

# Add k group back
fish_sal_combined$K_cluster<-km.out4$cluster

bad_years<-fish_sal_combined %>% filter(K_cluster==1&Month %in% c(2,3,4)&FirstFlushOccurred=="No")

#####################
# K-mean clustering JUST FOR response variables
fish_kmean_data <- fish_sal_combined %>% 
  select(dist_km_swp,dist_km_cvp,dist_km_cache,Sal_surf_mean) %>%
  mutate(dist_km_swp_z = (dist_km_swp - mean(dist_km_swp))/sd(dist_km_swp),
                dist_km_cvp_z = (dist_km_cvp - mean(dist_km_cvp))/sd(dist_km_cvp),
                dist_km_cache_z = (dist_km_cache - mean(dist_km_cache))/sd(dist_km_cache),
                Sal_surf_mean_z = (Sal_surf_mean - mean(Sal_surf_mean))/sd(Sal_surf_mean)) %>%
  select(dist_km_swp_z,dist_km_cvp_z,dist_km_cache_z,Sal_surf_mean_z)


# Decide how many clusters to look at
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(fish_kmean_data, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

km.out5 <- kmeans(fish_kmean_data, centers = 5, nstart = 100)
km.out5
View(km.out5$centers)

km.out6 <- kmeans(fish_kmean_data, centers = 6, nstart = 100)
km.out6
View(km.out6$centers)

fish_sal_combined$Response_cluster_5<-km.out5$cluster
fish_sal_combined$Response_cluster_6<-km.out6$cluster

fish_sal_combined %>% group_by(Response_cluster_5) %>%
  summarise(dist_km_swp=mean(dist_km_swp),dist_km_cvp=mean(dist_km_cvp),
            dist_km_cache=mean(dist_km_cache),Sal_surf_mean=mean(Sal_surf_mean))

fish_sal_combined %>% group_by(Response_cluster_6) %>%
  summarise(dist_km_swp=mean(dist_km_swp),dist_km_cvp=mean(dist_km_cvp),
            dist_km_cache=mean(dist_km_cache),Sal_surf_mean=mean(Sal_surf_mean))

fish_sal_combined <- fish_sal_combined %>% mutate(
  Response_group = factor(case_when(Response_cluster_5 ==1 ~ "1 LowSal CloseCache",
                             Response_cluster_5 ==2 ~ "2 LowSal FarCache",
                             Response_cluster_5 ==3 ~ "3 VeryHighSal Central",
                             Response_cluster_5 ==4 ~ "4 0.9ppt Central",
                             Response_cluster_5 ==5 ~ "5 1.3ppt ClosetoFac NearCache")
))


trctrl_loo <- trainControl(
  method = 'LOOCV',                # k-fold cross validation 'cv'
  number = 1                    # number of folds
) 

tree_tc_loo <- train(Response_group ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month)+Secchi_South +Secchi_West+WY, 
                     data = fish_sal_combined, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo
#0.020


fit.tree.grp = rpart(Response_group ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West+WY, 
                       data=fish_sal_combined, 
                        cp=0.040)
rpart.plot(fit.tree.grp)
summary(fit.tree.grp)

tree_tc_loo <- train(as.factor(Response_cluster_6) ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month)+Secchi_South +Secchi_West+WY, 
                     data = fish_sal_combined, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo
#0.035


fit.tree.grp2 = rpart(as.factor(Response_cluster_6) ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West+WY, 
                     data=fish_sal_combined, 
                     cp=0.035)
rpart.plot(fit.tree.grp2)

fish_sal_combined %>% group_by(FirstFlushOccurred) %>%
  summarise(OUT_1=mean(OUT_1),OUT=mean(OUT),X2_1=mean(X2_1),X2=mean(X2))
