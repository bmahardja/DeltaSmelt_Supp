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

# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

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

set.seed(49)

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

fish_sal_combined$Response_cluster_5<-km.out5$cluster

fish_sal_combined %>% group_by(Response_cluster_5) %>%
  summarise(dist_km_swp=mean(dist_km_swp),dist_km_cvp=mean(dist_km_cvp),
            dist_km_cache=mean(dist_km_cache),Sal_surf_mean=mean(Sal_surf_mean))

fish_sal_combined <- fish_sal_combined %>% mutate(
  Response_group = factor(case_when(Response_cluster_5 ==1 ~ "1 0.5ppt CloseCache",
                             Response_cluster_5 ==2 ~ "2 3.8ppt Central",
                             Response_cluster_5 ==3 ~ "3 0.5ppt FarCache",
                             Response_cluster_5 ==4 ~ "4 0.9ppt Central",
                             Response_cluster_5 ==5 ~ "5 1.3ppt ClosetoFac NearCache")
))


trctrl_loo <- trainControl(
  method = 'LOOCV',                # k-fold cross validation 'cv'
  number = 1                    # number of folds
) 

cp_grid <- data.frame(cp = seq(0.005, .2, .005))

tree_tc_loo <- train(Response_group ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month)+Secchi_South +Secchi_West+WY, 
                     data = fish_sal_combined, method = 'rpart',
                     trControl = trctrl_loo,      
                     tuneGrid = cp_grid)
tree_tc_loo
#About 0.040


fit.tree.grp = rpart(Response_group ~ SAC + SJR + TOT + OUT + EXPORTS + X2 + SAC_1 + SJR_1 + TOT_1 + OUT_1 + EXPORTS_1 + X2_1 + mean_OMR + FirstFlushOccurred + as.factor(Month) +Secchi_South +Secchi_West+WY, 
                       data=fish_sal_combined, 
                        cp=0.040)
rpart.plot(fit.tree.grp)
summary(fit.tree.grp)