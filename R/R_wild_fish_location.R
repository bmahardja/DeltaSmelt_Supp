# Purpose: Survey data to assess swimming movement in SKT
# Author: Brian Mahardja
# Date: 2023-09-22

# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_SwimDistance"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

# Load Package Libraries
library(gdistance)
library(raster)
library(sf)
library(sp)
library(rgdal)
library(tidyverse)
library(readxl)
library(ggspatial)
library(lubridate)
library(viridis)
library(readr)
library(rgeos)
library(deltafish)
library(mapview)


# Read the pre-processed raster file from Ryan since the conversion takes 2 hours 
raster_baydelta_tr <- readRDS(file.path(data_root,"raster_baydelta_tr.rds"))

# open our two data files
surv <- open_survey()
fish <- open_fish()

# filter for sources and taxa of interest
fish_smelt <- fish %>% 
  filter(Taxa %in% c("Hypomesus transpacificus"))

# do a join and collect the resulting data frame
# collect executes the sql query and gives you a table
df <- left_join(surv, fish_smelt) %>% 
  collect() 

# Subset to just December to April
# open our two data files
surv <- open_survey()
fish <- open_fish()

# filter for sources and taxa of interest
fish_smelt <- fish %>% 
    filter(Taxa %in% c("Hypomesus transpacificus"))


# do a join and collect the resulting data frame
# collect executes the sql query and gives you a table
df <- left_join(surv, fish_smelt) %>% 
    collect() %>% filter(Count>0)

# Subset to just December-April
df_subset <- df %>% mutate(WY=ifelse(month(Date)>9,year(Date)+1,year(Date))) %>% filter(month(Date) %in% c(12,1:4)) %>%
  # Filter out NA latitude and longitude
  filter(!is.na(Longitude)&!is.na(Latitude)&!is.na(Sal_surf))

hist(df_subset$Length)
# For take, smelt is considered "adults" at 58 mm
# Histogram shows most fish in these months were from 50mm and above, so we will use that
df_subset_v2 <- df_subset %>% filter(Length>=50) %>% mutate(Month=month(Date),Latitude_mult=Latitude*Count,Longitude_mult=Longitude*Count) %>% group_by(WY, Month) %>%
  summarise(Count=sum(Count),Date_median=median(Date),Latitude_sum=sum(Latitude_mult),Longitude_sum=sum(Longitude_mult),Sal_surf_sum=sum(Sal_surf),Sal_max=max(Sal_surf,na.rm=T),Sal_min=min(Sal_surf,na.rm=T)) %>%
  mutate(Latitude_c = Latitude_sum/Count, Longitude_c= Longitude_sum/Count,Sal_surf_mean=Sal_surf_sum/Count)


df_subset_v3 <- df_subset %>% filter(Length>=50)

df_subset_v3<-df_subset_v3[rep(row.names(df_subset_v3),df_subset_v3$Count),1:28]
  
df_subset_v3 <- df_subset_v3 %>% mutate(Month=month(Date),Count=1) %>% group_by(WY, Month) %>%
  summarise(Count=sum(Count),Date_median=median(Date),Sal_surf_mean=mean(Sal_surf),Sal_max=max(Sal_surf,na.rm=T),Sal_min=min(Sal_surf,na.rm=T),Sal_surf_sd=sd(Sal_surf)) %>%
  mutate(Sal_surf_cv=Sal_surf_sd/Sal_surf_mean)

ggplot(data=df_subset_v3) + geom_point(aes(x=Month,y=Sal_surf_mean))
ggplot(data=df_subset_v3) + geom_point(aes(x=Month,y=Sal_surf_sd))

df_subset_v3_2 <- df_subset_v3 %>% group_by(Month) %>% summarise(Sal_surf_mean=mean(Sal_surf_mean),Sal_surf_sd=mean(Sal_surf_sd,na.rm=T))

##NMDS
df_skt <- left_join(surv, fish_smelt) %>% 
  collect() %>% filter(Source=="SKT") %>% filter(!is.na(Latitude)&year(Date)<2015)

skt_mean_vol <- df_skt %>% filter(Tow_volume>0)  
mean(skt_mean_vol$Tow_volume)

df_skt_sum <- df_skt %>% mutate(Year=year(Date),Month=month(Date)) %>% group_by(Station,Year,Month) %>%
  mutate(Tow_volume=ifelse(Tow_volume<0,mean(skt_mean_vol$Tow_volume),Tow_volume)) %>% summarise(Count=sum(Count),Tow_volume=mean(Tow_volume)) %>%
  mutate(CPUE=Count/Tow_volume*1000)

df_skt_vegan_v1 <- df_skt_sum %>% select(-Tow_volume,-Count) %>% spread(Station,CPUE) %>%
  filter(Month %in% c(1:4)) 

df_skt_vegan_v1_1 <- df_skt_vegan_v1[complete.cases(df_skt_vegan_v1),]
df_skt_vegan_v1_2 <- df_skt_vegan_v1[ , colSums(is.na(df_skt_vegan_v1))==0]


df_skt_vegan_v2 <- df_skt_vegan_v1_2 %>% ungroup() 
df_skt_vegan_v2 <- df_skt_vegan_v2 %>% select_if(colSums(.) != 0) 

df_skt_vegan_v2_2 <- df_skt_vegan_v2%>% select(-Month,-Year)
#df_skt_vegan_v2$sum <- rowSums(df_skt_vegan_v2[])



library(vegan)

skt.nmds=metaMDS(df_skt_vegan_v2_2,dist="bray",k=2, autotransform=F)
skt.nmds

plot(skt.nmds)
orditorp(skt.nmds,display="species",col="red",air=0.01)

score_data<-scores(skt.nmds)
score_data<-as.data.frame(score_data[["sites"]])
score_data$Month<-df_skt_vegan_v2$Month
score_data$Year<-df_skt_vegan_v2$Year
hull_cyl <- score_data %>%
  group_by(Month) %>%
  slice(chull(NMDS1,NMDS2))

ggplot(score_data) + geom_point(aes(x=NMDS1,y=NMDS2,color=factor(Month))) +
  geom_polygon(data=hull_cyl,aes(x=NMDS1,y=NMDS2,fill=factor(Month)),alpha = 0.5)

                                