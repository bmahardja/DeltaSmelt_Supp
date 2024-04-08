# Purpose: Prepare swim distance data for analysis for wild fish/historical data
# Author: Brian Mahardja
# Date: 2023-11-22

# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_Supp"
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



# Read the TIF file just to get the crs
raster_baydelta<-raster(file.path(data_root,"baydelta.tif"))
# Read the pre-processed raster file from Ryan since the conversion takes 2 hours 
raster_baydelta_tr <- readRDS(file.path(data_root,"raster_baydelta_tr.rds"))

# Enter facilities data (long and lat)
# Added Cache Slough Complex point to triangulate positoon 2/5/24
facility_data<- data.frame(Latitude=c(37.825957,37.816984,38.232488), Longitude= c(-121.595535,-121.556881,-121.674608), Facility=c("Skinner","TFCF","Cache"))

facility_data_sf <- facility_data %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))

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

# Subset to just November-April
df_subset <- df %>% mutate(WY=ifelse(month(Date)>9,year(Date)+1,year(Date))) %>% filter(month(Date) %in% c(11,12,1:4)) %>%
  # Filter out NA latitude and longitude
  filter(!is.na(Longitude)&!is.na(Latitude))

# Filter out juveniles and replicate row by fish count
df_subset_rep <- df_subset %>% filter(Length>=50)

df_subset_rep<-df_subset_rep[rep(row.names(df_subset_rep),df_subset_rep$Count),1:28] %>% mutate(Count=1)

# Extract stations at which Delta Smelt were observed (for checking on maps later)
df_subset_rep_station <- df_subset_rep %>% group_by(Station) %>% 
  summarise(Longitude=median(Longitude),Latitude=median(Latitude))

# Fix station 20-34-LSSC02 due to site potentially being too close to shore and giving infinite value to Tracy Fish Facility
df_subset_rep_station$Latitude[df_subset_rep_station$Station == "20-34-LSSC02"] <- 38.374748
df_subset_rep_station$Longitude[df_subset_rep_station$Station == "20-34-LSSC02"] <- -121.629602

# Station 208 was also a problem, where distance from 208 to Skinner was infinite
df_subset_rep %>% filter(Station=="208")
# Upon further inspection, this seems erroneous and probably a mis-identification of Longfin Smelt given that
# Salinity at surface was almost seawater (32 ppt)
# Upon inspection of map, it also seems to be incorrect, given that it is the only site past San Pablo Bay where Delta Smelt
# were ever recorded

# Extract stations at which Delta Smelt were observed and make them spatial
df_subset_rep_sf <- df_subset_rep_station %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(df_subset_rep_sf)

# Remove site 208
df_subset_rep <- df_subset_rep %>% filter(Station != "208")
df_subset_rep_sf <- df_subset_rep_sf %>% filter(Station != "208")

# Calculate all possible distances with least cost analysis
dist_facil <- costDistance(raster_baydelta_tr,
                          fromCoords = as(as_Spatial(df_subset_rep_sf), "SpatialPoints"),
                          toCoords = as(as_Spatial(facility_data_sf), "SpatialPoints"))
rownames(dist_facil)<-as.vector(df_subset_rep_sf$Station)
colnames(dist_facil)<-as.vector(facility_data_sf$Facility)

# Turn back to dataframe
Least_Cost_data_wild <-data.frame(dist_meters = as.vector(as.matrix(dist_facil)),
                                  comparison = as.vector(outer(Y=colnames(as.matrix(dist_facil)) ,
                                                               X=rownames(as.matrix(dist_facil)) ,
                                                               paste, sep="xxx"))) %>%
  mutate(Station=sub("xxx.*","", comparison),Facility=sub(".*xxx","", comparison),dist_km=dist_meters/1000) %>%
  dplyr::select(Station,Facility,dist_km)

# Export dataset out
write.csv(df_subset_rep_station,file=file.path(output_root,"wild_fish_station_data.csv"),row.names = F)
write.csv(Least_Cost_data_wild,file=file.path(output_root,"wild_fish_leastcost_distance_from_salvage_facilities.csv"),row.names = F)

# Reload data if not re-running the distance analysis
Least_Cost_data_wild<- read.csv(file=file.path(output_root,"wild_fish_leastcost_distance_from_salvage_facilities.csv"))

#######
# Start with Skinner facility
#Least_Cost_data_wild_avg <- Least_Cost_data_wild %>% group_by(Station) %>% summarise(dist_km_facility=mean(dist_km))
Least_Cost_data_wild_swp <- Least_Cost_data_wild %>% filter(Facility=="Skinner") %>% select(-Facility)
Least_Cost_data_wild_cvp <- Least_Cost_data_wild %>% filter(Facility=="TFCF") %>% select(-Facility)
Least_Cost_data_wild_cache <- Least_Cost_data_wild %>% filter(Facility=="Cache") %>% select(-Facility)

# Add skinner distance data back to catch data
df_subset_rep_swp <- df_subset_rep %>% left_join(Least_Cost_data_wild_swp)

# Summarize by month
df_subset_dist_swp <- df_subset_rep_swp %>% mutate(Month=month(Date)) %>% group_by(WY,Month) %>% 
  summarise(Date_median=median(Date),dist_km_facility=mean(dist_km),dist_km_sd=sd(dist_km,na.rm=T),Count=sum(Count))

# Add skinner distance data back to catch data
df_subset_rep_cvp <- df_subset_rep %>% left_join(Least_Cost_data_wild_cvp)

# Summarize by month
df_subset_dist_cvp <- df_subset_rep_cvp %>% mutate(Month=month(Date)) %>% group_by(WY,Month) %>% 
  summarise(Date_median=median(Date),dist_km_facility=mean(dist_km),dist_km_sd=sd(dist_km,na.rm=T),Count=sum(Count))

# Add cache distance data back to catch data
df_subset_rep_cache <- df_subset_rep %>% left_join(Least_Cost_data_wild_cache)

# Summarize by month
df_subset_dist_cache <- df_subset_rep_cache %>% mutate(Month=month(Date)) %>% group_by(WY,Month) %>% 
  summarise(Date_median=median(Date),dist_km_facility=mean(dist_km),dist_km_sd=sd(dist_km,na.rm=T),Count=sum(Count))

# Export out as csv
write.csv(df_subset_dist_swp,file=file.path(output_root,"wild_fish_average_distance_from_skinner.csv"),row.names = F)
write.csv(df_subset_dist_cvp,file=file.path(output_root,"wild_fish_average_distance_from_TFCF.csv"),row.names = F)
write.csv(df_subset_dist_cache,file=file.path(output_root,"wild_fish_average_distance_from_Cache.csv"),row.names = F)

##### Extra figures for data exploration
# Filter out just those with sample size >20
df_subset_dist <- df_subset_rep_swp %>% filter(Count>=20)

ggplot(data=df_subset_dist) + geom_point(aes(x=dist_km_facility,y=dist_km_sd))
summary(lm(data=df_subset_dist,dist_km_sd~dist_km_facility))
#Adjusted R-squared:  0.4213 

ggplot(data=df_subset_dist) + geom_boxplot(aes(x=as.factor(Month),y=dist_km_facility))

tiff(filename=file.path(output_root,"Figure_Distance_Month.tiff"),
     type="cairo",
     units="in", 
     width=6, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
ggplot(data=df_subset_dist) + geom_boxplot(aes(x=as.factor(Month),y=dist_km_facility))
dev.off()

ggplot(data=df_subset_dist) + geom_point(aes(x=Count,y=dist_km_facility))

tiff(filename=file.path(output_root,"Figure_Distance_Samplesize.tiff"),
     type="cairo",
     units="in", 
     width=6, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
ggplot(data=df_subset_dist) + geom_point(aes(x=Count,y=dist_km_facility))
dev.off()


ggplot(data=df_subset_dist) + geom_boxplot(aes(x=as.factor(WY),y=dist_km_facility))+
  theme(axis.text.x = element_text(angle = 45))

tiff(filename=file.path(output_root,"Figure_Distance_WY.tiff"),
     type="cairo",
     units="in", 
     width=12, #10*1, 
     height=5, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
ggplot(data=df_subset_dist) + geom_boxplot(aes(x=as.factor(WY),y=dist_km_facility))+
  theme(axis.text.x = element_text(angle = 45))
dev.off()

# Filter out just those with sample size >20
df_subset_dist <- df_subset_dist_cache %>% filter(Count>=20)

ggplot(data=df_subset_dist) + geom_point(aes(x=dist_km_facility,y=dist_km_sd))
summary(lm(data=df_subset_dist,dist_km_sd~dist_km_facility))
#Adjusted R-squared:  0.4213 
