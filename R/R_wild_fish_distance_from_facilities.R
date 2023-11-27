# Purpose: Prepare swim distance data for analysis for year 1
# Author: Brian Mahardja
# Date: 2023-11-22

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


# Read the TIF file just to get the crs
raster_baydelta<-raster(file.path(data_root,"baydelta.tif"))
# Read the pre-processed raster file from Ryan since the conversion takes 2 hours 
raster_baydelta_tr <- readRDS(file.path(data_root,"raster_baydelta_tr.rds"))


# Facilities data

facility_data<- data.frame(Latitude=c(37.825957,37.816984), Longitude= c(-121.595535,-121.556881), Facility=c("Skinner","TFCF"))

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
  filter(!is.na(Longitude)&!is.na(Latitude))

df_subset_rep <- df_subset %>% filter(Length>=50)

df_subset_rep<-df_subset_rep[rep(row.names(df_subset_rep),df_subset_rep$Count),1:28] %>% mutate(Count=1)

df_subset_rep_station <- df_subset_rep %>% group_by(Station) %>% 
  summarise(Longitude=median(Longitude),Latitude=median(Latitude))

df_subset_rep_sf <- df_subset_rep %>% group_by(Station) %>% 
  summarise(Longitude=median(Longitude),Latitude=median(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(df_subset_rep_sf)

# Calculate distances with least cost analysis
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