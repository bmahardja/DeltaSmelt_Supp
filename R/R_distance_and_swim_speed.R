# Purpose: Prepare swim distance data for analysis
# Author: Brian Mahardja
# Date: 2023-09-22

# Set working directory
root <- "~/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

# Load Package Libraries
library(gdistance)
library(raster)
library(sf)
library(sp)
#library(rgdal)
library(tidyverse)
library(readxl)
library(ggspatial)
library(lubridate)
library(viridis)
library(readr)
#library(rgeos)
library(deltafish)
library(mapview)

# Read the TIF file so that we can get the proper crs
raster_baydelta<-raster(file.path(data_root,"baydelta.tif"))

# Read the pre-processed raster file from Ryan since the conversion takes 2 hours 
raster_baydelta_tr <- readRDS(file.path(data_root,"raster_baydelta_tr.rds"))

# Add release site and number information for 2022
release_data_WY22 <- read_csv(file.path(data_root,"2021-2022_Year1ReleaseEventData.csv")) %>%
  rename(ReleaseEvent='Release Event') %>% group_by(ReleaseEvent,FirstDayRelease) %>%
  summarise(Total=sum(Total)) %>% mutate(FirstDayRelease=as.Date(FirstDayRelease,"%m/%d/%Y"))

release_stations_WY22 <- read_csv(file.path(data_root,"2021-2022_release_locations.csv")) %>%
  rename(ReleaseEvent='Release Event',Longitude=ReleaseLongitude,Latitude=ReleaseLatitude) %>% group_by(ReleaseEvent) %>%
  summarise(Latitude=median(Latitude,na.rm=T),Longitude=median(Longitude,na.rm=T))

release_info_WY22 <- left_join(release_data_WY22,release_stations_WY22)

releases_WY22_sf <- release_info_WY22 %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),
         ReleaseEvent=as.numeric(ReleaseEvent)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(releases_WY22_sf)

# Add release site and number information for 2023
release_data_WY23 <- read_csv(file.path(data_root,"2022-2023_ReleaseInfo.csv")) %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),ReleaseEvent=as.character(ReleaseEvent),FirstDayRelease=as.Date(FirstDayRelease, "%m/%d/%Y"))

releases_WY23_sf <- release_data_WY23  %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(releases_WY23_sf)


# Add release site and number information for 2024
release_data_WY24 <- read_excel(file.path(data_root,"2023-2024_release_locations_dates.xlsx")) %>%
  mutate(Latitude = as.numeric(Latitude),
         Longitude = as.numeric(Longitude),FirstDayRelease=as.Date(Date, "%m/%d/%Y")) %>%
  #Change Release Event name to match the catch data
  mutate(ReleaseEvent = paste("BY2023",Event,sep=" "))

releases_WY24_sf <- release_data_WY24  %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(releases_WY24_sf)

# Load running Delta Smelt catch data from Denise, downloaded 4/5/2024

catch_data <- read_excel(file.path(data_root,"Running Delta Smelt Catch_2024-04-05.xlsx"),sheet="Delta Smelt Catch Data") 
###
catch_t<-catch_data %>% mutate(Year=year(SampleDate),Month=month(SampleDate),Count=1) %>%
  filter(LifeStage %in% c("Adult")) %>%
  group_by(Year,Month) %>% summarise(Count=sum(Count))

# Remove unmarked fish and subset for each WY
catch_data_subset <- catch_data %>% filter(MarkCode != "None") %>% mutate(SampleDate=as.Date(SampleDate)) %>%
  rename(Longitude=LongitudeStart,Latitude=LatitudeStart)

# Change latitude and longitude for TFCF and Skinner because they're currently on land

catch_data_subset$Latitude[catch_data_subset$Survey == "Skinner"] <- 37.825957
catch_data_subset$Longitude[catch_data_subset$Survey == "Skinner"] <- -121.595535

catch_data_subset$Latitude[catch_data_subset$Survey == "TFCF"] <- 37.816984
catch_data_subset$Longitude[catch_data_subset$Survey == "TFCF"] <- -121.556881

# Subset just WY22 because their origin is unknown
catch_data_WY22 <- catch_data_subset %>% filter(SampleDate < as.Date("2022-05-01")) %>%
  mutate()

# Add unique ID to unnamed fish in WY22
catch_data_WY22$SpecialStudyID[catch_data_WY22$StationCode == "22-20-RV01"] <- "FCCL_01"
catch_data_WY22$SpecialStudyID[catch_data_WY22$StationCode == "22-20-LSR02"&catch_data_WY22$SpecialStudy == "FCCL Broodstock"] <- "FCCL_02"
catch_data_WY22$SpecialStudyID[catch_data_WY22$Survey == "CDFW Bay Study"] <- "BS_01"

catch_data_WY22$SpecialStudyID <- paste0("WY22_", catch_data_WY22$SpecialStudyID, sep="")


catch_data_WY22_sf <- catch_data_WY22 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(catch_data_WY22_sf)

# Add unique ID to unnamed fish in WY23
catch_data_WY23 <- catch_data_subset %>% filter(SampleDate > as.Date("2022-05-01")&SampleDate < as.Date("2023-08-01"))

catch_data_WY23$SpecialStudyID[catch_data_WY23$StationCode == "Sherman Island"] <- "FCCL_01"

catch_data_WY23$SpecialStudyID <- paste0("WY23_", catch_data_WY23$SpecialStudyID, sep="")
catch_data_WY23$ReleaseEvent <- gsub(" ", "", catch_data_WY23$ReleaseEvent, fixed = TRUE)

catch_data_WY23_sf <- catch_data_WY23 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(catch_data_WY23_sf)

# Add unique ID to unnamed fish in WY24
catch_data_WY24 <- catch_data_subset %>% filter(SampleDate > as.Date("2023-10-01")&SampleDate < as.Date("2024-08-01"))

catch_data_WY24 <- catch_data_WY24 %>% mutate(newID=paste("Salvage",row_number(),sep="_")) %>%
  mutate(SpecialStudyID=ifelse(is.na(SpecialStudyID),newID,SpecialStudyID))

catch_data_WY24$newID<-NULL
catch_data_WY24$SpecialStudyID <- paste0("WY24_", catch_data_WY24$SpecialStudyID, sep="")

catch_data_WY24_sf <- catch_data_WY24 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(raster_baydelta))
mapview(catch_data_WY24_sf)

#########################
# Calculate 2022 distances with least cost analysis
dist_WY22 <- costDistance(raster_baydelta_tr,
                          fromCoords = as(as_Spatial(releases_WY22_sf), "SpatialPoints"),
                          toCoords = as(as_Spatial(catch_data_WY22_sf), "SpatialPoints"))

rownames(dist_WY22)<-as.vector(releases_WY22_sf$ReleaseEvent)
colnames(dist_WY22)<-as.vector(catch_data_WY22_sf$SpecialStudyID)

# Turn back to dataframe
Least_Cost_data_2022 <-data.frame(dist_meters = as.vector(as.matrix(dist_WY22)),
                                  comparison = as.vector(outer(Y=colnames(as.matrix(dist_WY22)) ,
                                                               X=rownames(as.matrix(dist_WY22)) ,
                                                               paste, sep="xxx"))) %>%
  mutate(ReleaseEvent=sub("xxx.*","", comparison),SpecialStudyID=sub(".*xxx","", comparison),dist_km=dist_meters/1000) %>%
  dplyr::select(SpecialStudyID,ReleaseEvent,dist_km) %>% mutate(ReleaseEvent=as.numeric(ReleaseEvent))

dist_data_2022 <- left_join(Least_Cost_data_2022,catch_data_WY22 %>% select(-ReleaseEvent,-ReleaseDate,-ReleaseSite,-ReleaseMethod), by="SpecialStudyID") %>%
  left_join(release_info_WY22 %>% rename(ReleaseLongitude=Longitude,ReleaseLatitude=Latitude)) %>%
  #Remove occasions where fish were caught before release date
  dplyr::filter(SampleDate>=FirstDayRelease) %>%
  mutate(DaySinceRelease=as.numeric(SampleDate-FirstDayRelease)) %>%
  mutate(swim_speed=dist_km/DaySinceRelease)


# Calculate 2023 distances with least cost analysis
dist_WY23 <- costDistance(raster_baydelta_tr,
                          fromCoords = as(as_Spatial(releases_WY23_sf), "SpatialPoints"),
                          toCoords = as(as_Spatial(catch_data_WY23_sf), "SpatialPoints"))

rownames(dist_WY23)<-as.vector(releases_WY23_sf$ReleaseEvent)
colnames(dist_WY23)<-as.vector(catch_data_WY23_sf$SpecialStudyID)

# Turn back to dataframe
Least_Cost_data_2023 <-data.frame(dist_meters = as.vector(as.matrix(dist_WY23)),
                                  comparison = as.vector(outer(Y=colnames(as.matrix(dist_WY23)) ,
                                                               X=rownames(as.matrix(dist_WY23)) ,
                                                               paste, sep="xxx"))) %>%
  mutate(ReleaseEvent=sub("xxx.*","", comparison),SpecialStudyID=sub(".*xxx","", comparison),dist_km=dist_meters/1000) %>%
  dplyr::select(SpecialStudyID,ReleaseEvent,dist_km) 


dist_data_2023 <- left_join(catch_data_WY23 %>% select(-ReleaseDate,-ReleaseSite,-ReleaseMethod),Least_Cost_data_2023, by=c("SpecialStudyID","ReleaseEvent")) %>%
  left_join(release_data_WY23 %>% rename(ReleaseLongitude=Longitude,ReleaseLatitude=Latitude)) %>%
  mutate(DaySinceRelease=as.numeric(SampleDate-FirstDayRelease)) %>%
  mutate(swim_speed=dist_km/DaySinceRelease)


# Calculate 2024 distances with least cost analysis
dist_WY24 <- costDistance(raster_baydelta_tr,
                          fromCoords = as(as_Spatial(releases_WY24_sf), "SpatialPoints"),
                          toCoords = as(as_Spatial(catch_data_WY24_sf), "SpatialPoints"))

rownames(dist_WY24)<-as.vector(releases_WY24_sf$ReleaseEvent)
colnames(dist_WY24)<-as.vector(catch_data_WY24_sf$SpecialStudyID)


# Turn back to dataframe
Least_Cost_data_2024 <-data.frame(dist_meters = as.vector(as.matrix(dist_WY24)),
                                  comparison = as.vector(outer(Y=colnames(as.matrix(dist_WY24)) ,
                                                               X=rownames(as.matrix(dist_WY24)) ,
                                                               paste, sep="xxx"))) %>%
  mutate(ReleaseEvent=sub("xxx.*","", comparison),SpecialStudyID=sub(".*xxx","", comparison),dist_km=dist_meters/1000) %>%
  dplyr::select(SpecialStudyID,ReleaseEvent,dist_km) 


dist_data_2024 <- left_join(catch_data_WY24 %>% dplyr::select(-ReleaseDate,-ReleaseSite,-ReleaseMethod),Least_Cost_data_2024) %>%
  left_join(release_data_WY24 %>% rename(ReleaseLongitude=Longitude,ReleaseLatitude=Latitude)) %>%
  mutate(DaySinceRelease=as.numeric(SampleDate-FirstDayRelease)) %>%
  mutate(swim_speed=dist_km/DaySinceRelease) 

# Export dataset out
write.csv(dist_data_2022,file=file.path(output_root,"WY2022_distance_from_release_site.csv"),row.names = F)
write.csv(dist_data_2023,file=file.path(output_root,"WY2023_distance_from_release_site.csv"),row.names = F)
write.csv(dist_data_2024,file=file.path(output_root,"WY2024_distance_from_release_site.csv"),row.names = F)


###
plot(dist_data_2022$swim_speed~dist_data_2022$DaySinceRelease)
plot(dist_data_2023$swim_speed~dist_data_2023$DaySinceRelease)
plot(dist_data_2024$swim_speed~dist_data_2024$DaySinceRelease)

plot(dist_data_2022$swim_speed~dist_data_2022$SampleDate)
plot(dist_data_2023$swim_speed~dist_data_2023$SampleDate)
plot(dist_data_2024$swim_speed~dist_data_2024$SampleDate)

plot(dist_data_2022$swim_speed~dist_data_2022$ForkLength)
plot(dist_data_2023$swim_speed~dist_data_2023$ForkLength)
plot(dist_data_2024$swim_speed~dist_data_2024$ForkLength)

plot(dist_data_2022$dist_km~dist_data_2022$DaySinceRelease)
plot(dist_data_2023$dist_km~dist_data_2023$DaySinceRelease)
plot(dist_data_2024$dist_km~dist_data_2024$DaySinceRelease)

