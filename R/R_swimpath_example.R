# Purpose: Example how to create swim paths
# Author: Brian Mahardja
# Date: 2024-01-24

# Set working directory
root <- "~/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

# Load packages
library(tidyverse)
library(readxl)
library(boot)
library(gdistance)
library(raster)
library(sf)
library(sp)
library(ggspatial)
library(deltamapr)

###################
# Load example fish data
smelt_data_join <- readRDS(file.path(output_root,"smelt_data_example.rds"))

# Read in spatial data
sp_releases <- readRDS(file.path(output_root,"spatial_release_data_example.rds"))
sp_station <- readRDS(file.path(output_root,"spatial_station_data_example.rds"))
raster_baydelta_tr <- readRDS(file.path(data_root,"raster_baydelta_tr.rds"))
raster_baydelta<-raster(file.path(data_root,"baydelta.tif"))

# Subset station to captured delta smelt
sp_station <- sp_station %>% filter(StationCode %in% unique(smelt_data_join$StationCode))

# This doesn't work
AtoB <- shortestPath(raster_baydelta_tr, sp_releases[1,], sp_station[2,], output="SpatialLines")

# But example data does seem to work with a single point to another single point
AtoB <- shortestPath(raster_baydelta_tr,c(611515.2,4217818), c(614331.5,4221416), output="SpatialLines")
AtoC <- shortestPath(raster_baydelta_tr,c(611515.2,4217818), c(626685.7, 4186284), output="SpatialLines")

# This can be plotted in sp
plot(raster_baydelta)
lines(AtoB, col="dark green", lwd=2)

# Convert line type to Sf
line1<-st_as_sf(AtoB) 
line2<-st_as_sf(AtoC) 

# Bind data 
linec<-bind_rows(line1,line2) %>% st_set_crs(st_crs(sp_station))


ggplot() + geom_sf(data=WW_Delta, fill="cadetblue1", color="cadetblue1") + 
  geom_sf(data=sp_station,color="purple")+
  geom_sf(data=sp_releases,color="red")+
  geom_sf(data=linec)

# But how do we automate this for all the fish in WY23?
