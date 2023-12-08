# Purpose: Survey data to assess swimming movement in SKT
# Author: Brian Mahardja
# Date: 2023-09-22

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
library(collapse)
library(deltamapr)

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
df_subset <- df %>% mutate(WY=ifelse(month(Date)>9,year(Date)+1,year(Date))) %>% filter(month(Date) %in% c(11:12,1:4)) %>%
  # Filter out NA salinity data
  filter(!is.na(Sal_surf))

# Plot histogram to screen off juveniles vs. adults
hist(df_subset$Length,breaks=30)
# For take, smelt is considered "adults" at 58 mm
# Histogram shows most fish in these months were from 50mm and above, so we will use that

# Multiply row by individuals
df_subset_rep<-df_subset[rep(row.names(df_subset),df_subset$Count),1:28] %>% mutate(Count=1)

# Filter out individuals under 50 mm, add Month column
df_subset_rep <- df_subset_rep %>% filter(Length>=50) %>% mutate(Month=month(Date))


# Summarize by month and WY
df_subset_month_sum <- df_subset_rep %>% group_by(WY, Month) %>%
  summarise(Count_Sal=sum(Count),Sal_surf_mean=mean(Sal_surf),Sal_surf_sd=sd(Sal_surf),
            Date_median_sal=median(Date),
            Sal_min=min(Sal_surf,na.rm=T),Sal_max=max(Sal_surf,na.rm=T),
            Survey_mostcommon_sal=fmode(as.factor(Source)))

# Filter out data with less than 20 fish
df_subset_month_sum<- df_subset_month_sum %>% filter(Count_Sal>=20)

# Export out data
write.csv(df_subset_month_sum,file=file.path(output_root,"wild_fish_salinity_data.csv"),row.names = F)

#Print figures
tiff(filename=file.path(output_root,"Figure_Salinity_mean_Month.tiff"),
     type="cairo",
     units="in", 
     width=6, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
ggplot(data=df_subset_month_sum) + geom_boxplot(aes(x=as.factor(Month),y=Sal_surf_mean))
dev.off()

tiff(filename=file.path(output_root,"Figure_Salinity_sd_Month.tiff"),
     type="cairo",
     units="in", 
     width=6, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
ggplot(data=df_subset_month_sum) + geom_boxplot(aes(x=as.factor(Month),y=Sal_surf_sd))
dev.off()

########################################################
# Summarize effort
df_subset_effort <- df <- left_join(surv, fish_smelt) %>% 
  collect() %>% mutate(WY=ifelse(month(Date)>9,year(Date)+1,year(Date))) %>% filter(month(Date) %in% c(11:12,1:4)) %>%
  filter(!is.na(Longitude)&!is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = st_crs(R_EDSM_Subregions_Mahardja))

df_subset_effort_spatial<- st_join(df_subset_effort,R_EDSM_Subregions_Mahardja)
unique(df_subset_effort_spatial$SubRegion)

df_subset_effort_spatial <- df_subset_effort_spatial %>%
  mutate(Region_new = case_when(
    SubRegion %in% c("Liberty Island","Lower Sacramento River Ship Channel","Cache Slough and Lindsey Slough",
                     "Upper Sacramento River Ship Channel","Upper Yolo Bypass") ~ "Cache Slough Complex",
    SubRegion %in% c("Grizzly Bay","Lower San Joaquin River","Lower Sacramento River","Confluence","Mid Suisun Bay","Honker Bay",
                     "Sacramento River near Rio Vista","Lower Cache Slough","West Suisun Bay",
                     "San Joaquin River at Prisoners Pt","Sacramento River near Ryde") ~ "Confluence and Suisun Bay",
    SubRegion %in% c("Suisun Marsh") ~ "Suisun Marsh",
    SubRegion %in% c("San Joaquin River at Twitchell Island","Holland Cut",
                     "Old River","Mildred Island","Lower Mokelumne River",
                     "Upper Mokelumne River","Steamboat and Miner Slough","San Joaquin River near Stockton",
                     "Middle Sacramento River","Franks Tract","Victoria Canal",
                     "Upper Sacramento River", "Disappointment Slough",
                     "Georgiana Slough","Middle River","Rock Slough and Discovery Bay",
                     "Grant Line Canal and Old River","Upper San Joaquin River" ) ~ "Interior Delta",
    SubRegion %in% c("Lower Napa River","Upper Napa River","Carquinez Strait","San Pablo Bay","San Francisco Bay",
                     "South Bay") ~ "Carquinez and Bays",
  ))

# Some dataset returned NA, this is to check what we're excluding, especially for when smelt are present
df_orphan<- df_subset_effort_spatial %>% filter(is.na(Region_new)&Count>0)
mapView(df_orphan)

# They are all very upstream on the Sac and SJ rivers, may be a misidentification with Wakasagi on the upper Sac side
# Seems safe to remove
df_subset_effort_spatial <- df_subset_effort_spatial %>% filter(!is.na(Region_new))

# Summarize by sampling occasion
df_subset_effort_spatial_sum <- df_subset_effort_spatial %>% st_drop_geometry() %>% group_by(WY, Date, Datetime, Station, SampleID, Method, Source,Region_new) %>%
  summarise(Count=sum(Count)) %>% mutate(SampleSize=1,Month=month(Date))

# Summarize by region and month
df_subset_effort_spatial_sum <- df_subset_effort_spatial_sum %>% ungroup() %>% group_by(WY,Month,Region_new) %>%
  summarise(Date_median=median(Date),SampleSize=sum(SampleSize),Survey_mostcommon=fmode(as.factor(Source))) %>% 
  mutate(Region_new = as.factor(Region_new),Year=year(Date_median)) %>%
  mutate(Year_month=floor_date(Date_median, "month"))

ggplot(data=df_subset_effort_spatial_sum, aes(Region_new, as.factor(Year_month), fill= SampleSize)) + 
  geom_tile()

ggplot(data=df_subset_effort_spatial_sum, aes(Region_new, as.factor(WY), fill= SampleSize)) + 
  geom_tile()


ggplot(data=df_subset_effort_spatial_sum %>% filter(Region_new=="Cache Slough Complex"), aes(x=as.factor(WY), y= SampleSize)) + 
  geom_line()

############
df_subset_v2 <- df_subset %>% filter(Length>=50) %>% mutate(Month=month(Date),Latitude_mult=Latitude*Count,Longitude_mult=Longitude*Count) %>% group_by(WY, Month) %>%
  summarise(Count=sum(Count),Date_median=median(Date),Latitude_sum=sum(Latitude_mult),Longitude_sum=sum(Longitude_mult),Sal_surf_sum=sum(Sal_surf),Sal_max=max(Sal_surf,na.rm=T),Sal_min=min(Sal_surf,na.rm=T)) %>%
  mutate(Latitude_c = Latitude_sum/Count, Longitude_c= Longitude_sum/Count,Sal_surf_mean=Sal_surf_sum/Count)


df_subset_v3 <- df_subset %>% filter(Length>=50)

df_subset_v3<-df_subset_v3[rep(row.names(df_subset_v3),df_subset_v3$Count),1:28]
  
df_subset_v3 <- df_subset_v3 %>% mutate(Month=month(Date),Count=1) %>% group_by(WY, Month) %>%
  summarise(Count=sum(Count),Date_median=median(Date),Sal_surf_mean=mean(Sal_surf),Sal_max=max(Sal_surf,na.rm=T),Sal_min=min(Sal_surf,na.rm=T),Sal_surf_sd=sd(Sal_surf)) %>%
  mutate(Sal_surf_cv=Sal_surf_sd/Sal_surf_mean)

hist(df_subset_v3$Count)
df_subset_v3_min<-df_subset_v3 %>% filter(Count<=100)
hist(df_subset_v3_min$Count)


ggplot(data=df_subset_v3) + geom_point(aes(x=Month,y=Sal_surf_mean))
ggplot(data=df_subset_v3) + geom_point(aes(x=Month,y=Sal_surf_sd))
ggplot(data=df_subset_v3) + geom_point(aes(x=Month,y=Sal_surf_cv))
ggplot(data=df_subset_v3) + geom_boxplot(aes(x=as.factor(Month),y=Sal_surf_mean))
ggplot(data=df_subset_v3) + geom_boxplot(aes(x=as.factor(Month),y=Sal_surf_sd))

ggplot(data=df_subset_v3) + geom_point(aes(x=Sal_surf_mean,y=Sal_surf_sd))

df_subset_v3_1 <- df_subset_v3 %>% filter(Count>=20)

ggplot(data=df_subset_v3_1) + geom_boxplot(aes(x=as.factor(Month),y=Sal_surf_mean))
ggplot(data=df_subset_v3_1) + geom_boxplot(aes(x=as.factor(Month),y=Sal_surf_sd))
ggplot(data=df_subset_v3_1) + geom_point(aes(x=Sal_surf_mean,y=Sal_surf_sd))

write.csv(df_subset_v3_1,file=file.path(output_root,"wild_fish_salinity_data.csv"),row.names = F)

df_subset_v3_2 <- df_subset_v3 %>% group_by(Month) %>% summarise(Sal_surf_mean=mean(Sal_surf_mean),Sal_surf_sd=mean(Sal_surf_sd,na.rm=T))


