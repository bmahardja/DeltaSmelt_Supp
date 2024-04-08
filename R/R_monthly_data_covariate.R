# Purpose: Monthly data covariate
# Author: Brian Mahardja
# Date: 2023-12-11

library(tidyverse)
library(rvest)
library(lubridate)
library(janitor)
library(sharpshootR)
library(deltamapr)
library(sf)
library(imputeTS)

# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

## Extract OMR data

#Code from Trinh Nguyen's OMR Data Pull 8/4/2022; fixed incorrect end_date url
calc_OMR <- function(dateStart, dateEnd = NULL, timing, extrap = F, proof = F, showQAQC = F) {
  # Do you have the two required packages?
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("zoo package is required", call. = FALSE)
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr package is required", call. = FALSE)
  }
  
  # The dateStart date has to dateStart 14 days before the specified dateStart date.
  # This is to account for NA padding during the rolling mean calculations.
  # The 13 is correct here.
  dateStartPad <- as.Date(dateStart) - lubridate::days(13)
  
  # If no dateEnd date is given, use today
  if (is.null(dateEnd)) {
    dateEnd <- lubridate::today()
  }
  
  # For hourly data, timing = "uv"; for daily it's "dv"
  if (!timing %in% c("daily", "hourly")) {
    stop("Timing variable should either be daily or hourly.")
  }
  
  if (timing %in% "daily") {
    timing <- "dv"
    columnNames <- c("agency", "site", "date", "flow", "qaqc")
  } else {
    if (timing %in% "hourly") {
      timing <- "uv"
      columnNames <- c("agency", "site", "date", "time","TZ", "flow", "qaqc")
    }
  }
  
  # Old R
  oldURL <- paste0("https://nwis.waterdata.usgs.gov/nwis/",
                   timing,
                   "?cb_72137=on&format=rdb&site_no=11313405&referred_module=sw&period=&begin_date=",
                   dateStartPad,
                   "&end_date=", 
                   dateEnd)
  
  oldR <- readr::read_table2(oldURL, 
                             col_names = columnNames,
                             skip = 31) %>%
    rename(oldR = site, oldFlow = flow)
  cat("Old River URL: \n", oldURL, "\n")
  
  # Middle R
  middleURL <- paste0("https://nwis.waterdata.usgs.gov/nwis/",
                      timing,
                      "?cb_72137=on&format=rdb&site_no=11312676&referred_module=sw&period=&begin_date=",
                      dateStartPad,
                      "&end_date=", 
                      dateEnd)
  
  middleR <- readr::read_table2(middleURL, 
                                col_names = columnNames,
                                skip = 31) %>%
    rename(middleR = site, middleFlow = flow)
  cat("Middle R URL: \n", middleURL, "\n")
  
  # Is there data for the time frame that was specified?
  if (nrow(middleR) == 0 | nrow(oldR) == 0) {
    cat("Middle R has", nrow(middleR), "rows. \n" )
    cat("Old R has", nrow(oldR), "rows. \n")
    stop("One of the gage had no data for this time frame")
  }
  
  # QAQC range
  if (isTRUE(showQAQC)) {
    
    cat("For Old River \n")
    
    oldR %>% 
      group_by(qaqc) %>% 
      slice(1, n()) %>% 
      print()
    
    cat("For Middle River \n")
    warning("Function will stop here. Disable QAQC to run further.", call. = F)
    return(middleR %>% 
             group_by(qaqc) %>% 
             slice(1, n()) %>% 
             print())
  }
  
  # Creating OMR table from the tidally filtered means
  if (timing %in% "uv") {
    omrDF <- full_join(oldR, middleR, by = c("agency", "date", "time", "TZ"))
  } else {
    if (timing %in% "dv") {
      omrDF <- full_join(oldR, middleR, by = c("agency", "date"))
    }
  }
  
  omrDF <- omrDF %>%
    transmute(date,
              week = lubridate::week(date),
              month = lubridate::month(date),
              year = lubridate::year(date),
              waterYear = year + (month > 9),
              oldFlow, middleFlow,
              # You want this to not have na.rm
              omrFlow = oldFlow + middleFlow) %>% 
    # If timing = daily, then this gives the same results
    group_by(date, week, month, year, waterYear) %>%
    summarise(oldFlow = mean(oldFlow),
              middleFlow = mean(middleFlow),
              omrFlow = mean(omrFlow),
              .groups = "drop")
  # mutate(omrFlow5 = zoo::rollmean(omrFlow, 5, na.pad = T, align = "right"),
  #        omrFlow14 = zoo::rollmean(omrFlow, 14, na.pad = T, align = "right")) %>%
  # # Filtering out the padded dates
  # filter(date >= dateStart)
  
  # Did you want to extrapolate missing values from the other?
  if (extrap) {
    if (sum(is.na(omrDF$oldFlow), is.na(omrDF$middleFlow)) == 0) {
      message("Extrapolation wasn't needed. Full data for both rivers.")
    }
    # Extrapolation occuring on DAILY levels
    oldLM <- lm(oldFlow ~ middleFlow, data = omrDF)
    middleLM <- lm(middleFlow ~ oldFlow, data = omrDF)
    # Need to predict at only the NA values.
    omrDF <- omrDF %>%
      mutate(oldRExtrap = case_when(is.na(oldFlow) ~ predict(oldLM, omrDF),
                                    TRUE ~ oldFlow),
             middleRExtrap = case_when(is.na(middleFlow) ~ predict(middleLM, omrDF),
                                       TRUE ~ middleFlow)) %>%
      mutate(omrFlowExtrap = oldRExtrap + middleRExtrap) %>%
      ungroup()
    # mutate(omrFlow5Extrap = zoo::rollmean(omrFlowExtrap, 5, na.pad = T, align = "right"),
    #        omrFlow14Extrap = zoo::rollmean(omrFlowExtrap, 14, na.pad = T, align = "right"))
    
    if (proof) {
      proofPlot <- ggplot2::ggplot(omrDF, ggplot2::aes(oldFlow, middleFlow)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", formula = "y ~ x") +
        ggplot2::labs(title = paste0("Regression of daily Middle ~ Old flow", ", adjR = ", 
                                     round(summary(middleLM)$adj.r.squared, 2)),
                      subtitle = paste0("From ", min(omrDF$date), " to ", max(omrDF$date)),
                      x = "Old River Flow (cfs)", y = "Middle River Flow (cfs)") +
        ggplot2::theme_classic(base_size = 18)
      
      print(proofPlot)
    }
  }
  
  return(omrDF)
}

###############
#Pull OMR Data

data_OMR <- calc_OMR(dateStart = "1987-10-01",dateEnd = "2023-12-01", timing = "daily", extrap = T, proof = T)
#Regression of daily middle vs old river flow - adj R = 0.96
#data from 1987-09-18 to 2023-12-01

# Export out OMR data
write.csv(data_OMR,file=file.path(output_root,"rawdata_OMR.csv"),row.names = F)
# Load data if necessary
#data_OMR <- read.csv(file=file.path(output_root,"rawdata_OMR.csv"))

# Summarize OMR data into monthly timestep
data_OMR_month <- data_OMR %>% rename(Month=month, WY=waterYear) %>%
  group_by(Month,WY) %>% summarise(mean_OMR = mean(omrFlowExtrap))


###############
# Pull Dayflow data

dayflow_post97<-read.csv(file=file.path(data_root,"dayflow-results-1997-2020.csv")) %>%
  mutate(Date=as.Date(Date,"%m/%d/%Y"))

dayflow_pre97<-read.csv(file=file.path(data_root,"dayflowCalculations_pre97.csv"))  %>%
  mutate(Date=as.Date(Date,"%m/%d/%Y")) %>% rename(Month=Mo)

dayflow_merged<- bind_rows(dayflow_post97 %>% dplyr::select(Year,Month,Date,SAC,SJR,TOT,YOLO,EXPORTS,OUT,X2),
                           dayflow_pre97 %>% dplyr::select(Year,Month,Date,SAC,SJR,TOT,YOLO,EXPORTS,OUT,X2))

# Summarize by month
dayflow_month <- dayflow_merged %>% mutate(WY= ifelse(Month>9,Year+1,Year)) %>% group_by(WY,Year,Month) %>%
  select(-Date) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) 

#############
# Combine OMR and Dayflow data
data_month_combined <- dayflow_month %>% left_join(data_OMR_month)

# Create model to predict OMR based on SJR and Export flows
summary(lm(data=data_month_combined, mean_OMR ~ SJR+EXPORTS))                                               
OMR_model <- lm(data=data_month_combined, mean_OMR ~ SJR+EXPORTS)

# Fill OMR data gaps with OMR model above
data_month_combined <- data_month_combined %>% mutate(
  mean_OMR = ifelse(is.na(mean_OMR),
               predict(OMR_model), mean_OMR)
)


############
# Order data by date and month, and add lag value (data from previous month)
data_month_combined<- data_month_combined %>% mutate(DateMonth=as.Date(paste(Year,Month,"01",sep="-")))  %>% ungroup()
data_month_combined<- data_month_combined[order(data_month_combined$DateMonth),] %>%
  mutate(SAC_1 = lag(SAC,1),
         SJR_1 = lag(SJR,1),
         TOT_1 = lag(TOT,1),
         YOLO_1 = lag(YOLO,1),
         EXPORTS_1 = lag(EXPORTS,1),
         OUT_1 = lag(OUT,1),
         X2_1 = lag(X2,1))

##########
# Pull Freeport data from CDEC

FPT_daily <- CDECquery(id='FPT', sensor=20, interval='D', start=as.Date("1976-10-01"), end=as.Date("2024-04-03"))

# Data contains some zeroes
#ggplot(data=FPT_daily) + geom_point(aes(x=datetime,y=value))

# Due to potentially erroneous data from CDEC gage, whenever CDEC FPT data has a difference of more than 10,000 cfs,
# we use the DAYFLOW data
FPT_data<- FPT_daily %>% mutate(Date=as.Date(datetime),Freeport_flow_CDEC=as.numeric(value)) %>%
  left_join((dayflow_merged %>%  dplyr::select(Date,SAC))) %>%
  mutate(Freeport_difference=SAC-Freeport_flow_CDEC) %>%
  mutate(Freeport_flow_final = ifelse(!is.na(SAC)&(Freeport_flow_CDEC==0|Freeport_difference>=10000|Freeport_difference<=-10000),SAC,Freeport_flow_CDEC)) 


# Calculate 3-day rolling average
FPT_data_rollavg <- FPT_data %>% dplyr::select(water_year,Date,Freeport_flow_final) %>%
  mutate(day_average_Freeport_flow = zoo::rollmean(Freeport_flow_final,3,align = "right", na.pad = TRUE)) %>%
  mutate(above_25000=ifelse(day_average_Freeport_flow>25000,1,0))

# Calculate when first flush occurred based on the 2019 PA
FPT_data_firstflush <- FPT_data_rollavg %>% group_by(water_year) %>% filter(month(Date) %in% c(11,12,1,2,3,4)) %>%
  filter(above_25000>0) %>% summarise(FirstFlush_firstDate=min(Date),FirstFlush_medianDate=median(Date)) %>%
  rename(WY=water_year)


# Add data from DAYFLOW to fill since there was no data from WY 1976 but it exists in the fish data
FPT_data_add <- dayflow_merged %>% mutate(WY= ifelse(Month>9,Year+1,Year)) %>%
  mutate(day_average_Freeport_flow = zoo::rollmean(SAC,3,align = "right", na.pad = TRUE)) %>%
  mutate(above_25000=ifelse(day_average_Freeport_flow>25000,1,0))%>%
  filter(above_25000>0) %>% group_by(WY) %>% summarise(FirstFlush_firstDate=min(Date),FirstFlush_medianDate=median(Date)) %>%
  filter(WY==1976)

FPT_data_firstflush<-rbind(FPT_data_firstflush,FPT_data_add)

# Export 3-day rolling average Freeport flow data for hatchery fish purposes
write.csv(FPT_data_rollavg,file=file.path(output_root,"Freeport_FirstFlush_data.csv"),row.names = F)

#########
# Pull water quality data for secchi depth
# https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=731


inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/731/7/6c5f35b1d316e39c8de0bfadfb3c9692" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "Source",     
                 "Station",     
                 "Latitude",     
                 "Longitude",     
                 "Field_coords",     
                 "Date",     
                 "Datetime",     
                 "Depth",     
                 "Sample_depth_surface",     
                 "Sample_depth_nutr_surface",     
                 "Sample_depth_bottom",     
                 "Tide",     
                 "Temperature",     
                 "Temperature_bottom",     
                 "Conductivity",     
                 "Conductivity_bottom",     
                 "Salinity",     
                 "Salinity_bottom",     
                 "Secchi",     
                 "Secchi_estimated",     
                 "TurbidityNTU",     
                 "TurbidityNTU_bottom",     
                 "TurbidityFNU",     
                 "TurbidityFNU_bottom",     
                 "DissolvedOxygen",     
                 "DissolvedOxygen_bottom",     
                 "DissolvedOxygenPercent",     
                 "DissolvedOxygenPercent_bottom",     
                 "pH",     
                 "pH_bottom",     
                 "Microcystis",     
                 "Chlorophyll_Sign",     
                 "Chlorophyll",     
                 "Pheophytin_Sign",     
                 "Pheophytin",     
                 "TotAmmonia_Sign",     
                 "TotAmmonia",     
                 "DissAmmonia_Sign",     
                 "DissAmmonia",     
                 "DissNitrateNitrite_Sign",     
                 "DissNitrateNitrite",     
                 "TotPhos_Sign",     
                 "TotPhos",     
                 "DissOrthophos_Sign",     
                 "DissOrthophos",     
                 "TON_Sign",     
                 "TON",     
                 "DON_Sign",     
                 "DON",     
                 "TKN_Sign",     
                 "TKN",     
                 "TotAlkalinity_Sign",     
                 "TotAlkalinity",     
                 "DissBromide_Sign",     
                 "DissBromide",     
                 "DissCalcium_Sign",     
                 "DissCalcium",     
                 "TotChloride_Sign",     
                 "TotChloride",     
                 "DissChloride_Sign",     
                 "DissChloride",     
                 "DissSilica_Sign",     
                 "DissSilica",     
                 "TOC_Sign",     
                 "TOC",     
                 "DOC_Sign",     
                 "DOC",     
                 "TDS_Sign",     
                 "TDS",     
                 "TSS_Sign",     
                 "TSS",     
                 "VSS_Sign",     
                 "VSS",     
                 "Notes"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$Source)!="factor") dt1$Source<- as.factor(dt1$Source)
if (class(dt1$Station)!="factor") dt1$Station<- as.factor(dt1$Station)
if (class(dt1$Latitude)=="factor") dt1$Latitude <-as.numeric(levels(dt1$Latitude))[as.integer(dt1$Latitude) ]               
if (class(dt1$Latitude)=="character") dt1$Latitude <-as.numeric(dt1$Latitude)
if (class(dt1$Longitude)=="factor") dt1$Longitude <-as.numeric(levels(dt1$Longitude))[as.integer(dt1$Longitude) ]               
if (class(dt1$Longitude)=="character") dt1$Longitude <-as.numeric(dt1$Longitude)
if (class(dt1$Field_coords)!="factor") dt1$Field_coords<- as.factor(dt1$Field_coords)                                   
# attempting to convert dt1$Date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1Date<-as.Date(dt1$Date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1Date) == length(tmp1Date[!is.na(tmp1Date)])){dt1$Date <- tmp1Date } else {print("Date conversion failed for dt1$Date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1Date)                                    
# attempting to convert dt1$Datetime dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d %H:%M:%S" 
tmp1Datetime<-as.POSIXct(dt1$Datetime,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1Datetime) == length(tmp1Datetime[!is.na(tmp1Datetime)])){dt1$Datetime <- tmp1Datetime } else {print("Date conversion failed for dt1$Datetime. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1Datetime) 
if (class(dt1$Depth)=="factor") dt1$Depth <-as.numeric(levels(dt1$Depth))[as.integer(dt1$Depth) ]               
if (class(dt1$Depth)=="character") dt1$Depth <-as.numeric(dt1$Depth)
if (class(dt1$Sample_depth_surface)=="factor") dt1$Sample_depth_surface <-as.numeric(levels(dt1$Sample_depth_surface))[as.integer(dt1$Sample_depth_surface) ]               
if (class(dt1$Sample_depth_surface)=="character") dt1$Sample_depth_surface <-as.numeric(dt1$Sample_depth_surface)
if (class(dt1$Sample_depth_nutr_surface)=="factor") dt1$Sample_depth_nutr_surface <-as.numeric(levels(dt1$Sample_depth_nutr_surface))[as.integer(dt1$Sample_depth_nutr_surface) ]               
if (class(dt1$Sample_depth_nutr_surface)=="character") dt1$Sample_depth_nutr_surface <-as.numeric(dt1$Sample_depth_nutr_surface)
if (class(dt1$Sample_depth_bottom)=="factor") dt1$Sample_depth_bottom <-as.numeric(levels(dt1$Sample_depth_bottom))[as.integer(dt1$Sample_depth_bottom) ]               
if (class(dt1$Sample_depth_bottom)=="character") dt1$Sample_depth_bottom <-as.numeric(dt1$Sample_depth_bottom)
if (class(dt1$Tide)!="factor") dt1$Tide<- as.factor(dt1$Tide)
if (class(dt1$Temperature)=="factor") dt1$Temperature <-as.numeric(levels(dt1$Temperature))[as.integer(dt1$Temperature) ]               
if (class(dt1$Temperature)=="character") dt1$Temperature <-as.numeric(dt1$Temperature)
if (class(dt1$Temperature_bottom)=="factor") dt1$Temperature_bottom <-as.numeric(levels(dt1$Temperature_bottom))[as.integer(dt1$Temperature_bottom) ]               
if (class(dt1$Temperature_bottom)=="character") dt1$Temperature_bottom <-as.numeric(dt1$Temperature_bottom)
if (class(dt1$Conductivity)=="factor") dt1$Conductivity <-as.numeric(levels(dt1$Conductivity))[as.integer(dt1$Conductivity) ]               
if (class(dt1$Conductivity)=="character") dt1$Conductivity <-as.numeric(dt1$Conductivity)
if (class(dt1$Conductivity_bottom)=="factor") dt1$Conductivity_bottom <-as.numeric(levels(dt1$Conductivity_bottom))[as.integer(dt1$Conductivity_bottom) ]               
if (class(dt1$Conductivity_bottom)=="character") dt1$Conductivity_bottom <-as.numeric(dt1$Conductivity_bottom)
if (class(dt1$Salinity)=="factor") dt1$Salinity <-as.numeric(levels(dt1$Salinity))[as.integer(dt1$Salinity) ]               
if (class(dt1$Salinity)=="character") dt1$Salinity <-as.numeric(dt1$Salinity)
if (class(dt1$Salinity_bottom)=="factor") dt1$Salinity_bottom <-as.numeric(levels(dt1$Salinity_bottom))[as.integer(dt1$Salinity_bottom) ]               
if (class(dt1$Salinity_bottom)=="character") dt1$Salinity_bottom <-as.numeric(dt1$Salinity_bottom)
if (class(dt1$Secchi)=="factor") dt1$Secchi <-as.numeric(levels(dt1$Secchi))[as.integer(dt1$Secchi) ]               
if (class(dt1$Secchi)=="character") dt1$Secchi <-as.numeric(dt1$Secchi)
if (class(dt1$Secchi_estimated)!="factor") dt1$Secchi_estimated<- as.factor(dt1$Secchi_estimated)
if (class(dt1$TurbidityNTU)=="factor") dt1$TurbidityNTU <-as.numeric(levels(dt1$TurbidityNTU))[as.integer(dt1$TurbidityNTU) ]               
if (class(dt1$TurbidityNTU)=="character") dt1$TurbidityNTU <-as.numeric(dt1$TurbidityNTU)
if (class(dt1$TurbidityNTU_bottom)=="factor") dt1$TurbidityNTU_bottom <-as.numeric(levels(dt1$TurbidityNTU_bottom))[as.integer(dt1$TurbidityNTU_bottom) ]               
if (class(dt1$TurbidityNTU_bottom)=="character") dt1$TurbidityNTU_bottom <-as.numeric(dt1$TurbidityNTU_bottom)
if (class(dt1$TurbidityFNU)=="factor") dt1$TurbidityFNU <-as.numeric(levels(dt1$TurbidityFNU))[as.integer(dt1$TurbidityFNU) ]               
if (class(dt1$TurbidityFNU)=="character") dt1$TurbidityFNU <-as.numeric(dt1$TurbidityFNU)
if (class(dt1$TurbidityFNU_bottom)=="factor") dt1$TurbidityFNU_bottom <-as.numeric(levels(dt1$TurbidityFNU_bottom))[as.integer(dt1$TurbidityFNU_bottom) ]               
if (class(dt1$TurbidityFNU_bottom)=="character") dt1$TurbidityFNU_bottom <-as.numeric(dt1$TurbidityFNU_bottom)
if (class(dt1$DissolvedOxygen)=="factor") dt1$DissolvedOxygen <-as.numeric(levels(dt1$DissolvedOxygen))[as.integer(dt1$DissolvedOxygen) ]               
if (class(dt1$DissolvedOxygen)=="character") dt1$DissolvedOxygen <-as.numeric(dt1$DissolvedOxygen)
if (class(dt1$DissolvedOxygen_bottom)=="factor") dt1$DissolvedOxygen_bottom <-as.numeric(levels(dt1$DissolvedOxygen_bottom))[as.integer(dt1$DissolvedOxygen_bottom) ]               
if (class(dt1$DissolvedOxygen_bottom)=="character") dt1$DissolvedOxygen_bottom <-as.numeric(dt1$DissolvedOxygen_bottom)
if (class(dt1$DissolvedOxygenPercent)=="factor") dt1$DissolvedOxygenPercent <-as.numeric(levels(dt1$DissolvedOxygenPercent))[as.integer(dt1$DissolvedOxygenPercent) ]               
if (class(dt1$DissolvedOxygenPercent)=="character") dt1$DissolvedOxygenPercent <-as.numeric(dt1$DissolvedOxygenPercent)
if (class(dt1$DissolvedOxygenPercent_bottom)=="factor") dt1$DissolvedOxygenPercent_bottom <-as.numeric(levels(dt1$DissolvedOxygenPercent_bottom))[as.integer(dt1$DissolvedOxygenPercent_bottom) ]               
if (class(dt1$DissolvedOxygenPercent_bottom)=="character") dt1$DissolvedOxygenPercent_bottom <-as.numeric(dt1$DissolvedOxygenPercent_bottom)
if (class(dt1$pH)=="factor") dt1$pH <-as.numeric(levels(dt1$pH))[as.integer(dt1$pH) ]               
if (class(dt1$pH)=="character") dt1$pH <-as.numeric(dt1$pH)
if (class(dt1$pH_bottom)=="factor") dt1$pH_bottom <-as.numeric(levels(dt1$pH_bottom))[as.integer(dt1$pH_bottom) ]               
if (class(dt1$pH_bottom)=="character") dt1$pH_bottom <-as.numeric(dt1$pH_bottom)
if (class(dt1$Microcystis)!="factor") dt1$Microcystis<- as.factor(dt1$Microcystis)
if (class(dt1$Chlorophyll_Sign)!="factor") dt1$Chlorophyll_Sign<- as.factor(dt1$Chlorophyll_Sign)
if (class(dt1$Chlorophyll)=="factor") dt1$Chlorophyll <-as.numeric(levels(dt1$Chlorophyll))[as.integer(dt1$Chlorophyll) ]               
if (class(dt1$Chlorophyll)=="character") dt1$Chlorophyll <-as.numeric(dt1$Chlorophyll)
if (class(dt1$Pheophytin_Sign)!="factor") dt1$Pheophytin_Sign<- as.factor(dt1$Pheophytin_Sign)
if (class(dt1$Pheophytin)=="factor") dt1$Pheophytin <-as.numeric(levels(dt1$Pheophytin))[as.integer(dt1$Pheophytin) ]               
if (class(dt1$Pheophytin)=="character") dt1$Pheophytin <-as.numeric(dt1$Pheophytin)
if (class(dt1$TotAmmonia_Sign)!="factor") dt1$TotAmmonia_Sign<- as.factor(dt1$TotAmmonia_Sign)
if (class(dt1$TotAmmonia)=="factor") dt1$TotAmmonia <-as.numeric(levels(dt1$TotAmmonia))[as.integer(dt1$TotAmmonia) ]               
if (class(dt1$TotAmmonia)=="character") dt1$TotAmmonia <-as.numeric(dt1$TotAmmonia)
if (class(dt1$DissAmmonia_Sign)!="factor") dt1$DissAmmonia_Sign<- as.factor(dt1$DissAmmonia_Sign)
if (class(dt1$DissAmmonia)=="factor") dt1$DissAmmonia <-as.numeric(levels(dt1$DissAmmonia))[as.integer(dt1$DissAmmonia) ]               
if (class(dt1$DissAmmonia)=="character") dt1$DissAmmonia <-as.numeric(dt1$DissAmmonia)
if (class(dt1$DissNitrateNitrite_Sign)!="factor") dt1$DissNitrateNitrite_Sign<- as.factor(dt1$DissNitrateNitrite_Sign)
if (class(dt1$DissNitrateNitrite)=="factor") dt1$DissNitrateNitrite <-as.numeric(levels(dt1$DissNitrateNitrite))[as.integer(dt1$DissNitrateNitrite) ]               
if (class(dt1$DissNitrateNitrite)=="character") dt1$DissNitrateNitrite <-as.numeric(dt1$DissNitrateNitrite)
if (class(dt1$TotPhos_Sign)!="factor") dt1$TotPhos_Sign<- as.factor(dt1$TotPhos_Sign)
if (class(dt1$TotPhos)=="factor") dt1$TotPhos <-as.numeric(levels(dt1$TotPhos))[as.integer(dt1$TotPhos) ]               
if (class(dt1$TotPhos)=="character") dt1$TotPhos <-as.numeric(dt1$TotPhos)
if (class(dt1$DissOrthophos_Sign)!="factor") dt1$DissOrthophos_Sign<- as.factor(dt1$DissOrthophos_Sign)
if (class(dt1$DissOrthophos)=="factor") dt1$DissOrthophos <-as.numeric(levels(dt1$DissOrthophos))[as.integer(dt1$DissOrthophos) ]               
if (class(dt1$DissOrthophos)=="character") dt1$DissOrthophos <-as.numeric(dt1$DissOrthophos)
if (class(dt1$TON_Sign)!="factor") dt1$TON_Sign<- as.factor(dt1$TON_Sign)
if (class(dt1$TON)=="factor") dt1$TON <-as.numeric(levels(dt1$TON))[as.integer(dt1$TON) ]               
if (class(dt1$TON)=="character") dt1$TON <-as.numeric(dt1$TON)
if (class(dt1$DON_Sign)!="factor") dt1$DON_Sign<- as.factor(dt1$DON_Sign)
if (class(dt1$DON)=="factor") dt1$DON <-as.numeric(levels(dt1$DON))[as.integer(dt1$DON) ]               
if (class(dt1$DON)=="character") dt1$DON <-as.numeric(dt1$DON)
if (class(dt1$TKN_Sign)!="factor") dt1$TKN_Sign<- as.factor(dt1$TKN_Sign)
if (class(dt1$TKN)=="factor") dt1$TKN <-as.numeric(levels(dt1$TKN))[as.integer(dt1$TKN) ]               
if (class(dt1$TKN)=="character") dt1$TKN <-as.numeric(dt1$TKN)
if (class(dt1$TotAlkalinity_Sign)!="factor") dt1$TotAlkalinity_Sign<- as.factor(dt1$TotAlkalinity_Sign)
if (class(dt1$TotAlkalinity)=="factor") dt1$TotAlkalinity <-as.numeric(levels(dt1$TotAlkalinity))[as.integer(dt1$TotAlkalinity) ]               
if (class(dt1$TotAlkalinity)=="character") dt1$TotAlkalinity <-as.numeric(dt1$TotAlkalinity)
if (class(dt1$DissBromide_Sign)!="factor") dt1$DissBromide_Sign<- as.factor(dt1$DissBromide_Sign)
if (class(dt1$DissBromide)=="factor") dt1$DissBromide <-as.numeric(levels(dt1$DissBromide))[as.integer(dt1$DissBromide) ]               
if (class(dt1$DissBromide)=="character") dt1$DissBromide <-as.numeric(dt1$DissBromide)
if (class(dt1$DissCalcium_Sign)!="factor") dt1$DissCalcium_Sign<- as.factor(dt1$DissCalcium_Sign)
if (class(dt1$DissCalcium)=="factor") dt1$DissCalcium <-as.numeric(levels(dt1$DissCalcium))[as.integer(dt1$DissCalcium) ]               
if (class(dt1$DissCalcium)=="character") dt1$DissCalcium <-as.numeric(dt1$DissCalcium)
if (class(dt1$TotChloride_Sign)!="factor") dt1$TotChloride_Sign<- as.factor(dt1$TotChloride_Sign)
if (class(dt1$TotChloride)=="factor") dt1$TotChloride <-as.numeric(levels(dt1$TotChloride))[as.integer(dt1$TotChloride) ]               
if (class(dt1$TotChloride)=="character") dt1$TotChloride <-as.numeric(dt1$TotChloride)
if (class(dt1$DissChloride_Sign)!="factor") dt1$DissChloride_Sign<- as.factor(dt1$DissChloride_Sign)
if (class(dt1$DissChloride)=="factor") dt1$DissChloride <-as.numeric(levels(dt1$DissChloride))[as.integer(dt1$DissChloride) ]               
if (class(dt1$DissChloride)=="character") dt1$DissChloride <-as.numeric(dt1$DissChloride)
if (class(dt1$DissSilica_Sign)!="factor") dt1$DissSilica_Sign<- as.factor(dt1$DissSilica_Sign)
if (class(dt1$DissSilica)=="factor") dt1$DissSilica <-as.numeric(levels(dt1$DissSilica))[as.integer(dt1$DissSilica) ]               
if (class(dt1$DissSilica)=="character") dt1$DissSilica <-as.numeric(dt1$DissSilica)
if (class(dt1$TOC_Sign)!="factor") dt1$TOC_Sign<- as.factor(dt1$TOC_Sign)
if (class(dt1$TOC)=="factor") dt1$TOC <-as.numeric(levels(dt1$TOC))[as.integer(dt1$TOC) ]               
if (class(dt1$TOC)=="character") dt1$TOC <-as.numeric(dt1$TOC)
if (class(dt1$DOC_Sign)!="factor") dt1$DOC_Sign<- as.factor(dt1$DOC_Sign)
if (class(dt1$DOC)=="factor") dt1$DOC <-as.numeric(levels(dt1$DOC))[as.integer(dt1$DOC) ]               
if (class(dt1$DOC)=="character") dt1$DOC <-as.numeric(dt1$DOC)
if (class(dt1$TDS_Sign)!="factor") dt1$TDS_Sign<- as.factor(dt1$TDS_Sign)
if (class(dt1$TDS)=="factor") dt1$TDS <-as.numeric(levels(dt1$TDS))[as.integer(dt1$TDS) ]               
if (class(dt1$TDS)=="character") dt1$TDS <-as.numeric(dt1$TDS)
if (class(dt1$TSS_Sign)!="factor") dt1$TSS_Sign<- as.factor(dt1$TSS_Sign)
if (class(dt1$TSS)=="factor") dt1$TSS <-as.numeric(levels(dt1$TSS))[as.integer(dt1$TSS) ]               
if (class(dt1$TSS)=="character") dt1$TSS <-as.numeric(dt1$TSS)
if (class(dt1$VSS_Sign)!="factor") dt1$VSS_Sign<- as.factor(dt1$VSS_Sign)
if (class(dt1$VSS)=="factor") dt1$VSS <-as.numeric(levels(dt1$VSS))[as.integer(dt1$VSS) ]               
if (class(dt1$VSS)=="character") dt1$VSS <-as.numeric(dt1$VSS)
if (class(dt1$Notes)!="factor") dt1$Notes<- as.factor(dt1$Notes)

# Convert Missing Values to NA for non-dates

dt1$Latitude <- ifelse((trimws(as.character(dt1$Latitude))==trimws("NA")),NA,dt1$Latitude)               
suppressWarnings(dt1$Latitude <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Latitude))==as.character(as.numeric("NA"))),NA,dt1$Latitude))
dt1$Longitude <- ifelse((trimws(as.character(dt1$Longitude))==trimws("NA")),NA,dt1$Longitude)               
suppressWarnings(dt1$Longitude <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Longitude))==as.character(as.numeric("NA"))),NA,dt1$Longitude))
dt1$Depth <- ifelse((trimws(as.character(dt1$Depth))==trimws("NA")),NA,dt1$Depth)               
suppressWarnings(dt1$Depth <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Depth))==as.character(as.numeric("NA"))),NA,dt1$Depth))
dt1$Sample_depth_surface <- ifelse((trimws(as.character(dt1$Sample_depth_surface))==trimws("NA")),NA,dt1$Sample_depth_surface)               
suppressWarnings(dt1$Sample_depth_surface <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Sample_depth_surface))==as.character(as.numeric("NA"))),NA,dt1$Sample_depth_surface))
dt1$Sample_depth_nutr_surface <- ifelse((trimws(as.character(dt1$Sample_depth_nutr_surface))==trimws("NA")),NA,dt1$Sample_depth_nutr_surface)               
suppressWarnings(dt1$Sample_depth_nutr_surface <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Sample_depth_nutr_surface))==as.character(as.numeric("NA"))),NA,dt1$Sample_depth_nutr_surface))
dt1$Sample_depth_bottom <- ifelse((trimws(as.character(dt1$Sample_depth_bottom))==trimws("NA")),NA,dt1$Sample_depth_bottom)               
suppressWarnings(dt1$Sample_depth_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Sample_depth_bottom))==as.character(as.numeric("NA"))),NA,dt1$Sample_depth_bottom))
dt1$Tide <- as.factor(ifelse((trimws(as.character(dt1$Tide))==trimws("NA")),NA,as.character(dt1$Tide)))
dt1$Temperature <- ifelse((trimws(as.character(dt1$Temperature))==trimws("NA")),NA,dt1$Temperature)               
suppressWarnings(dt1$Temperature <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Temperature))==as.character(as.numeric("NA"))),NA,dt1$Temperature))
dt1$Temperature_bottom <- ifelse((trimws(as.character(dt1$Temperature_bottom))==trimws("NA")),NA,dt1$Temperature_bottom)               
suppressWarnings(dt1$Temperature_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Temperature_bottom))==as.character(as.numeric("NA"))),NA,dt1$Temperature_bottom))
dt1$Conductivity <- ifelse((trimws(as.character(dt1$Conductivity))==trimws("NA")),NA,dt1$Conductivity)               
suppressWarnings(dt1$Conductivity <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Conductivity))==as.character(as.numeric("NA"))),NA,dt1$Conductivity))
dt1$Conductivity_bottom <- ifelse((trimws(as.character(dt1$Conductivity_bottom))==trimws("NA")),NA,dt1$Conductivity_bottom)               
suppressWarnings(dt1$Conductivity_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Conductivity_bottom))==as.character(as.numeric("NA"))),NA,dt1$Conductivity_bottom))
dt1$Salinity <- ifelse((trimws(as.character(dt1$Salinity))==trimws("NA")),NA,dt1$Salinity)               
suppressWarnings(dt1$Salinity <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Salinity))==as.character(as.numeric("NA"))),NA,dt1$Salinity))
dt1$Salinity_bottom <- ifelse((trimws(as.character(dt1$Salinity_bottom))==trimws("NA")),NA,dt1$Salinity_bottom)               
suppressWarnings(dt1$Salinity_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Salinity_bottom))==as.character(as.numeric("NA"))),NA,dt1$Salinity_bottom))
dt1$Secchi <- ifelse((trimws(as.character(dt1$Secchi))==trimws("NA")),NA,dt1$Secchi)               
suppressWarnings(dt1$Secchi <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Secchi))==as.character(as.numeric("NA"))),NA,dt1$Secchi))
dt1$Secchi_estimated <- as.factor(ifelse((trimws(as.character(dt1$Secchi_estimated))==trimws("NA")),NA,as.character(dt1$Secchi_estimated)))
dt1$TurbidityNTU <- ifelse((trimws(as.character(dt1$TurbidityNTU))==trimws("NA")),NA,dt1$TurbidityNTU)               
suppressWarnings(dt1$TurbidityNTU <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TurbidityNTU))==as.character(as.numeric("NA"))),NA,dt1$TurbidityNTU))
dt1$TurbidityNTU_bottom <- ifelse((trimws(as.character(dt1$TurbidityNTU_bottom))==trimws("NA")),NA,dt1$TurbidityNTU_bottom)               
suppressWarnings(dt1$TurbidityNTU_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TurbidityNTU_bottom))==as.character(as.numeric("NA"))),NA,dt1$TurbidityNTU_bottom))
dt1$TurbidityFNU <- ifelse((trimws(as.character(dt1$TurbidityFNU))==trimws("NA")),NA,dt1$TurbidityFNU)               
suppressWarnings(dt1$TurbidityFNU <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TurbidityFNU))==as.character(as.numeric("NA"))),NA,dt1$TurbidityFNU))
dt1$TurbidityFNU_bottom <- ifelse((trimws(as.character(dt1$TurbidityFNU_bottom))==trimws("NA")),NA,dt1$TurbidityFNU_bottom)               
suppressWarnings(dt1$TurbidityFNU_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TurbidityFNU_bottom))==as.character(as.numeric("NA"))),NA,dt1$TurbidityFNU_bottom))
dt1$DissolvedOxygen <- ifelse((trimws(as.character(dt1$DissolvedOxygen))==trimws("NA")),NA,dt1$DissolvedOxygen)               
suppressWarnings(dt1$DissolvedOxygen <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissolvedOxygen))==as.character(as.numeric("NA"))),NA,dt1$DissolvedOxygen))
dt1$DissolvedOxygen_bottom <- ifelse((trimws(as.character(dt1$DissolvedOxygen_bottom))==trimws("NA")),NA,dt1$DissolvedOxygen_bottom)               
suppressWarnings(dt1$DissolvedOxygen_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissolvedOxygen_bottom))==as.character(as.numeric("NA"))),NA,dt1$DissolvedOxygen_bottom))
dt1$DissolvedOxygenPercent <- ifelse((trimws(as.character(dt1$DissolvedOxygenPercent))==trimws("NA")),NA,dt1$DissolvedOxygenPercent)               
suppressWarnings(dt1$DissolvedOxygenPercent <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissolvedOxygenPercent))==as.character(as.numeric("NA"))),NA,dt1$DissolvedOxygenPercent))
dt1$DissolvedOxygenPercent_bottom <- ifelse((trimws(as.character(dt1$DissolvedOxygenPercent_bottom))==trimws("NA")),NA,dt1$DissolvedOxygenPercent_bottom)               
suppressWarnings(dt1$DissolvedOxygenPercent_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissolvedOxygenPercent_bottom))==as.character(as.numeric("NA"))),NA,dt1$DissolvedOxygenPercent_bottom))
dt1$pH <- ifelse((trimws(as.character(dt1$pH))==trimws("NA")),NA,dt1$pH)               
suppressWarnings(dt1$pH <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$pH))==as.character(as.numeric("NA"))),NA,dt1$pH))
dt1$pH_bottom <- ifelse((trimws(as.character(dt1$pH_bottom))==trimws("NA")),NA,dt1$pH_bottom)               
suppressWarnings(dt1$pH_bottom <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$pH_bottom))==as.character(as.numeric("NA"))),NA,dt1$pH_bottom))
dt1$Microcystis <- as.factor(ifelse((trimws(as.character(dt1$Microcystis))==trimws("NA")),NA,as.character(dt1$Microcystis)))
dt1$Chlorophyll_Sign <- as.factor(ifelse((trimws(as.character(dt1$Chlorophyll_Sign))==trimws("NA")),NA,as.character(dt1$Chlorophyll_Sign)))
dt1$Chlorophyll <- ifelse((trimws(as.character(dt1$Chlorophyll))==trimws("NA")),NA,dt1$Chlorophyll)               
suppressWarnings(dt1$Chlorophyll <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Chlorophyll))==as.character(as.numeric("NA"))),NA,dt1$Chlorophyll))
dt1$Pheophytin_Sign <- as.factor(ifelse((trimws(as.character(dt1$Pheophytin_Sign))==trimws("NA")),NA,as.character(dt1$Pheophytin_Sign)))
dt1$Pheophytin <- ifelse((trimws(as.character(dt1$Pheophytin))==trimws("NA")),NA,dt1$Pheophytin)               
suppressWarnings(dt1$Pheophytin <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$Pheophytin))==as.character(as.numeric("NA"))),NA,dt1$Pheophytin))
dt1$TotAmmonia_Sign <- as.factor(ifelse((trimws(as.character(dt1$TotAmmonia_Sign))==trimws("NA")),NA,as.character(dt1$TotAmmonia_Sign)))
dt1$TotAmmonia <- ifelse((trimws(as.character(dt1$TotAmmonia))==trimws("NA")),NA,dt1$TotAmmonia)               
suppressWarnings(dt1$TotAmmonia <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TotAmmonia))==as.character(as.numeric("NA"))),NA,dt1$TotAmmonia))
dt1$DissAmmonia_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissAmmonia_Sign))==trimws("NA")),NA,as.character(dt1$DissAmmonia_Sign)))
dt1$DissAmmonia <- ifelse((trimws(as.character(dt1$DissAmmonia))==trimws("NA")),NA,dt1$DissAmmonia)               
suppressWarnings(dt1$DissAmmonia <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissAmmonia))==as.character(as.numeric("NA"))),NA,dt1$DissAmmonia))
dt1$DissNitrateNitrite_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissNitrateNitrite_Sign))==trimws("NA")),NA,as.character(dt1$DissNitrateNitrite_Sign)))
dt1$DissNitrateNitrite <- ifelse((trimws(as.character(dt1$DissNitrateNitrite))==trimws("NA")),NA,dt1$DissNitrateNitrite)               
suppressWarnings(dt1$DissNitrateNitrite <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissNitrateNitrite))==as.character(as.numeric("NA"))),NA,dt1$DissNitrateNitrite))
dt1$TotPhos_Sign <- as.factor(ifelse((trimws(as.character(dt1$TotPhos_Sign))==trimws("NA")),NA,as.character(dt1$TotPhos_Sign)))
dt1$TotPhos <- ifelse((trimws(as.character(dt1$TotPhos))==trimws("NA")),NA,dt1$TotPhos)               
suppressWarnings(dt1$TotPhos <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TotPhos))==as.character(as.numeric("NA"))),NA,dt1$TotPhos))
dt1$DissOrthophos_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissOrthophos_Sign))==trimws("NA")),NA,as.character(dt1$DissOrthophos_Sign)))
dt1$DissOrthophos <- ifelse((trimws(as.character(dt1$DissOrthophos))==trimws("NA")),NA,dt1$DissOrthophos)               
suppressWarnings(dt1$DissOrthophos <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissOrthophos))==as.character(as.numeric("NA"))),NA,dt1$DissOrthophos))
dt1$TON_Sign <- as.factor(ifelse((trimws(as.character(dt1$TON_Sign))==trimws("NA")),NA,as.character(dt1$TON_Sign)))
dt1$TON <- ifelse((trimws(as.character(dt1$TON))==trimws("NA")),NA,dt1$TON)               
suppressWarnings(dt1$TON <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TON))==as.character(as.numeric("NA"))),NA,dt1$TON))
dt1$DON_Sign <- as.factor(ifelse((trimws(as.character(dt1$DON_Sign))==trimws("NA")),NA,as.character(dt1$DON_Sign)))
dt1$DON <- ifelse((trimws(as.character(dt1$DON))==trimws("NA")),NA,dt1$DON)               
suppressWarnings(dt1$DON <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DON))==as.character(as.numeric("NA"))),NA,dt1$DON))
dt1$TKN_Sign <- as.factor(ifelse((trimws(as.character(dt1$TKN_Sign))==trimws("NA")),NA,as.character(dt1$TKN_Sign)))
dt1$TKN <- ifelse((trimws(as.character(dt1$TKN))==trimws("NA")),NA,dt1$TKN)               
suppressWarnings(dt1$TKN <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TKN))==as.character(as.numeric("NA"))),NA,dt1$TKN))
dt1$TotAlkalinity_Sign <- as.factor(ifelse((trimws(as.character(dt1$TotAlkalinity_Sign))==trimws("NA")),NA,as.character(dt1$TotAlkalinity_Sign)))
dt1$TotAlkalinity <- ifelse((trimws(as.character(dt1$TotAlkalinity))==trimws("NA")),NA,dt1$TotAlkalinity)               
suppressWarnings(dt1$TotAlkalinity <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TotAlkalinity))==as.character(as.numeric("NA"))),NA,dt1$TotAlkalinity))
dt1$DissBromide_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissBromide_Sign))==trimws("NA")),NA,as.character(dt1$DissBromide_Sign)))
dt1$DissBromide <- ifelse((trimws(as.character(dt1$DissBromide))==trimws("NA")),NA,dt1$DissBromide)               
suppressWarnings(dt1$DissBromide <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissBromide))==as.character(as.numeric("NA"))),NA,dt1$DissBromide))
dt1$DissCalcium_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissCalcium_Sign))==trimws("NA")),NA,as.character(dt1$DissCalcium_Sign)))
dt1$DissCalcium <- ifelse((trimws(as.character(dt1$DissCalcium))==trimws("NA")),NA,dt1$DissCalcium)               
suppressWarnings(dt1$DissCalcium <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissCalcium))==as.character(as.numeric("NA"))),NA,dt1$DissCalcium))
dt1$TotChloride_Sign <- as.factor(ifelse((trimws(as.character(dt1$TotChloride_Sign))==trimws("NA")),NA,as.character(dt1$TotChloride_Sign)))
dt1$TotChloride <- ifelse((trimws(as.character(dt1$TotChloride))==trimws("NA")),NA,dt1$TotChloride)               
suppressWarnings(dt1$TotChloride <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TotChloride))==as.character(as.numeric("NA"))),NA,dt1$TotChloride))
dt1$DissChloride_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissChloride_Sign))==trimws("NA")),NA,as.character(dt1$DissChloride_Sign)))
dt1$DissChloride <- ifelse((trimws(as.character(dt1$DissChloride))==trimws("NA")),NA,dt1$DissChloride)               
suppressWarnings(dt1$DissChloride <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissChloride))==as.character(as.numeric("NA"))),NA,dt1$DissChloride))
dt1$DissSilica_Sign <- as.factor(ifelse((trimws(as.character(dt1$DissSilica_Sign))==trimws("NA")),NA,as.character(dt1$DissSilica_Sign)))
dt1$DissSilica <- ifelse((trimws(as.character(dt1$DissSilica))==trimws("NA")),NA,dt1$DissSilica)               
suppressWarnings(dt1$DissSilica <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DissSilica))==as.character(as.numeric("NA"))),NA,dt1$DissSilica))
dt1$TOC_Sign <- as.factor(ifelse((trimws(as.character(dt1$TOC_Sign))==trimws("NA")),NA,as.character(dt1$TOC_Sign)))
dt1$TOC <- ifelse((trimws(as.character(dt1$TOC))==trimws("NA")),NA,dt1$TOC)               
suppressWarnings(dt1$TOC <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TOC))==as.character(as.numeric("NA"))),NA,dt1$TOC))
dt1$DOC_Sign <- as.factor(ifelse((trimws(as.character(dt1$DOC_Sign))==trimws("NA")),NA,as.character(dt1$DOC_Sign)))
dt1$DOC <- ifelse((trimws(as.character(dt1$DOC))==trimws("NA")),NA,dt1$DOC)               
suppressWarnings(dt1$DOC <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$DOC))==as.character(as.numeric("NA"))),NA,dt1$DOC))
dt1$TDS_Sign <- as.factor(ifelse((trimws(as.character(dt1$TDS_Sign))==trimws("NA")),NA,as.character(dt1$TDS_Sign)))
dt1$TDS <- ifelse((trimws(as.character(dt1$TDS))==trimws("NA")),NA,dt1$TDS)               
suppressWarnings(dt1$TDS <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TDS))==as.character(as.numeric("NA"))),NA,dt1$TDS))
dt1$TSS_Sign <- as.factor(ifelse((trimws(as.character(dt1$TSS_Sign))==trimws("NA")),NA,as.character(dt1$TSS_Sign)))
dt1$TSS <- ifelse((trimws(as.character(dt1$TSS))==trimws("NA")),NA,dt1$TSS)               
suppressWarnings(dt1$TSS <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$TSS))==as.character(as.numeric("NA"))),NA,dt1$TSS))
dt1$VSS_Sign <- as.factor(ifelse((trimws(as.character(dt1$VSS_Sign))==trimws("NA")),NA,as.character(dt1$VSS_Sign)))
dt1$VSS <- ifelse((trimws(as.character(dt1$VSS))==trimws("NA")),NA,dt1$VSS)               
suppressWarnings(dt1$VSS <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$VSS))==as.character(as.numeric("NA"))),NA,dt1$VSS))
dt1$Notes <- as.factor(ifelse((trimws(as.character(dt1$Notes))==trimws("NA")),NA,as.character(dt1$Notes)))

# Add EDSM subregion polygon from deltamapr package
EDSM_regions <- R_EDSM_Subregions_Mahardja

# Add region to water quality data

dt_spatial <- dt1 %>% filter(!is.na(Latitude)&!is.na(Longitude)) %>%
  st_as_sf(coords=c("Longitude", "Latitude"),crs=st_crs(WW_Delta))

dt_spatial_EDSM <- st_transform(dt_spatial, crs = st_crs(EDSM_regions))

dt_spatial_EDSM<- st_join(dt_spatial_EDSM,EDSM_regions)

st_geometry(dt_spatial_EDSM) <- NULL # remove geometry, coerce to data.frame

# Create summarized data for Secchi
dt_secchi <- dt_spatial_EDSM %>% filter(Region %in% c("West","South","North") & !is.na(Secchi)) %>% mutate(Month=month(Date)) %>%
  mutate(WY=ifelse(Month>9,year(Date)+1,year(Date)),Count=1) %>% group_by(WY,Month,Region) %>%
  summarise(Secchi=mean(Secchi),SampleSize_secchi=sum(Count))

write.csv(dt_secchi,file=file.path(output_root,"raw_secchi_data.csv"),row.names = F)

# Remove data with <5 secchi data points
dt_secchi_subset <- dt_secchi %>% filter(SampleSize_secchi>=5) %>% select(-SampleSize_secchi)
#### Lowest sample size is Jan 1997, where there were only 5 sampling points for secchi in the South Delta

# Adjust data to fit fish data frames
dt_secchi_wide <- dt_secchi_subset %>% spread(Region,Secchi) %>%
  rename(Secchi_South=South,Secchi_North=North,Secchi_West=West) %>%
  # Remove Secchi_North for now because of lack of data
  select(-Secchi_North)


#############
# Pull fish salinity data and merge

fish_sal_data<-read.csv(file.path(output_root,"wild_fish_salinity_data.csv")) %>% mutate(Date_median=as.Date(Date_median_sal)) %>%
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred)) %>%
  # Add secchi data
  left_join(dt_secchi_wide)


# Export out data again
write.csv(fish_sal_data,file=file.path(output_root,"wild_fish_salinity_data_with_covariate.csv"),row.names = F)

#############
# Pull fish distance data and merge

fish_dist_data_swp<-read.csv(file.path(output_root,"wild_fish_average_distance_from_skinner.csv")) %>% 
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred)) %>%
  # Add secchi data
  left_join(dt_secchi_wide) %>%
  # Impute secchi data
  arrange(Date_median) %>%
  mutate(Secchi_South= na_interpolation(Secchi_South),Secchi_West= na_interpolation(Secchi_West))

fish_dist_data_cvp<-read.csv(file.path(output_root,"wild_fish_average_distance_from_TFCF.csv")) %>% 
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred)) %>%
  # Add secchi data
  left_join(dt_secchi_wide) %>%
  # Impute secchi data
  arrange(Date_median) %>%
  mutate(Secchi_South= na_interpolation(Secchi_South),Secchi_West= na_interpolation(Secchi_West))

fish_dist_data_cache<-read.csv(file.path(output_root,"wild_fish_average_distance_from_Cache.csv")) %>% 
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred)) %>%
  # Add secchi data
  left_join(dt_secchi_wide) %>%
  # Impute secchi data
  arrange(Date_median) %>%
  mutate(Secchi_South= na_interpolation(Secchi_South),Secchi_West= na_interpolation(Secchi_West))

# Export out data again
write.csv(fish_dist_data_swp,file=file.path(output_root,"wild_fish_SWP_distance_data_with_covariate.csv"),row.names = F)
write.csv(fish_dist_data_cvp,file=file.path(output_root,"wild_fish_CVP_distance_data_with_covariate.csv"),row.names = F)
write.csv(fish_dist_data_cache,file=file.path(output_root,"wild_fish_Cache_distance_data_with_covariate.csv"),row.names = F)
