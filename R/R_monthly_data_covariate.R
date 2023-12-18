# Purpose: Monthly data covariate
# Author: Brian Mahardja
# Date: 2023-12-11

setwd("~/GitHub/salmon_negbinmodel")

library(tidyverse)
library(rvest)
library(lubridate)
library(rgdal)
library(janitor)
library(sharpshootR)


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

FPT_daily <- CDECquery(id='FPT', sensor=20, interval='D', start=as.Date("1976-10-01"), end=as.Date("2022-12-31"))

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
#############
# Pull fish salinity data and merge

fish_sal_data<-read.csv(file.path(output_root,"wild_fish_salinity_data.csv")) %>% mutate(Date_median=as.Date(Date_median_sal)) %>%
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred))


# Export out data again
write.csv(fish_sal_data,file=file.path(output_root,"wild_fish_salinity_data_with_covariate.csv"),row.names = F)

#############
# Pull fish distance data and merge

fish_dist_data_swp<-read.csv(file.path(output_root,"wild_fish_average_distance_from_skinner.csv")) %>% 
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred))

fish_dist_data_cvp<-read.csv(file.path(output_root,"wild_fish_average_distance_from_TFCF.csv")) %>% 
  left_join(FPT_data_firstflush) %>%
  mutate(FirstFlushOccurred=ifelse(Date_median>=FirstFlush_firstDate,"Yes","No")) %>%
  left_join(data_month_combined) %>%
  # There was never first flush in WY 1977, 2021, and 2022
  mutate(FirstFlushOccurred=ifelse(is.na(FirstFlushOccurred),"No",FirstFlushOccurred))

# Export out data again
write.csv(fish_dist_data_swp,file=file.path(output_root,"wild_fish_SWP_distance_data_with_covariate.csv"),row.names = F)
write.csv(fish_dist_data_cvp,file=file.path(output_root,"wild_fish_CVP_distance_data_with_covariate.csv"),row.names = F)
