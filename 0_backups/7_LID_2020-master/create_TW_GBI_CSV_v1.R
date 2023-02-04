## create_TW_GBI_CSV.R
## Caleb C. Vogt, Cornell University

# LOAD PACKAGES
library(tidyverse)
library(readxl)
library(plyr)
library(dplyr)
library(data.table)
library(readr)
library(reshape)
library(lubridate)
library(scales)
library(plot.matrix)
library(asnipe)
library(transformr)

# SET WD, OUTPUT PATH, LOAD DATA, FORMAT TIME SERIES
wd <- setwd("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_DATA/RFID_analysis_v3")

# LOAD ALL RFID DATA 
data <- as.data.frame(fread("LID_2020_ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE))

# CONVERT TIME FORMAT
data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS")

## filter data down to 10 nights between 6pm and 6 am. 
data2 <- data %>%
  filter(noon_to_noon_day == c(1:10),
         Time <= "06:00:00" | Time >= "18:00:00")

unique(data2$noon_to_noon_day)

# CREATE GBI CSVS FOR ALL TRIALS ------------------------------------

# LOOP THROUGH ALL TRIALS. THIS DOESNT ADDRESS THE ISSUE WITH ANIMALS CROSSING OVER. 
# IF A T002 MALE CROSSED INTO T003, THIS WOULDNT BE CAPTURED HERE

trials <- unique(data2$Trial)
failed_dataset <- list()
fail_flag=1
i=trials[1]
for(i in trials[1:length(trials)]){
  df <- data2 %>% 
    filter(Trial == i) 
  
  first_read <- df[1,]
  first_read$Field_Time <- as.POSIXct(first_read$Field_Time, format="%Y-%m-%d %H:%M:%OS")
  
  df1 <- df %>% 
    select(noon_to_noon_day, Time_sec, Name, Days_Antenna)
  
  
  ## RENAME COLUMNS
  colnames(df1) = c("Date",
                    "Time", 
                    "ID",
                    "Location")
  
  df1$Location <- as.character(df1$Location)
  df1 <- unique(df1)
  ids <- sort(unique(df1$ID))
  
  # create GBI for entire trial across all locations and days
  group_by_individual <- get_associations_points_tw(df1, time_window = 60, which_days = NULL, which_locations = NULL) 
  
  gbi <- group_by_individual[[1]]
  times <- group_by_individual[[2]]
  days <- group_by_individual[[3]]
  locations <- group_by_individual[[4]]
  
  ## unfortunately does not provide a start/stop time. 
  gbi_combo <- as.data.frame(cbind(days, times, locations, gbi))
  gbi_combo$times <- as.numeric(gbi_combo$times)
  
  ## TAKE FIRST TIME SEC COUNT AND SUBTRACT FROM FIELD TIME TO GET FIELD TIME OF THE START OF THE TRIAL. 
  gbi_combo$Field_Time <- as.POSIXct(NA)
  first_time_sec <- first_read$Time_sec
  origin <-  first_read$Field_Time[1] - first_time_sec
  gbi_combo$Field_Time <- gbi_combo$times + origin
  gbi_combo$Field_Time <- as.POSIXct(gbi_combo$Field_Time, format = "%Y-%m-%d %H:%M:%OS")
  
  # Split up location and date data
  tmp <- strsplit(gbi_combo$locations,"_")
  tmp <- do.call("rbind",tmp)
  gbi_combo$Zone <- tmp[,2]
  
  
  gbi_combo <- gbi_combo %>% relocate(c(Field_Time, Zone), .before = locations)
  
  write.csv(gbi_combo, paste0(i, "_TW_GBI_DATA.csv"))

}
