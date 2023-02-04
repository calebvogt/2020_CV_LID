## create_alltrial_data_csv.R
## Caleb C. Vogt, Cornell University

# LOAD PACKAGES
library(tidyverse)
library(readxl)
library(plyr)
library(gganimate)
library(data.table)
library(readr)
library(reshape)
library(lubridate)
library(scales)
library(plot.matrix)
library(asnipe)
library(igraph)
library(transformr)

# SET WD AND LOAD DATA
wd <- setwd("C:/Users/caleb/Desktop/Liddell_2020_RFID_full_10day_data")
filenames <- list.files(wd, pattern = "*RFID_full_data*")


all_trial_list <- list()
for(i in 1:7){
  df <- fread(filenames[i])
  
  # CHANGE FIELD.TIME COL TO POSIXct
  df$Field.Time <- as.POSIXct(df$Field.Time, format="%Y-%m-%d %H:%M:%OS")
  
  # EXTRACT DESIRED COLUMNS (LEGACY CODE, NOT PARTICULARLY EFFICIENT)
  Trial <- df$trial
  Paddock <- df$paddock
  Date <-df$scan.date
  Time <-df$scan.time
  Field_Time <- df$Field.Time
  Antenna <- df$antenna.id
  Strain <- df$strain
  Sex <- df$sex
  Full_ID <- df$full_ids
  
  # CREATE RFID DATA FRAME
  rfid <- cbind.data.frame(Trial, 
                           Paddock, 
                           Date, 
                           Time, 
                           Field_Time, 
                           Antenna,
                           Strain,
                           Sex,
                           Full_ID)
  
  class(rfid$Field_Time)
  rfid$Time_sec <- as.numeric(difftime(Field_Time,min(Field_Time),units="secs")+1)
  
  origin <- as.POSIXct(paste(rfid$Date[1], "12:00:00", sep =" "), format="%m/%d/%Y %H:%M:%OS")
  rfid$noon_to_noon_day <- ceiling(difftime(rfid$Field_Time, origin,  units = "days"))
  rfid$Days_Antenna <- paste(rfid$noon_to_noon_day,rfid$Antenna,sep="_")
  
  # SUBSET ONLY FOR 10 DAY PERIOD, NOON TO NOON. 
  rfid <- subset(rfid, noon_to_noon_day <= 10) # single = sign, else deletes repeated time rows. 
  unique(rfid$noon_to_noon_day)
  
  all_trial_list[[i]] <- rfid
  
}

all_trial_data = do.call(rbind, all_trial_list)

write.csv(all_trial_data, "Liddell_2020_ALLTRIAL_DATA.csv")


