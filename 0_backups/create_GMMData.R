## create_GMMData.R
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
wd <- setwd("C:/Users/Caleb Vogt/Desktop/Liddell_2020_RFID_full_10day_data")
metadata <- read_excel("Liddell.2020.xlsx", sheet = 1, skip = 1)
filenames <- list.files(wd, pattern = "*RFID_full_data*")

# MEMORY CHECK AND CLEAR
memory.size()  # CHECK
gc(verbose = TRUE)
memory.size()  # CHECK
memory.limit()
memory.limit(8000000)

# SET FILE 
for(i in 6:6){
  df <- fread(filenames[i])
  
  # CHANGE FIELD.TIME COL TO POSIXct
  df$Field.Time <- as.POSIXct(df$Field.Time, format="%Y-%m-%d %H:%M:%OS")
  
  # PULL OUT CURRENT TRIAL 
  trial <- substr(filenames[i],1,4)
  id <- unique(df$full_ids)
  Tag <- df$full_ids
  Date <-df$scan.date
  Time <-df$scan.time
  Field_Time <- df$Field.Time
  Antenna <- df$antenna.id
  memory.size()  # CHECK
  
  # CREATE RFID DATA FRAME
  rfid <- cbind.data.frame(Tag, Date, Time, Field_Time, Antenna)
  class(rfid$Field_Time)
  rfid$Time_sec <- as.numeric(difftime(Field_Time,min(Field_Time),units="secs")+1)
  
  origin <- as.POSIXct(paste(rfid$Date[1], "12:00:00", sep =" "), format="%m/%d/%Y %H:%M:%OS")
  rfid$noon_to_noon_day <- ceiling(difftime(rfid$Field_Time, origin,  units = "days"))
  rfid$Days_Antenna <- paste(rfid$noon_to_noon_day,rfid$Antenna,sep="_")
  global_ids <- unique(rfid$Tag)
  days <- unique(rfid$noon_to_noon_day)
  memory.size()  # CHECK
  
  #REMOVE UNNECESSARY COLUMNs FROM RFID
  rfid <- rfid[, c("Tag", "Field_Time", "Time_sec", "noon_to_noon_day", "Days_Antenna")]
  object.size(rfid)
  memory.size()  # CHECK
  
  
  #REMOVE UNNECESSARY OBJECTS FROM ENVIRONMENTAL MEMORY IF NECESSARY!
  rm(df, Antenna, Date, Field_Time, Tag, Time)
  memory.size()  # CHECK
  gc(verbose = TRUE)
  memory.size()  # CHECK
  
  
  ## I FUCKING DID IT BIATCH>>
  # for(aa in 3:length(days)){
  for(aa in 3:3){
    # SPLIT INTO DAILY CHUNKS
    daily_df <- subset(rfid, rfid$noon_to_noon_day == aa)
    memory.size()  # CHECK
    
    # GENERATE GMM DATA
    gmm_data <- gmmevents(time=daily_df$Time_sec,
                          identity=daily_df$Tag,
                          location=daily_df$Days_Antenna,
                          global_ids= global_ids,
                          verbose = TRUE,
                          splitGroups = TRUE) # Try this with FALSE
    
    # SAVE GMM DATA
    save(gmm_data, file=paste0(trial, "_GMM_RFID_Day_",aa,".RData"))
    rm(daily_df, gmm_data) #TEST. REMOVE THE DAILY DF BEFORE STARTING WRITING OVER IT!!!!
    gc(verbose=TRUE)
  }
}

