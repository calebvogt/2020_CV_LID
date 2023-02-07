## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(beepr)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")

data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE))
data$field_time <- as.POSIXct(data$field_time, format="%Y-%m-%d %H:%M:%OS")
is.POSIXct(data$field_time) 
df1 <- data %>% 
  select(-c(V1)) %>% 
  filter(noon_day %in% 1:10)

## STOP!!! Ctrl+F and Change the all values to the 95% capture threshold. 

mice <- unique(df1$name)
mylist <- list()
aa = mice[1]
# aa = "Jimmy"
flag <- 1
for (aa in mice[1:length(mice)]){
  df2 <- subset(df1, name == aa)
  df2 <- unique(df2)
  df2 <- df2[order(df2$field_time),] 
  df2$bout_status <- NA
  days <- unique(df2$noon_day)
  daylist <- list()
  zz = days[9]
  for (zz in days[1:length(days)]){
    print(paste("Processing ",aa, ", mouse ",flag, " of ", length(mice), " for day ",zz, sep=''))
    df3 <- subset(df2, noon_day == zz)
    # Check #1: animal only has one read
    if(nrow(df3) == 1){
      df3$bout_status[1] <- "SINGLE_READ"  
    }
    
    # Check #2: animal only has two reads    
    if((nrow(df3) == 2) & (df3$zone[1] == df3$zone[2]) && (difftime(df3$field_time[2], df3$field_time[1], units="secs") <= 139)) {
      df3$bout_status[1] <- "START"
      df3$bout_status[2] <- "STOP"
    } 
    
    # Check #3: animal only has two reads
    if((nrow(df3) == 2) & (difftime(df3$field_time[2], df3$field_time[1], units="secs") >= 139)) {    #### CHANGE     
      df3$bout_status[1] <- "SINGLE_READ"
      df3$bout_status[2] <- "SINGLE_READ"
    }
    
    # Check #4: if none of the others apply, move to the following
    if(nrow(df3) > 2) {
      
      #check first row
      if((difftime(df3$field_time[2], df3$field_time[1], units="secs") >= 139) | (df3$zone[2] != df3$zone[1])) {    #### CHANGE     
        df3$bout_status[1] <- "SINGLE_READ"
      }
      
      #check last row
      if( (difftime(df3$field_time[nrow(df3)], df3$field_time[nrow(df3)-1], units="secs") >= 139) | (df3$zone[nrow(df3)] != df3$zone[nrow(df3)-1]) ){    #### CHANGE     
        df3$bout_status[nrow(df3)] <- "SINGLE_READ"
      }
      
      #loop 1: populate single reads based on time
      bb=2
      for(bb in 2:(nrow(df3)-1)){
        if((difftime(df3$field_time[bb+1], df3$field_time[bb], units="secs") >= 139) && (difftime(df3$field_time[bb], df3$field_time[bb-1], units="secs") >= 139)) { #### CHANGE    
          df3$bout_status[bb] <- "SINGLE_READ"
        }
      }
      
      #loop 2: populate single reads based on zone switch
      for(cc in 2:(nrow(df3)-1)){
        if((df3$zone[cc-1] != df3$zone[cc]) && (df3$zone[cc+1] != df3$zone[cc])) {     
          df3$bout_status[cc] <- "SINGLE_READ"
        }
      }
    }
    
    # PULL OUT SINGLE_READ ROWS TO NEW SINGLE_READS DF. WILL BE ADDED BACK TO DATAFRAME LATER. DURATION OF SINGLE_READS = 1 SEC. 
    SINGLE_READs <- df3 %>% 
      filter(bout_status == "SINGLE_READ")
    
    if(nrow(SINGLE_READs) >= 1) {
      SINGLE_READs$duration_s <- 1
      SINGLE_READs$field_time_STOP <- SINGLE_READs$field_time + SINGLE_READs$duration_s
      SINGLE_READs$visit_STOP <- "STOP"
    }
    
    # REMOVE SINGLE_READ ROWS 
    df4 <- df3[!grepl("SINGLE_READ", df3$bout_status),]
    
    if(nrow(df4) == 0 ){
      df7 <- SINGLE_READs
    } else{
      
      # ADD START AND STOP TO ROWS. 
      df4$bout_status[1] <- "START"
      df4$bout_status[nrow(df4)] <- "STOP"
      
      
      ### add starts and stops to all of the NA rows if they fit the criteria. 
      for (i in 1:nrow(df4)) {
        if(is.na(df4$bout_status[i])) {
          ## if the zone i is the same as zone i+1, check time difference and add start stops
          if((df4$zone[i+1] == df4$zone[i]) && (difftime(df4$field_time[i+1], df4$field_time[i], units="secs") >= 139)) { ## change
            df4$bout_status[i] <- "STOP"
            df4$bout_status[i+1] <- "START"
          }
          ## if zone i is different from zone i+1, STOP the bout in zone i and START the bout in zone i+1. Assumes dataframe is sorted by increasing time (which it is)
          if(df4$zone[i+1]!= df4$zone[i]) {
            df4$bout_status[i] <- "STOP"
            df4$bout_status[i+1] <- "START"
          }
        } 
      }
      
      # DROP ALL ROWS WITH NAS (ALL ROWS "IN-BETWEEN" START/STOP EVENTS)
      df5 <- na.omit(df4)
      
      # PULL START ROWS. 
      START_rows <- df5[df5$bout_status == 'START', ] 
      
      # PULL STOP ROWS AND RENAME COLUMNS. 
      STOP_rows <- df5[df5$bout_status == 'STOP', ]
      names(STOP_rows)[names(STOP_rows) == 'bout_status'] <- 'visit_STOP'
      names(STOP_rows)[names(STOP_rows) == 'field_time'] <- 'field_time_STOP'
      
      # REMOVE UNECESSARY COLUMNS FROM STOP ROWS
      STOP_rows <- STOP_rows %>% 
        select(field_time_STOP, visit_STOP)
      
      #IF START_ROW LENGTH IS +1 GREATER THAN STOP ROW, DELETE LAST LINE
      if(nrow(START_rows) == nrow(STOP_rows)+1) {
        START_rows <- START_rows[-nrow(START_rows),]
      }
      
      #IF STOP_ROW LENGTH IS +1 GREATER THAN START ROW, DELETE LAST LINE
      if(nrow(STOP_rows) == nrow(START_rows)+1) {
        STOP_rows <- STOP_rows[-nrow(STOP_rows),]
      }
      
      # CHECK TO SEE THAT NROWS FOR STARTS = NROWS FOR STOPS. IF NOT, PRINT MOUSE NAME AND MOVE TO NEXT ITERATION OF THE LOOP
      if(nrow(START_rows) != nrow(STOP_rows)){
        print(aa)
        next 
      }
      
      # CBIND START AND STOP ROWS INTO DUR_DUR DF
      df6 <- cbind(START_rows, STOP_rows)
      
      df6$field_time <- as.POSIXct(df6$field_time, format="%Y-%m-%d %H:%M:%OS")
      df6$field_time_STOP <- as.POSIXct(df6$field_time_STOP, format="%Y-%m-%d %H:%M:%OS")
      
      # CALCULATE DURATION OF VISITS ### check duration min or sec. 
      df6$duration_s <- df6$field_time_STOP - df6$field_time
      df6$duration_s <- as.numeric(df6$duration_s)
      # replace any 0's with 1s duration. 
      df6$duration_s[df6$duration_s == 0] <- 1
      
      #FINALLY, ADD BACK IN YOUR SINGLE_READ ROWS. ACCOUNT FOR WHEN THERE ARE NO SINGLE_READS
      if(nrow(SINGLE_READs) == 0){
        df7 <- df6
      } else{
        df7 <- bind_rows(df6, SINGLE_READs)
      }
    }
    
    daylist[[zz]] <- df7
    print(paste("Processing done!"))
    
  } 
  allmousedays <- do.call("rbind",daylist)
  mylist[[aa]] <- allmousedays
  flag<- flag+1
}

master_class <- do.call("rbind",mylist)
master_class$field_time <- as.POSIXct(master_class$field_time, format="%Y-%m-%d %H:%M:%OS")
master_class$field_time_STOP <- as.POSIXct(master_class$field_time_STOP, format="%Y-%m-%d %H:%M:%OS")
write.csv(master_class, "data/ALLTRIAL_MOVEBOUT.csv")
