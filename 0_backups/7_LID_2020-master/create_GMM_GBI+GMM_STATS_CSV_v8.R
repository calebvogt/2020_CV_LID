## create_GMM_GBI_CSV.R
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
# wd <- setwd("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_analysis_v5")
wd <- setwd("C:/Users/caleb/Box/CV_Shared_Analyses/7_LID_2020/RFID_analysis_v5")
# output_fp <- paste("C:/Users/Caleb Vogt/Desktop")
output_fp <- paste("C:/Users/caleb/Desktop")

# load metadata
metadata <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)

# LOAD CLEAN ALL RFID DATA FOR GMM ANALYSIS (CSV or EXCEL)
# data <- as.data.frame(fread("LID_2020_CLEAN_ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE))

# CONVERT TIME FORMAT
data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS")

# CREATE GMM DATA FILES FOR ALL TRIALS ------------------------------------
## NOTE!!! THIS CODE THROWS ERRORS FOR SOME ZONE_DAY LOCATIONS WHEN ONLY A SINGLE ANIMAL IS DETECTED AT THE ANTENNA DURING THE NIGHT. 

# LOOP THROUGH ALL TRIALS. 
trials <- unique(data$Trial)
failed_dataset <- list()
fail_flag=1
i=trials[7]
for(i in trials[1:length(trials)]){
  df <- data %>% 
    filter(Trial == i) %>% 
    select(Name, Field_Time, Time_sec, noon_to_noon_day, Days_Antenna)
  
  ids <- unique(df$Name)
  ids <- append(ids, c("dud_mouse"))
  day_zone <- unique(df$Days_Antenna)
  print(day_zone)
  aa=day_zone[8]
  for(aa in day_zone[1:length(day_zone)]){
    # SPLIT INTO DAILY CHUNKS
    daily_df <- subset(df, df$Days_Antenna == aa)
    #INSERT TRY CATCH AND WRITE WHICH DATA FRAME DIDNT WORK OUT
    tryCatch({
      
      # GENERATE GMM DATA
      gmm_data <- gmmevents(time=daily_df$Time_sec,
                            identity=daily_df$Name,
                            location=daily_df$Days_Antenna,
                            global_ids= ids, ## LEAVE AS TRUE TO GET CONSISTENT MATRICES WITH SAME COLUMN STRUCTURE
                            verbose = TRUE,
                            splitGroups = TRUE) 
      # SAVE GMM DATA
      save(gmm_data, file=paste0(output_fp,"/", i, "_GMM_RFID_",aa,".RData"))
      rm(daily_df, gmm_data) #TEST. REMOVE THE DAILY DF BEFORE STARTING WRITING OVER IT!!!!
      gc(verbose=TRUE)
    }, error=function(e){
      # ADD NAME OF FAILED DATA TO FAILED_DATASET LIST. WEIRDLY, IN THE ERROR FUNCTION ASSIGNMENT OPERATOR IS <<- INSTEAD OF <-!!!
      failed_dataset[[fail_flag]] <<- paste0(i, "_GMM_RFID_",aa,".RData")
      fail_flag <<- fail_flag+1
      
      ## ADD A SINGLE DUD MOUSE READ TO THE DATA FRAME
      dud_row <- daily_df[nrow(daily_df),]
      dud_row$Name[1] <- "dud_mouse"
      
      ## ADD 60 MINUTES TO THE LAST OBSERVED MOUSE FOR THE DUD MOUSE, WHICH SHOULD GIVE IT ITS OWN FLOCKING EVENT PRESUMABLY. 
      dud_row$Time_sec <- dud_row$Time_sec[1] + 3600
      
      ## BIND THE DUD ROW TO THE DAILY DF
      daily_dud <- rbind(daily_df, dud_row)
      
      # GENERATE GMM DATA FROM DAILY_DUD
      gmm_data <- gmmevents(time=daily_dud$Time_sec,
                            identity=daily_dud$Name,
                            location=daily_dud$Days_Antenna,
                            global_ids= ids, ## LEAVE AS TRUE TO GET CONSISTENT MATRICES WITH SAME COLUMN STRUCTURE
                            verbose = TRUE,
                            splitGroups = TRUE) 
      # SAVE GMM DATA GENERATED FROM DAILY DUD!!!! THIS WORKS!!!
      save(gmm_data, file=paste0(output_fp,"/", i, "_GMM_RFID_",aa,".RData"))
      rm(daily_df, gmm_data, daily_dud) #TEST. REMOVE THE DAILY DF BEFORE STARTING WRITING OVER IT!!!!
      gc(verbose=TRUE)
      
      
    }) 
  }
}

# WHICH DAY_ANTENNA COMBINATION DATA DOES THE GMM_CODE FAIL? PRINT THE LIST OF ERRORS. 
master_fail <- do.call(rbind,failed_dataset)
print(master_fail)

## AS OF 4/30/2021, THESE ARE THE "FAILED" DATA SETS
# [1,] "T001_GMM_RFID_1_6.RData" 
# [2,] "T002_GMM_RFID_2_3.RData" 
# [3,] "T002_GMM_RFID_2_4.RData" 
# [4,] "T002_GMM_RFID_2_1.RData" 
# [5,] "T002_GMM_RFID_3_3.RData" 
# [6,] "T002_GMM_RFID_3_4.RData" 
# [7,] "T002_GMM_RFID_6_5.RData" 
# [8,] "T003_GMM_RFID_2_7.RData" 
# [9,] "T003_GMM_RFID_2_2.RData" 
# [10,] "T003_GMM_RFID_2_3.RData" 
# [11,] "T003_GMM_RFID_3_7.RData" 
# [12,] "T003_GMM_RFID_3_3.RData" 
# [13,] "T003_GMM_RFID_6_3.RData" 
# [14,] "T003_GMM_RFID_10_3.RData"
# [15,] "T006_GMM_RFID_2_5.RData" 
# [16,] "T006_GMM_RFID_5_5.RData" 


# COMBINE TRIAL .RDATA FILES AND MERGE INTO TRIAL_GMM_GBI_DATA.CSV ----------------------------------------------

# SET TRIAL WD AND LOAD DATA. NEED TO SET NEW FOR EACH TRIAL. 
wd <- setwd(output_fp) 
trials <- unique(data$Trial)
aa = trials[2]
for(aa in trials[1:length(trials)]) {
  df <- data %>% 
    filter(Trial == aa)
  first_read <- df[1,]
  first_read$Field_Time <- as.POSIXct(first_read$Field_Time, format="%Y-%m-%d %H:%M:%OS")
  filenames <- list.files(wd, pattern = aa)
  filelist <- list()
  i = 1
  for(i in 1:length(filenames)){
    # for(i in 1:2){
    load(filenames[i])
    gbi <- gmm_data[[1]]
    metadata <- gmm_data[[2]]
    df <- cbind(gbi, metadata)
    filelist[[i]] <- df
  }
  ## MERGE
  trial_master <- do.call(rbind,filelist)
  trial_master <- trial_master[order(trial_master$Start),]
  ## ADD IN EMPTY COLUMNS FOR MATCHING
  trial_master$Field_Time_START <- as.POSIXct(NA)
  trial_master$Field_Time_STOP <- as.POSIXct(NA)
  ## TAKE FIRST TIME SEC COUNT AND SUBTRACT FROM FIELD TIME TO GET FIELD TIME OF THE START OF THE TRIAL. 
  first_time_sec <- first_read$Time_sec
  origin <-  first_read$Field_Time[1] - first_time_sec
  ## ADD TIME_SEC TO CALCULATED ORIGIN FIELD TIME TO GET TRUE FIELD TIMES OF START AND STOP OF FLOCKING EVENTS. 
  trial_master$Field_Time_START <- trial_master$Start + origin
  trial_master$Field_Time_STOP <- trial_master$End + origin
  trial_master$Field_Time_START <- as.POSIXct(trial_master$Field_Time_START,format="%Y-%m-%d %H:%M:%OS")
  trial_master$Field_Time_STOP <-  as.POSIXct(trial_master$Field_Time_STOP,format="%Y-%m-%d %H:%M:%OS")
  # Split up location and date data
  tmp <- strsplit(trial_master$Location,"_")
  tmp <- do.call("rbind",tmp)
  trial_master$Day <- tmp[,1]
  trial_master$Zone <- tmp[,2]
  trial_master$Duration <- trial_master$End - trial_master$Start
  trial_master <- trial_master %>% relocate(c(Day, Zone, Location, Start, End, Duration, Field_Time_START,  Field_Time_STOP), .before = trial_master[,1])
  write.csv(trial_master, paste0(aa, "_GMM_GBI_DATA.csv"))
  
}
#MOVE FILES TO PARENT DIRECTORY.



# SET WD, OUTPUT PATH, LOAD DATA, FORMAT TIME SERIES
wd <- setwd("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_DATA/RFID_analysis_v3")
output_fp <- paste("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_DATA/RFID_analysis_v3/0_output_plots")
metadata <- read_excel("Liddell_2020_metadata.xlsx", sheet = 1, skip = 1)

# CREATE GMM_MOUSE_STATS --------------------------------------------------

meta_short <- metadata %>% 
  select(trial, paddock, strain, sex, name, code, family_group)

## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
males <- meta_short %>% 
  filter(sex == "M") %>% 
  select(name) %>% 
  filter(!is.na(name))
male_list <- dplyr::pull(males, name)

## female names list
females <- meta_short %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  select(name) %>% 
  filter(!is.na(name))
female_list <- dplyr::pull(females, name)

## READ IN DATA
filenames <- list.files(wd, pattern = "*MOVEBOUT_GBI.csv")

## READ IN ALL FILES
myfiles = lapply(filenames, fread)

trial_stats <- list()
aa = 1
for(aa in 1:length(myfiles)){
  ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- myfiles[[aa]]
  df <- subset(df, select = -c(V1))
  
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum)
  
  colnames(df2)<-gsub("C57-M-","",colnames(df2))
  colnames(df2)<-gsub("C57-F-","",colnames(df2))
  colnames(df2)<-gsub("NYOB-M-","",colnames(df2))
  colnames(df2)<-gsub("NYOB-F-","",colnames(df2))

  ## get mouse column names starting at col 13
  col_ids <- colnames(df2[,10:ncol(df2)])
  ## get rid of dud_mouse column
  col_ids <- col_ids[!grepl('dud_mouse', col_ids)]
  
  bb = col_ids[1]
  all_mouse_list <- list()
  first_flag = 1
  for(bb in col_ids[1:length(col_ids)]) {
    df3 <- df2 %>% 
      filter((!!as.symbol(bb)) == 1) %>% 
      mutate(Name = bb) %>% 
      relocate(Name)
    
    rows <- list()
    i = col_ids[2]
    
    second_flag = 1
    for(i in col_ids[1:length(col_ids)]) {
      df4 <- df3 %>% 
        filter((!!as.symbol(i)) == 1)
      
      rows[[second_flag]] <- df4[1,]
      second_flag = second_flag + 1
    }
    
    all_mouse_list[[first_flag]] <- do.call("rbind", rows)
    first_flag <- first_flag + 1
  }
  master <- do.call("rbind", all_mouse_list)
      
    
    # select(Day, Zone, Name, Field_Time_Start, Field_Time_Stop, m_sum, f_sum, mf_sum)
    
    
    
    
    
    df4 <- df3[match(df3$Hebe == "1")]
    df3
    
    if(df3$)
    
    df4 <- df3 %>% 
      distinct(!!! syms(col_ids), .keep_all = TRUE)
    
    df4 <- df3 %>% 
      distinct(Hathor, Hebe, .keep_all = TRUE)
    
    
    z <- df3[!duplicated(df3$Hebe),]
    
   x <- apply(df3,2,function(x) (x[x>0])[1])
    
    apply(df3,1, function(x) head(x[x!=0],1))
    x<- lapply(unique(df3$ID),function(a) head(subset(df3,ID==a & Sales>0),1))
    
    
    df4 <- df3 %>% 
      filter()
    
    
    ## ADD RELEVANT METADATA INFORMATION. 
    df4 <- merge(df3, meta_short, by.x = "Name", by.y = "name")
    df5 <- df4 %>% 
      relocate(trial, paddock, Day, Zone, strain, sex, Name, code, family_group, Field_Time_Start, Field_Time_Stop, m_sum, f_sum, mf_sum)
    
    
    
    
    stats[[bb]] <- df5
  }
  trial_stats[[aa]] <- do.call("rbind", stats)
  
}

master_stats <- do.call("rbind", trial_stats)
write.csv(master_stats, "LID_2020_GMM_MOUSE_STATS.csv")
