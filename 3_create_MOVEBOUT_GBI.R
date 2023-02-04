## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## MOVEBOUT_GBI creates flocking events based on overlapping movebouts from every individual mouse within a zone. 
## thus,  

library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")

data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE))
data$field_time <- as.POSIXct(data$field_time, format="%Y-%m-%d %H:%M:%OS")
is.POSIXct(data$field_time) 
data$field_time_STOP <- as.POSIXct(data$field_time_STOP, format="%Y-%m-%d %H:%M:%OS")
data$field_time <- with_tz(data$field_time, "UTC")
data$field_time_STOP <- with_tz(data$field_time_STOP, "UTC")
trials <- unique(data$trial)
df1 <- data %>% 
  select(-c(V1)) %>% 
  filter(noon_day %in% 1:10)

metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

## Currently, this GBI matrix puts 1s in mouse columns if the mouse is participating in the flocking event based off its overlapping estimated movebout. 
trial_list <- list()
aa = trials[1] 
for(aa in trials[3:length(trials)]){
  df <- data %>%
    filter(trial == aa) %>% 
    mutate(mouse = paste(strain,sex,name, sep = "-"))
  
  mice <- unique(df$mouse)  
  zone_list <- list()
  zones <- unique(df$zone)
  bb = zones[1]
  for(bb in zones[1:length(zones)]){
    df1 <- df %>% 
      filter(zone == bb)
    df1 <- df1[order(df1$field_time, na.last=FALSE), ]
    df1$field_time <- as.POSIXct(df1$field_time, format="%Y-%m-%d %H:%M:%OS")
    df1$field_time_STOP <- as.POSIXct(df1$field_time_STOP, format="%Y-%m-%d %H:%M:%OS")
    START <- as.data.frame(df1$field_time)
    colnames(START) <- "x"
    STOP <- as.data.frame(df1$field_time_STOP)
    colnames(STOP) <- "x"
    start_stop <- as.data.frame(rbind(START, STOP)) # create long row of starts and stops and order them
    colnames(start_stop) <- "x"
    start_stop <- as.data.frame(start_stop[order(start_stop$x, na.last=FALSE),])
    colnames(start_stop) <- "x"
    
    list <- list()
    cc=1
    for(cc in 1:(nrow(start_stop)-1)) {
      field_time_start <- as.character(start_stop$x[cc], format="%Y-%m-%d %H:%M:%OS")
      field_time_stop <- as.character(start_stop$x[cc+1], format="%Y-%m-%d %H:%M:%OS")
      list[[cc]] <- cbind(field_time_start, field_time_stop)
    }
    gbi <- as.data.frame(do.call("rbind", list)) # Create group by individual matrix
    gbi$field_time_start <- as.POSIXct(gbi$field_time_start, format="%Y-%m-%d %H:%M:%OS") #convert to posixct
    gbi$field_time_stop <- as.POSIXct(gbi$field_time_stop, format="%Y-%m-%d %H:%M:%OS")
    gbi$field_time_start <- force_tz(gbi$field_time_start, tzone = "UTC") #force convert to UTC without changing time
    gbi$field_time_stop <- force_tz(gbi$field_time_stop, tzone = "UTC")
    
    gbi$duration_s <- gbi$field_time_stop - gbi$field_time_start #create duration
    gbi_cols <- c("trial", "day", "zone")
    gbi[gbi_cols] <- NA
    gbi$center_time <- gbi$field_time_start + (gbi$duration_s/2) # create center_time which will be compared later to determine participation in grouping event
    gbi[mice] <- NA     ## add columns for all trial mice
    
    dd = 1
    for(dd in 1:nrow(gbi)) {
      center <- as.POSIXct(gbi$center_time[dd]) # find the center point of the visitation event. 
      print(paste("Processing row ",dd," out of ",nrow(gbi), " for zone ",bb," in trial ",aa, sep=''))
      ff =1
      for(ff in 1:nrow(df1)){
        int <- interval(df1$field_time[ff], df1$field_time_STOP[ff])
        int
        if(center %within% int) {
          cool_mouse <- df1$mouse[ff]
          # add a 1 to the mouse column present in the interaction
          gbi[[cool_mouse]][[dd]] <- 1 ## critical step. 
          gbi$day[dd] <- df1$noon_day[ff]
        }
      }
    }
    
    gbi2 <- gbi[!is.na(gbi$day),] # now, remove grouping bouts with no detected animals. 
    gbi2$trial <- paste(aa) # add other details
    gbi2$zone <- paste(bb)
    gbi2[is.na(gbi2)] <- 0 #replace na's with 0s
    zone_list[[bb]] <- subset(gbi2, select = -c(center_time)) #drop center time
  }
  
  ## new code
  meta_short <- metadata %>% 
    mutate(mouse = paste(strain, sex, name, sep = "-")) %>% 
    select(trial, strain, sex, name, code, mouse)
  
  ## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
  males <- meta_short %>% 
    filter(sex == "M") %>% 
    select(mouse) %>% 
    filter(!is.na(mouse))
  male_list <- dplyr::pull(males, mouse)
  
  ## female names list
  females <- meta_short %>% 
    filter(sex == "F", na.rm = TRUE) %>% 
    select(mouse) %>% 
    filter(!is.na(mouse))
  female_list <- dplyr::pull(females, mouse)
  
  df <- do.call("rbind", zone_list) 
  current_trial <- unique(df$trial)
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum)
  
  df3 <- df2 %>% 
    relocate(trial, day, zone, field_time_start, field_time_stop, duration_s, m_sum, f_sum, mf_sum)
  
  write.csv(df3, paste0("data/", current_trial, "_MOVEBOUT_GBI.csv"))
  
}
