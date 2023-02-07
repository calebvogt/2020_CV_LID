## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)
social_data = lapply(filenames, fread) ## READ IN ALL FILES

# clean social data for triaged mice from social interaction bouts. 
# note that this cleaning step merely deletes columns and adds 0s where appropriate. Does not adjust GBI summary information (m_sum, f_sum, mf_sum)
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  
  df2 <- df %>% 
    #T004: George only mouse to cross between trials on Day 3. triage completely (drop column)
    dplyr::select(-one_of(c( "V1", "NYOB-M-George"))) %>% 
    #T003: Anubis visually confirmed dead by seizure on day 5.  
    mutate_at(vars(one_of(c("C57-M-Anubis"))), ~ ifelse(day >= 5, 0, .)) %>%
    #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
    mutate_at(vars(one_of(c("C57-F-Rae"))), ~ ifelse(day >= 2, 0, .)) %>%
    #T004: Hare only appears day 1. Not recovered, presumed dead. 
    mutate_at(vars(one_of(c("NYOB-M-Hare"))), ~ ifelse(day >= 2, 0, .)) %>%
    #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
    mutate_at(vars(one_of(c("NYOB-F-Isis"))), ~ ifelse(day >= 3, 0, .)) %>%
    #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
    mutate_at(vars(one_of(c("C57-F-Rose"))), ~ ifelse(day >= 10, 0, .))
  
  social_data[[aa]] <- df2
}

males <- metadata %>% 
  filter(sex == "M") %>% 
  dplyr::select(name) %>% 
  filter(!is.na(name))
male_list <- dplyr::pull(males, name)

females <- metadata %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  dplyr::select(name) %>% 
  filter(!is.na(name))
female_list <- dplyr::pull(females, name)


trial_list <- list()
aa = 1
for(aa in 1:length(social_data)){
  df0 <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- df0 %>% 
    select(-contains("-F-"))

  colnames(df)<-gsub("C57-M-","",colnames(df))
  colnames(df)<-gsub("C57-F-","",colnames(df))
  colnames(df)<-gsub("NYOB-M-","",colnames(df))
  colnames(df)<-gsub("NYOB-F-","",colnames(df))
  
  zone_list <- list()
  bb=1
  for(bb in 1:8){
    #filter for >1 males
    df2 <- df %>% 
      filter(zone == bb, m_sum >0)
    
    win_list <- list()
    cc=1
    for(cc in 1:(nrow(df2)-2)) {
      if(df2$m_sum[cc]==1 & df2$m_sum[cc+1]>1 & df2$m_sum[cc+2]==1) {
        print(paste("Processing row ",cc," out of ",nrow(df2), " for zone ",bb," in trial ",aa, sep=''))
        
        trial<- df2$trial[cc+1]
        strain <- df2$strain[cc+1]
        day <- df2$day[cc+1]
        zone <- df2$zone[cc+1]
        field_time_start <- as.character(df2$field_time_start[cc+1], format="%Y-%m-%d %H:%M:%OS")
        field_time_stop <- as.character(df2$field_time_stop[cc+1], format="%Y-%m-%d %H:%M:%OS")
        duration_s <- df2$duration_s[cc+1]
        mouse1 <- names(df2[cc,10:ncol(df2)])[which(df2[cc,10:ncol(df2)] == 1, arr.ind=T)[, "col"]] # extract the column name containing 1 from row cc
        fighters <- names(df2[cc+1,10:ncol(df2)])[which(df2[cc+1,10:ncol(df2)] == 1, arr.ind=T)[, "col"]]
        mouse2 <- fighters[! fighters %in% mouse1]
        winner <- names(df2[cc+2,10:ncol(df2)])[which(df2[cc+2,10:ncol(df2)] == 1, arr.ind=T)[, "col"]] 
        loser <- fighters[! fighters %in% winner]
        winner_order <- ifelse(winner == mouse1, "mouse1", "mouse2")
        
        win_list[[cc]] <- cbind(trial,
                            strain,
                            day,
                            zone,
                            field_time_start, 
                            field_time_stop, 
                            duration_s, 
                            mouse1,
                            mouse2, 
                            winner_order, 
                            winner,
                            loser)
      }
    }
    zone_list[[bb]] <- do.call("rbind", win_list)
  }
  trial_list[[aa]] <- do.call("rbind", zone_list)
}
terr_disputes <- do.call("rbind", trial_list)
terr_disputes <- as.data.frame(terr_disputes)
terr_disputes$field_time_start <- as.POSIXct(terr_disputes$field_time_start, format="%Y-%m-%d %H:%M:%OS")
terr_disputes$field_time_stop <- as.POSIXct(terr_disputes$field_time_stop, format="%Y-%m-%d %H:%M:%OS")
displace <- terr_disputes %>% 
  mutate_at(vars(contains("field_time")), as.POSIXct) %>% 
  mutate_at(vars(contains(c("day", "zone", "duration_s"))), as.numeric) %>% 
  arrange(trial,zone,field_time_start) %>% 
  mutate(trial_zone = paste(trial, zone, sep="_"))

## Load terr ownership
terr <- as.data.frame(fread("data/ALLTRIAL_RFID_zone_ownership.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
terr1 <- terr %>% 
  filter(rank_order==1) %>% 
  # filter(noon_day == 10) %>% ## deprecated: get day 10 rank 1s as opposed to daily rank 1s. 
  mutate(trial_zone = paste(trial, zone, sep="_")) %>% 
  mutate(trial_zone_day = paste(trial, zone, noon_day, sep="_"))

## Load PAS data
pas <- as.data.frame(fread("data/priority_access_scores_males.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
pas <- pas %>% 
  group_by(name) %>% 
  filter(noon_day == max(noon_day)) ## get max day observed. For Anubis, this is day 4 (he then dies)

## Merge winner PAS attributes
data1 <- merge(displace, pas[,c("strain","name", "csum_daily_capture_penalty")], by.x = "winner", by.y = "name", all.x=T)
data2 <- data1 %>% 
  rename(winner_PAS = csum_daily_capture_penalty) %>% 
  mutate(winner_PAS_status = ifelse(winner_PAS>0, "high", "low")) %>% 
  mutate(trial_zone_day = paste(trial, zone, day, sep="_"))

## Merge loser PAS attributes
data3 <- merge(data2, pas[,c("strain","name", "csum_daily_capture_penalty")], by.x = "loser", by.y = "name", all.x=T)
data4 <- data3 %>% 
  rename(loser_PAS = csum_daily_capture_penalty) %>% 
  mutate(loser_PAS_status = ifelse(loser_PAS>0, "high", "low"))

## deprecated: Merge Day 10 zone owner attributes 
# df <- merge(data4, terr1[,c("trial_zone", "name", "perc_zone_reads")], by.x = "trial_zone", by.y = "trial_zone", all.x=T)

# Merge evolving zone owner attributes over days
df <- merge(data4, terr1[,c("trial_zone_day", "name", "perc_zone_reads")], by.x = "trial_zone_day", by.y = "trial_zone_day", all.x=T)

## fix any missing strain values
df$strain <- ifelse(df$trial == "T001", "C57", 
                    ifelse(df$trial == "T002", "C57", 
                           ifelse(df$trial == "T003", "C57",
                                  ifelse(df$trial == "T004", "Outbred",
                                         ifelse(df$trial == "T005", "Outbred",
                                                ifelse(df$trial == "T006", "C57",
                                                       ifelse(df$trial == "T007", "Outbred", NA)))))))
## add home/away info
df1 <- df %>% 
  select(!c(trial_zone)) %>% 
  rename(zone_owner = name, zone_owner_perc_reads = perc_zone_reads) %>% 
  mutate(pre_dispute_type = ifelse(mouse1 == zone_owner, "R", "I"),
         dispute_type = ifelse(mouse1 == zone_owner | mouse2 == zone_owner, "RI", "II"), 
         post_dispute_type = ifelse(winner == zone_owner, "R", "I"),
         interaction_type = paste(pre_dispute_type, dispute_type, post_dispute_type, sep="_"), 
         winner_loc = ifelse(winner == zone_owner, "home", "away"),
         loser_loc = ifelse(loser == zone_owner, "home", "away")) %>% 
  select(trial, strain, day, field_time_start, field_time_stop, duration_s, 
         zone, zone_owner, zone_owner_perc_reads, mouse1, mouse2, 
         winner, winner_loc, winner_PAS, winner_PAS_status, winner_order,
         loser, loser_loc, loser_PAS, loser_PAS_status, 
         pre_dispute_type, dispute_type, post_dispute_type, interaction_type)


displace_final <- df1 %>% 
  merge(., terr2[, c("name","day", "zones_owned")], by.x = c("winner", "day"), by.y = c("name", "day"), all.x = T) %>% 
  select(zones_owned, everything()) %>% 
  rename(winner_zones_owned = zones_owned) %>% 
  merge(., terr2[, c("name","day", "zones_owned")], by.x = c("loser", "day"), by.y = c("name", "day"), all.x = T) %>% 
  select(zones_owned, everything()) %>% 
  rename(loser_zones_owned = zones_owned) %>% 
  select(trial, strain, day, field_time_start, field_time_stop, duration_s, 
         zone, zone_owner, zone_owner_perc_reads, mouse1, mouse2, winner_order,
         winner, winner_loc, winner_PAS, winner_PAS_status, winner_zones_owned,
         loser, loser_loc, loser_PAS, loser_PAS_status, loser_zones_owned,
         pre_dispute_type, dispute_type, post_dispute_type, interaction_type)

## check for NAs
colSums(is.na(displace_final))

write.csv(displace_final, "data/ALLTRIAL_MOVEBOUT_GBI_displace.csv")
