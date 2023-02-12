## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## creates edge list of total time any two animals spent with each other per day (A-B and B-A are not repeated in the dataset, represented only once in random order)
## Testing: create edge list for every continuous interaction between mice A and B, regardless of whether mouse C enters or leaves. 
## Testing: Only break an A-B interaction bout if either A or B leaves the interaction. 
## Goal: Get time spent tolerating any other mouse. 

library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)
social_data = lapply(filenames, fread) ## READ IN ALL FILES

# clean social data for triaged mice from social interaction bouts. 
# note that this cleaning step merely deletes columns and adds 0s where appropriate. Does not adjust GBI summary information (m_sum, f_sum, mf_sum)
aa = 1
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


## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
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

trial_stats <- list()
aa = 6
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  colnames(df)<-gsub("C57-M-","",colnames(df))
  colnames(df)<-gsub("C57-F-","",colnames(df))
  colnames(df)<-gsub("NYOB-M-","",colnames(df))
  colnames(df)<-gsub("NYOB-F-","",colnames(df))
  
  mice_names <- colnames(df[,10:ncol(df)]) ## get mouse column names starting at col 10
  focal = mice_names[1]
  focal ="Amy"
  all_mouse_list <- list()
  done_mice <- c()
  flag1 = 1
  for(focal in mice_names[1:length(mice_names)-1]) { ## choose mouse focal column and compare to each successive mouse column. 
    print(focal)
    done_mice <- c(done_mice, focal) ## keep track of which mice have been used as a focal and do not use them as a comparison so as not to repeat observations. 
    df2 <- df %>% 
      filter((!!as.symbol(focal)) == 1) %>% # pull all rows where focal mouse is present.
      mutate(ID1 = focal) %>% 
      relocate(ID1)
    
    dyad_interactions <- list()
    partner_names <- mice_names[! mice_names %in% done_mice] ## remove columns of mice that have already been used as a focal. 
    partner <- partner_names[1] ## select current partner to look at. 
    partner = "Aphrodite"
    flag2 = 1
    for(partner in partner_names[1:length(partner_names)]) {
      print(paste(focal, partner))
      df3 <- df2 %>% 
        filter((!!as.symbol(partner)) == 1) %>% ## select rows where partner ==1. This line screws up field time order! fixed with arrange
        mutate(ID2 = partner) %>% 
        relocate(ID2) %>% 
        arrange(field_time_start) ## critical step
      
      ## troubleshooting tools
      # x <- df7 %>% 
      #   group_by(ID1,ID2,day) %>% 
      #   summarize(total = sum(duration_s))
      # write.table(df6, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
      ## troubleshooting tools
      
      df4 <- df3 %>% 
        mutate(start = as.POSIXct(field_time_start, origin="1970-01-01"), 
               stop = as.POSIXct(field_time_stop, origin="1970-01-01")) %>% 
        select(trial, day, zone, start, stop, ID1, ID2) %>% ## all good
        filter(!(start==stop))## remove rows with exact same start/stop time
        
      ## concatenate continuous flocking events for focal-partner pair into single row. i.e. the time those two mice ACTUALLY spent together, independent of other mouse behavior. 
      df5 <- df4 %>% 
        select(start,stop)
        
      temp <- unlist(df5) ## put all values into single vector
      myDupeVec <- unique(temp[duplicated(temp)]) ## get unique elements of vector
      noDupesList <- lapply(df5, function(i) i[!(i %in% myDupeVec)]) ## use the unique elements to remove non unique elements
       
      df6 <- do.call(data.frame, noDupesList) ## slap that baby back together. 
      
      df7 <- df4 %>% 
        merge(., df6, by="start") %>% 
        rename(stop=stop.y) %>% ## y because you want the stop from df6!!!
        select(trial, day, zone, ID1, ID2, start, stop) %>% 
        mutate(duration_s = as.numeric(stop)-as.numeric(start)) ## change to numeric so you dont get value as minutes
       
      dyad_interactions[[flag2]] <- df7
      flag2 = flag2 + 1
    }
    
    all_mouse_list[[flag1]] <- do.call("rbind", dyad_interactions)
    flag1 <- flag1 + 1
  }
  trial_stats[[aa]] <- do.call("rbind", all_mouse_list)

}
df8 <- do.call("rbind", trial_stats)
write.csv(df8, "data/ALLTRIAL_MOVEBOUT_GBI_edgelist.csv")

