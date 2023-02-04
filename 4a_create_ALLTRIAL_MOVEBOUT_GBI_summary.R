## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## ALLTRIAL_MOVEBOUT_GBI_summary contains the flocking events for all individual animals in a large list
## thus, if animal A and animal B are in a flock together, there will be two rows in the summary data, one for animal A and animal B that will
## identical except for the column designating the focal animal. 

library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)

##

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

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

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
trial_list = lapply(filenames, fread) ## READ IN ALL FILES
trial_stats <- list()
aa = 1
for(aa in 1:length(trial_list)){
  df <- trial_list[[aa]] ## PULL OUT EACH trial'S DATAFRAME
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum) %>% 
    select(!("V1"))
  
  ## get mouse column names starting at col 10 (Check this and confirm)
  col_ids <- colnames(df2[,10:ncol(df2)]) ## CHANGE
  col_ids
  stats <- list()
  bb = col_ids[1]
  for(bb in col_ids[1:length(col_ids)]) {
    df3 <- df2 %>% 
      # select(Day, zone, Start, End, Duration, field_time_START,field_time_STOP, (bb)) %>% 
      filter((!!as.symbol(bb)) == 1) %>% 
      mutate(name = bb) %>% 
      select(day, field_time_start, field_time_stop, zone, duration_s, name,  m_sum, f_sum, mf_sum)
    
    df4 <- merge(df3, meta_short, by.x = "name", by.y = "mouse") ## ADD RELEVANT METADATA INFORMATION. 
    df5 <- df4 %>% 
      relocate(trial, day, zone, strain, sex, name, code, name, 
               field_time_start, field_time_stop, duration_s, m_sum, f_sum, mf_sum) %>% 
      select(!(name.y))
    stats[[bb]] <- df5
  }
  trial_stats[[aa]] <- do.call("rbind", stats)
}
master_stats <- do.call("rbind", trial_stats)
write.csv(master_stats, "data/ALLTRIAL_MOVEBOUT_GBI_summary.csv")


