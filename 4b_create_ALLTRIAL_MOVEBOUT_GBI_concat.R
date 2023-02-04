## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## ALLTRIAL_MOVEBOUT_GBI_concat contains the flocking events for each trial and the sums of male and female participants. 
## if animal A and animal B are in a flock together, there will be a single rows in the concat data. Thus this is different than the 
## ALLTRIAL_MOVEBOUT_GBI_summary data, and is used for different analyses

library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
trial_list = lapply(filenames, fread) ## READ IN ALL FILES
trial_stats <- list()
aa = 1
for(aa in 1:length(trial_list)){
  df <- trial_list[[aa]] ## PULL OUT EACH trial'S DATAFRAME
  trial_stats[[aa]] <- df %>% 
    select(trial, day, zone, zone, field_time_start, field_time_stop, duration_s, m_sum, f_sum, mf_sum)
}
master_stats <- do.call("rbind", trial_stats)
write.csv(master_stats, "data/ALLTRIAL_MOVEBOUT_GBI_concat.csv")


