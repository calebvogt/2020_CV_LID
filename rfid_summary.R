## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(plotrix)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))


# Table S1: Numbers and Descriptive Stats  -------------------------------------
# Mean RFID reads per trial

df <- rfid_data %>% 
  group_by(trial) %>%  
  tally()
mean(df$n)
plotrix::std.error(df$n)

# Mean RFID reads per mouse per night 
df <- rfid_data %>% 
  group_by(name, noon_day) %>%  
  tally()
mean(df$n)
std.error(df$n)

## 
df <- rfid_data %>% 
  group_by(trial) %>% 
  tally()

##
df <- rfid_data %>% 
  filter(trial == "T001") %>% 
  group_by(noon_day, zone) %>% 
  tally()
