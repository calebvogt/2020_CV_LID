## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(lme4)
library(lmerTest)

metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

## Load terr ownership
terr <- as.data.frame(fread("data/ALLTRIAL_RFID_zone_ownership.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# Explore terr ownership dataset ------------------------------------------
terr1 <- terr %>% 
  filter(rank_order==1) %>% 
  mutate(trial_zone = paste(trial, zone, sep="_")) %>% 
  mutate(trial_zone_day = paste(trial, zone, noon_day, sep="_"))

## Determine which males owned one or two territories by day, print out which ones they own
terr2 <- terr %>%
  filter(rank_order==1) %>% 
  group_by(trial, strain, noon_day, name) %>%
  summarise(zones_owned = toString(unique(zone)), .groups = 'drop') %>% 
  rename(day = noon_day)


## Determine which males owned one or two territories by day, get counts of zones owned
terr3 <- terr %>%
  filter(rank_order==1) %>% ## , noon_day == 10
  group_by(trial, strain, noon_day, name) %>%
  count() %>% 
  rename(zones_owned = n, day = noon_day)





