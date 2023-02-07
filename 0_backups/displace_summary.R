## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(lme4)
library(lmerTest)

metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

## Load displace data
displace <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT_GBI_displace.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

## find 95%CI, drop outliers overall. 

df <- displace %>% 
  filter(day %in% 4:10)
hist(df$duration_s)
mean(displace$duration_s)


# Explore displace dataset ------------------------------------------------
## check for NAs
colSums(is.na(displace))

## count of dispute events involving territory holders. 
df <- displace %>% 
  filter(dispute_type=="RI") %>% 
  count()

## count of total dispute events and mean/median/95% CI per territory owner
df <- displace %>% 
  filter(dispute_type=="RI") %>% 
  group_by(dispute_type, zone_owner) %>% 
  count()
mean(df$n)
median(df$n)
t.test(df$n,conf.level = 0.95) ## get 95% confidence interval

## Counts of dispute types
df <- displace %>% 
    group_by(strain, dispute_type) %>% 
    count()
  
## Counts of interaction types by strain
df <- displace %>% 
  group_by(strain, interaction_type) %>% 
  count()

## Counts of male wins at home and away 
df <- displace %>% 
  group_by(strain, winner, winner_loc) %>% 
  count()

## number of displacement events involving territory holders
df <- displace %>% 
  filter(!is.na(winner_zones_owned)) %>% ## only get away wins from mice who own zones
  filter(winner_loc == "away") %>% 
  mutate(count = n()) %>% 
  group_by(day, count) %>% 
  tally() %>% 
  mutate(daily_perc = n/count*100)


## number of away wins by territory holders
away_wins <- displace %>% 
  filter(!is.na(winner_zones_owned)) %>% ## only get away wins from mice who own zones
  filter(winner_loc == "away") %>% 
  mutate(count = n()) %>% 
  group_by(day, count) %>% 
  tally() %>% 
  mutate(daily_perc = n/count*100)


