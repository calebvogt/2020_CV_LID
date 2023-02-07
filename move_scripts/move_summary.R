## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(data.table)
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(plotrix)

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, 
                                 fill = TRUE, header = TRUE, check.names = TRUE))

df <- move_data %>% 
  filter(trial == "T001") %>% 
  group_by(noon_day, zone) %>% 
  tally()

## Mean estimated hours in the tubs per trial
df <- move_data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(trial) %>% 
  tally(sum(duration_hours)) 
sum(df$n)
summary(df$n)
mean(df$n)
std.error(df$n)

## Mean estimated hours in the tub per mouse per night
df <- move_data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(name, noon_day) %>% 
  tally(sum(duration_hours)) 
sum(df$n)
summary(df$n)
mean(df$n)
std.error(df$n)
t.test(df$n,conf.level = 0.95) ## get 95% confidence interval


## Mean resource zones visited per mouse per night 
df <- move_data %>% # 
  group_by(name, noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(name, noon_day) %>%
  tally() # get number of distinct zone visits
mean(df$n)
std.error(df$n)


