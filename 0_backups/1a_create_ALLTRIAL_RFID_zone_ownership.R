## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(lubridate)
library(plotrix)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

df0 <- rfid_data %>% 
  filter(!(read_tag == "982.126057708741" & scan.date >= "07/24/2020")) %>%  #T005: Jeeves tag falls out of body on top of antenna)
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

# options(scipen = 3)
df1 <- df0 %>% 
  filter(sex=="M") %>% 
  group_by(trial, zone, noon_day) %>% 
  mutate(total_zone_reads=n()) %>% 
  group_by(trial, strain, name, zone, noon_day, total_zone_reads) %>% 
  tally() %>% 
  rename(mus_zone_reads = n) %>% 
  mutate(mus_perc_zone_reads = (mus_zone_reads/total_zone_reads)*100) %>% 
  group_by(trial, zone, noon_day) %>% 
  mutate(rank_order = rank(-mus_perc_zone_reads)) %>% ## create rank order
  mutate(trial_zone_day = paste(trial,zone, noon_day, sep = "_")) %>% 
  arrange(trial_zone_day, rank_order)

write.csv(df1, "data/ALLTRIAL_RFID_zone_ownership.csv")