## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
meta <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

# PAS_males ---------------------------------------------------------------
df <- move_data %>% 
  filter(sex  == "M") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, noon_day, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(trial_day_zone = paste0(trial, "_",noon_day, "_", zone)) %>% 
  dplyr::rename(total_duration_min = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  dplyr::select(trial_day_zone, total_duration_min)

## individual level zone duration usage per day. 
df1 <- move_data %>% 
  filter(sex  == "M") %>%
  mutate(duration_min = duration_s / 60) %>% 
  mutate(trial_day_zone = paste0(trial, "_", noon_day, "_", zone)) %>% 
  group_by(name, trial_day_zone, noon_day, zone) %>% 
  tally(sum(duration_min)) %>%
  dplyr::rename(mus_duration_min = n)#specify dplyr due to conflict

df2 <- merge(df1, df, by = "trial_day_zone", all = TRUE) # bring in rest of males that did not win any days. 

df3 <- merge(df2, meta, by = "name", all = FALSE) # bring in metadata

df4 <- df3 %>%
  dplyr::select(trial, name, noon_day, zone, mus_duration_min, total_duration_min)

#summing daily adjusted scores and taking single number to avoid guessing of slice_max when scores vacillate in the negative range. 
df5 <- df4 %>% 
  mutate(mus_percent_capture = (mus_duration_min / total_duration_min)) %>% 
  group_by(trial, name, noon_day) %>% 
  mutate(penalty = if(any(mus_percent_capture > 0.5)) 0 else -1) %>% # PENALTY #1: If on any day you dont capture greater than 50% for any zone, take off 1 point 
  group_by(trial, name) %>%
  complete(noon_day = 1:10, fill = list(penalty = -1, mus_percent_capture = 0)) %>% # PENALTY #2: Not observed at all penalty. Doesnt effect anyone after filtering.  
  arrange(name, noon_day) %>% 
  group_by(name, noon_day) %>% 
  mutate(sum_daily_capture = sum(mus_percent_capture)) %>% # for each day sum percent capture pre-penalty application. 
  group_by(name,noon_day) %>% 
  mutate(disc_col = paste0(name, "_", noon_day, "_", sum_daily_capture)) %>% #Need to do this to drop repeated rows. 
  distinct(disc_col, .keep_all = TRUE) %>% # drop repeated daily sums on days with two zone rows. 
  mutate(sum_daily_capture_penalty = sum(sum_daily_capture+penalty)) %>% 
  group_by(name) %>% 
  mutate(csum_daily_capture_penalty = cumsum(sum_daily_capture_penalty))

## Cleaning to remove filled in entries for triaged mice. 
df6 <- df5 %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

df7 <- merge(df6, meta, by = "name", all = FALSE) # bring in metadata

df8 <- df7 %>%
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  dplyr::select(trial.x, strain_sex, strain,sex, name, noon_day, sum_daily_capture, penalty, sum_daily_capture_penalty, csum_daily_capture_penalty) %>% 
  dplyr::rename(trial = trial.x) %>% 
  mutate(label = if_else(noon_day == max(noon_day), as.character(name), NA_character_))

# output csv
write.csv(df8, "data/priority_access_scores_males.csv")

## lineplot
(p <- ggplot(df8, aes(x=noon_day, y=csum_daily_capture_penalty, group = name, color = name)) + #y=csum_adj_mus_percent_capture_score
    geom_line(size =1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(1,10,by=1), limits = c(1,10)) +
    scale_y_continuous(breaks = seq(-10,20, by = 5), limits = c(-10,20)) +
    xlab("Day") +
    ylab("Cumulative zone priority access score") +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "none") 
)
# ggsave(p, filename = "output/move_data_M_priority_access_scores.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_M_priority_access_scores.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# Day 10 Cumulative Sum Priority Access Score and testing for multi-modality
(p <- ggplot(filter(df8, noon_day == 10), aes(x=csum_daily_capture_penalty, group = strain_sex, fill = strain_sex)) + 
    geom_density(adjust = 0.25,alpha = 0.8) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_y_continuous(breaks = seq(0,0.2,by=0.05), limits = c(0, 0.18)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day 10 Priority Access Score Distribution") +
    ylab("Density") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.75) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7), 
          legend.key.size = unit(0.25, 'cm'),
          legend.position = "top")
)
# ggsave(p, filename = "output/move_data_M_priority_access_scores_density.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_M_priority_access_scores_density.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

#multimode
library(multimode)

df9 <- df8 %>% ## change filtering above for males and females. 
  # filter(noon_day == 10, strain_sex == "C57-M")
  filter(noon_day == 10, strain_sex == "NYOB-M")
  
modetest(df9$csum_daily_capture_penalty, mod0 = 1, method = "ACR", B=500,submethod=1,n=NULL,tol=NULL) ## run once. 
# modetree(df9$csum_daily_capture_penalty)
# modeforest(df9$csum_daily_capture_penalty)

#Harington's Diptest for multimodality. outdated. 
library(diptest)
df9 <- df8 %>% 
  filter(noon_day == 10) 
dip.test(df9$csum_daily_capture_penalty)
dS <- (dip(df9$csum_daily_capture_penalty, full.result = TRUE))
# plot(dS)

# PAS_females ---------------------------------------------------------------
df <- move_data %>% 
  filter(sex  == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, noon_day, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(trial_day_zone = paste0(trial, "_",noon_day, "_", zone)) %>% 
  dplyr::rename(total_duration_min = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  dplyr::select(trial_day_zone, total_duration_min)

## individual level zone duration usage per day. 
df1 <- move_data %>% 
  filter(sex  == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  mutate(trial_day_zone = paste0(trial, "_", noon_day, "_", zone)) %>% 
  group_by(name, trial_day_zone, noon_day, zone) %>% 
  tally(sum(duration_min)) %>%
  dplyr::rename(mus_duration_min = n)#specify dplyr due to conflict

df2 <- merge(df1, df, by = "trial_day_zone", all = TRUE) # bring in rest of males that did not win any days. 

df3 <- merge(df2, meta, by = "name", all = FALSE) # bring in metadata

df4 <- df3 %>%
  dplyr::select(trial, name, noon_day, zone, mus_duration_min, total_duration_min)

#summing daily adjusted scores and taking single number to avoid guessing of slice_max when scores vacillate in the negative range. 
df5 <- df4 %>% 
  mutate(mus_percent_capture = (mus_duration_min / total_duration_min)) %>% 
  group_by(trial, name, noon_day) %>% 
  mutate(penalty = if(any(mus_percent_capture > 0.5)) 0 else -1) %>% # PENALTY #1: If on any day you dont capture greater than 50% for any zone, take off 1 point 
  group_by(trial, name) %>%
  complete(noon_day = 1:10, fill = list(penalty = -1, mus_percent_capture = 0)) %>% # PENALTY #2: Not observed at all penalty. Doesnt effect anyone after filtering.  
  arrange(name, noon_day) %>% 
  group_by(name, noon_day) %>% 
  mutate(sum_daily_capture = sum(mus_percent_capture)) %>% # for each day sum percent capture pre-penalty application. 
  group_by(name,noon_day) %>% 
  mutate(disc_col = paste0(name, "_", noon_day, "_", sum_daily_capture)) %>% #Need to do this to drop repeated rows. 
  distinct(disc_col, .keep_all = TRUE) %>% # drop repeated daily sums on days with two zone rows. 
  mutate(sum_daily_capture_penalty = sum(sum_daily_capture+penalty)) %>% 
  group_by(name) %>% 
  mutate(csum_daily_capture_penalty = cumsum(sum_daily_capture_penalty))

## Cleaning to remove filled in entries for triaged mice. 
df6 <- df5 %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

df7 <- merge(df6, meta, by = "name", all = FALSE) # bring in metadata

df8 <- df7 %>%
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  dplyr::select(trial.x, strain_sex, strain,sex, name, noon_day, sum_daily_capture, penalty, sum_daily_capture_penalty, csum_daily_capture_penalty) %>% 
  dplyr::rename(trial = trial.x) %>% 
  mutate(label = if_else(noon_day == max(noon_day), as.character(name), NA_character_))

# output csv
write.csv(df8, "data/priority_access_scores_females.csv")

## lineplot
(p <- ggplot(df8, aes(x=noon_day, y=csum_daily_capture_penalty, group = name, color = name)) + #y=csum_adj_mus_percent_capture_score
    geom_line(size =1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(1,10,by=1), limits = c(1,10)) +
    scale_y_continuous(breaks = seq(-10,20, by = 5), limits = c(-10,20)) +
    xlab("Day") +
    ylab("Cumulative zone priority access score") +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "none") 
)
# ggsave(p, filename = "output/move_data_F_priority_access_scores.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_F_priority_access_scores.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# Day 10 Cumulative Sum Priority Access Score and testing for multi-modality
(p <- ggplot(filter(df8, noon_day == 10), aes(x=csum_daily_capture_penalty, group = strain_sex, fill = strain_sex)) + 
    geom_density(adjust = 0.25,alpha = 0.8) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_y_continuous(breaks = seq(0,0.2,by=0.05), limits = c(0, 0.18)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day 10 Priority Access Score Distribution") +
    ylab("Density") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.75) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7), 
          legend.key.size = unit(0.25, 'cm'),
          legend.position = "top")
)
# ggsave(p, filename = "output/move_data_F_priority_access_scores_density.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_F_priority_access_scores_density.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

#multimode
library(multimode)

df9 <- df8 %>% ## change filtering above for males and females. 
  
  # filter(noon_day == 10, strain_sex == "C57-M")
  # filter(noon_day == 10, strain_sex == "NYOB-M")
  # filter(noon_day == 10, strain_sex == "C57-F") 
  # filter(noon_day == 10, strain_sex == "NYOB-F")
  
  modetest(df9$csum_daily_capture_penalty, mod0 = 1, method = "ACR", B=500,submethod=1,n=NULL,tol=NULL) ## run once. 
# modetree(df9$csum_daily_capture_penalty)
# modeforest(df9$csum_daily_capture_penalty)

#Harington's Diptest for multimodality. outdated. 
library(diptest)
df9 <- df8 %>% 
  filter(noon_day == 10) 
dip.test(df9$csum_daily_capture_penalty)
dS <- (dip(df9$csum_daily_capture_penalty, full.result = TRUE))
# plot(dS)

# PAS_choose ---------------------------------------------------------------
# Get priority access scores from move_data
## Choose your target population from which PAS values will be calculated through filtering. 
## population duration time in zones per day. 
df <- move_data %>% 
  filter(sex  == "M") %>%
  # filter(sex  == "F") %>%
  # filter(trial == "T001") %>%
  # filter(strain == "C57") %>%
  # filter(strain == "NYOB") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, noon_day, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(trial_day_zone = paste0(trial, "_",noon_day, "_", zone)) %>% 
  dplyr::rename(total_duration_min = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  dplyr::select(trial_day_zone, total_duration_min)

## individual level zone duration usage per day. 
df1 <- move_data %>% 
  filter(sex  == "M") %>%
  # filter(sex  == "F") %>%
  # filter(trial == "T001") %>%
  # filter(strain == "C57") %>%
  # filter(strain == "NYOB") %>%
  mutate(duration_min = duration_s / 60) %>% 
  mutate(trial_day_zone = paste0(trial, "_", noon_day, "_", zone)) %>% 
  group_by(name, trial_day_zone, noon_day, zone) %>% 
  tally(sum(duration_min)) %>%
  dplyr::rename(mus_duration_min = n)#specify dplyr due to conflict

df2 <- merge(df1, df, by = "trial_day_zone", all = TRUE) # bring in rest of males that did not win any days. 

df3 <- merge(df2, meta, by = "name", all = FALSE) # bring in metadata

df4 <- df3 %>%
  dplyr::select(trial, name, noon_day, zone, mus_duration_min, total_duration_min)

#summing daily adjusted scores and taking single number to avoid guessing of slice_max when scores vacillate in the negative range. 
df5 <- df4 %>% 
  mutate(mus_percent_capture = (mus_duration_min / total_duration_min)) %>% 
  group_by(trial, name, noon_day) %>% 
  mutate(penalty = if(any(mus_percent_capture > 0.5)) 0 else -1) %>% # PENALTY #1: If on any day you dont capture greater than 50% for any zone, take off 1 point 
  group_by(trial, name) %>%
  complete(noon_day = 1:10, fill = list(penalty = -1, mus_percent_capture = 0)) %>% # PENALTY #2: Not observed at all penalty. Doesnt effect anyone after filtering.  
  arrange(name, noon_day) %>% 
  group_by(name, noon_day) %>% 
  mutate(sum_daily_capture = sum(mus_percent_capture)) %>% # for each day sum percent capture pre-penalty application. 
  group_by(name,noon_day) %>% 
  mutate(disc_col = paste0(name, "_", noon_day, "_", sum_daily_capture)) %>% #Need to do this to drop repeated rows. 
  distinct(disc_col, .keep_all = TRUE) %>% # drop repeated daily sums on days with two zone rows. 
  mutate(sum_daily_capture_penalty = sum(sum_daily_capture+penalty)) %>% 
  group_by(name) %>% 
  mutate(csum_daily_capture_penalty = cumsum(sum_daily_capture_penalty))

## Cleaning to remove filled in entries for triaged mice. 
df6 <- df5 %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

df7 <- merge(df6, meta, by = "name", all = FALSE) # bring in metadata

df8 <- df7 %>%
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  dplyr::select(trial.x, strain_sex, strain,sex, name, noon_day, sum_daily_capture, penalty, sum_daily_capture_penalty, csum_daily_capture_penalty) %>% 
  dplyr::rename(trial = trial.x) %>% 
  mutate(label = if_else(noon_day == max(noon_day), as.character(name), NA_character_))

# output csv
# write.csv(df8, "data/priority_access_scores_males.csv")
# write.csv(df8, "data/priority_access_scores_females.csv")

# Priority Access Score Line Plot
# df8 <- read_excel("SM_data.xlsx", sheet = "Fig2e-f") #note that this is combined male and female data. you will need to filter by sex to produce the graphs below.

(p <- ggplot(df8, aes(x=noon_day, y=csum_daily_capture_penalty, group = name, color = name)) + #y=csum_adj_mus_percent_capture_score
    geom_line(size =1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(1,10,by=1), limits = c(1,10)) +
    scale_y_continuous(breaks = seq(-10,20, by = 5), limits = c(-10,20)) +
    xlab("Day") +
    ylab("Cumulative zone priority access score") +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "none") 
)
# ggsave(p, filename = "output/move_data_M_priority_access_scores.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_M_priority_access_scores.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")
# ggsave(p, filename = "output/move_data_F_priority_access_scores.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_F_priority_access_scores.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")



# Day 10 Cumulative Sum Priority Access Score and testing for multi-modality
(p <- ggplot(filter(df8, noon_day == 10), aes(x=csum_daily_capture_penalty, group = strain_sex, fill = strain_sex)) + 
    geom_density(adjust = 0.25,alpha = 0.8) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_y_continuous(breaks = seq(0,0.2,by=0.05), limits = c(0, 0.18)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day 10 Priority Access Score Distribution") +
    ylab("Density") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.75) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7), 
          legend.key.size = unit(0.25, 'cm'),
          legend.position = "top")
)
# ggsave(p, filename = "output/move_data_M_priority_access_scores_density.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_M_priority_access_scores_density.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")
# ggsave(p, filename = "output/move_data_F_priority_access_scores_density.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/move_data_F_priority_access_scores_density.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

#multimode
library(multimode)

df9 <- df8 %>% ## change filtering above for males and females. 
  
  # filter(noon_day == 10, strain_sex == "C57-M")
  # filter(noon_day == 10, strain_sex == "NYOB-M")
  # filter(noon_day == 10, strain_sex == "C57-F") 
  # filter(noon_day == 10, strain_sex == "NYOB-F")
  
  modetest(df9$csum_daily_capture_penalty, mod0 = 1, method = "ACR", B=500,submethod=1,n=NULL,tol=NULL) ## run once. 
# modetree(df9$csum_daily_capture_penalty)
# modeforest(df9$csum_daily_capture_penalty)

#Harington's Diptest for multimodality. outdated. 
library(diptest)
df9 <- df8 %>% 
  filter(noon_day == 10) 
dip.test(df9$csum_daily_capture_penalty)
dS <- (dip(df9$csum_daily_capture_penalty, full.result = TRUE))
# plot(dS)


