# 4_analyze_ALLTRIAL_MOVEBOUT
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# LOAD PACKAGES, FUNCTIONS, DATA -----------------------------------------------------------
library(tidyverse)
library(data.table)
library(readxl)
library(gganimate)
library(ggpubr)
library(Hmisc)
library(hrbrthemes)
library(gifski)
library(av)
library(reshape)
library(lubridate)
library(plot.matrix)
library(transformr)
library(rstatix)
library(psych)
library(asnipe)
library(igraph)
library(emmeans)
library(sjstats)
library(lme4)
library(lmerTest)
library(MuMIn)
library(broom)
library(knitr)
library(ggfortify)
library(svglite)
library(plotrix)
library(sjPlot)
library(jtools)
library(sjPlot)
library(effects)
library(ggeffects)
library(ggrepel)

wd <- setwd("C:/Users/Caleb/Box/0_CV_Shared_Analyses/7_LID_2020/RFID_analysis_v10")
output_fp <- paste("C:/Users/caleb/Desktop")
meta <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
# weather <- read_excel("LID_2020_metadata.xlsx", sheet = 2, skip = 0)
# weather$Weather_Time <- as.POSIXct(weather$Weather_Time, format = "%m/%d/%Y %H:%M") 
data <- as.data.frame(fread("LID_2020_ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE))
data <- subset(data, select = -c(V1))
data$field_time <- as.POSIXct(data$field_time, format="%Y-%m-%d %H:%M:%OS") # CONVERT TIME FORMAT
data$field_time_STOP <- as.POSIXct(data$field_time_STOP, format="%Y-%m-%d %H:%M:%OS")
data$strain_sex <- paste0(data$strain, "-", data$sex)

#Data Cleaning: removed george from all analyses due to movement between trials
data <- data %>% 
  filter(!(name == "George")) # george crosses from T004 to T005 on Day 2. swaps back and forth before settling as floater in T005 From Day 3 to day 10. Triage completely

# General summary Statistics -----------------------------------------------

## Mean estimated hours in the tubs per trial
df <- data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(trial) %>% 
  tally(sum(duration_hours)) 
sum(df$n)
summary(df$n)
mean(df$n)
stderror(df$n)
#coeffcient of variation
sd(df$n) / mean(df$n) * 100 

## Mean estimated hours in the tub per mouse per night
df <- data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(name, noon_to_noon_day) %>% 
  tally(sum(duration_hours)) 
sum(df$n)
summary(df$n)
mean(df$n)
stderror(df$n)
#coeffcient of variation
sd(df$n) / mean(df$n) * 100 

## Mean resource zones visited per mouse per night 
df <- data %>% # 
  group_by(name, noon_to_noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(name, noon_to_noon_day) %>%
  tally() # get number of distinct zone visits
mean(df$n)
stderror(df$n)

##extra
##quickly graph the quantiative variables
qplot(antenna, n, data = df, colour = sex) +
  stat_summary(fun.data = mean_cl_normal) +
  geom_smooth(method="lm")



# FIG2A_FigS2A_FigS2B_Unique zones visited per day -------------------------------------------------------------------
df <- data %>% # 
  group_by(trial, strain_sex, sex, strain, name, noon_to_noon_day, antenna) %>% 
  # filter(trial == "T007") %>% #Used for supplementary figure S2. 
  tally() %>%  # get number of visits to each zone
  group_by(trial, strain_sex, sex, strain, name, noon_to_noon_day, .drop = FALSE) %>%
  tally() %>% # get number of distinct zone visits
  complete(noon_to_noon_day = 1:10, fill = list(n = 0)) %>% # fill in days where mouse doesnt appear with 0s
  dplyr::rename(unique_zones_visited = n)

# data cleaning. 
df <- df %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & noon_to_noon_day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & noon_to_noon_day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & noon_to_noon_day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & noon_to_noon_day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & noon_to_noon_day >= 10))  

# output csv
write.csv(df, paste(wd, "Figure2A_data.csv", sep = '/'), row.names = FALSE)

mean(df$unique_zones_visited)
std.error(df$unique_zones_visited)

df1 <- df %>% 
  group_by(strain_sex, noon_to_noon_day) %>%
  summarise(mean_zone = mean(unique_zones_visited), sd_zone = sd(unique_zones_visited),count = n(), se_zone = (sd_zone/(sqrt(count))))

# PLOT
(p <- ggplot(df1, aes(x = noon_to_noon_day, y = mean_zone, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_zone - se_zone, ymax = mean_zone + se_zone), width = 0.2) +
    scale_x_continuous(limits = c(1,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab("Day") +
    ylab("# of Unique Zones") +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.21,0.8))
  # legend.position = "top")
  
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent") 

# STATS
summary(df)
df %>% 
  distinct(name) %>%
  group_by(strain_sex) %>%
  tally()

df$trial <- as.factor(df$trial)
df$sex <- as.factor(df$sex)
df$strain <- as.factor(df$strain)
df$noon_to_noon_day <- as.numeric(df$noon_to_noon_day)

mod1 = lmer(unique_zones_visited ~ strain*sex*log(noon_to_noon_day) + (1|trial) + (1+log(noon_to_noon_day)|name), data = df) 
qqnorm(resid(mod1))
qqline(resid(mod1))
anova(mod1)
summary(mod1)
emmeans(mod1, pairwise ~ strain*sex*log(noon_to_noon_day))
AIC(mod1)
ggpredict(mod1, terms = c("strain", "sex", "noon_to_noon_day"), type = "fe") %>% 
  plot()

# options(scipen = 0)
write.table(summary(mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# daily models: change for days 1 through 10. 
d = lmer(n ~ strain*sex + (1|trial), data = subset(df, day == 1)) 
anova(d)
summary(d)
#get daily contrast estimates
emmeans(d, pairwise ~ strain*sex)

#excluding C57 females does negate effect of day for other strains... barely. 
m3 = lmer(unique_zones_visited ~ sex*strain*day + (1|name), data = subset(df, df$strain_sex != "C57-F"))
summary(m3)
anova(m3)



# Fig2B_Cumulative unique zones visited over trial period --------
df <- data %>% 
  select(trial, strain_sex, strain, sex, name, noon_to_noon_day, zone) %>%
  group_by(trial, strain_sex, strain, sex, name, zone) %>% 
  distinct(zone, .keep_all = TRUE) %>% # get distinct zones ever visited, keep associated day it was visited
  group_by(trial, strain_sex,strain, sex, name, noon_to_noon_day) %>% #regroup by day
  tally() %>% #tally unique zones visited per day
  mutate(csum_novel_zones = cumsum(n)) %>%  
  complete(name, noon_to_noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, noon_to_noon_day) %>% 
  fill(csum_novel_zones) %>% ## fill cumulative sum data from last observed day
  select(trial, strain_sex, strain, sex, name, noon_to_noon_day, csum_novel_zones)

# output csv
write.csv(df, paste(wd, "Figure2B_data.csv", sep = '/'), row.names = FALSE)

df1 <- df %>% 
  group_by(strain_sex, noon_to_noon_day) %>% 
  summarise(mean_n = mean(csum_novel_zones), sd_n = sd(csum_novel_zones),count = n(), se_n = (sd_n/(sqrt(count))))

# PLOT
(p <- ggplot(df1, aes(x=noon_to_noon_day, y=mean_n, color = strain_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(breaks = seq(1,11,by=1), limits = c(1,10)) +
    scale_y_continuous(breaks = seq(1,8,by=1), limits =c(1,8)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    xlab("Day") +
    ylab("Cumulative Unique Zones") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.position = "none") +
    guides(color=guide_legend(override.aes=list(fill=NA)))
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent")

# STATS
summary(df)
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
df$csum_novel_zones <- as.numeric(df$csum_novel_zones)

#complex way of doing t-tests but allows you to control for the random effect of trial. 
mod2 = lmer(csum_novel_zones ~ sex*strain + (1|trial), data = subset(df, noon_to_noon_day == 10))
summary(mod2)
anova(mod2)
write.table(summary(mod2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
em.mod2 = emmeans(mod2, pairwise ~ sex*strain)
write.table(em.mod2, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


#logistic regression model of differences in all zones visited by strain*sex interaction
df2 <- df %>% 
  filter(noon_to_noon_day ==10)
df2$every_zone <- ifelse(df2$csum_novel_zones == 8, 1, 0)

mod3 = glmer(every_zone ~ strain_sex + (1|trial), data = df2, family ="binomial") #family binomial because response variable is 0 or 1. 
summary(mod3)
write.table(summary(mod3)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# Fig2C_Estimated minimum distance traveled -------------------------------------------------------------------
df <- data %>% 
  select(trial, strain_sex,sex,strain,name, noon_to_noon_day, zone, zone_x, zone_y, antenna)

df$zone_x[df$zone_x =="A"] <- 3.75
df$zone_x[df$zone_x =="B"] <- 11.25
df$zone_y[df$zone_y =="A"] <- 7.6
df$zone_y[df$zone_y =="B"] <- 15.2
df$zone_y[df$zone_y =="C"] <- 22.8
df$zone_y[df$zone_y =="D"] <- 30.4
df$zone_x <- as.numeric(df$zone_x)
df$zone_y <- as.numeric(df$zone_y)
ids <- unique(df$name)
data_list <- list()
aa = ids[60]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) 
  # delete consecutive repeat antenna hits, only keep rows where antennas change. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$zone)]
  df2$dist <- NA
  n <- nrow(df2)
  if(n==1){
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  } else{
    df2$dist[2:n] <- sqrt((df2$zone_x[2:n] - df2$zone_x[1:n-1]) ^ 2 + (df2$zone_y[2:n] - df2$zone_y[1:n-1]) ^ 2)
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  }
}
df3 <- do.call("rbind", data_list)
summary(df3)
df3$noon_to_noon_day <- as.numeric(df3$noon_to_noon_day)
df4 <- df3 %>% 
  group_by(trial, strain_sex,strain, sex, name, noon_to_noon_day) %>% 
  tally(sum(dist)) %>% 
  complete(noon_to_noon_day = 1:10, fill = list(n = 0)) %>%  # fill in days where mouse doesnt appear with 0s.
  dplyr::rename(dist = n)

# data cleaning. 
df5 <- df4 %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & noon_to_noon_day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & noon_to_noon_day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & noon_to_noon_day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & noon_to_noon_day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & noon_to_noon_day >= 10))  

df6 <- df5 %>% 
  group_by(strain_sex, noon_to_noon_day) %>%
  summarise(mean_n = mean(dist), sd_n = sd(dist),count = n(), se_n = (sd_n/(sqrt(count))))

(p <- ggplot(df6, aes(x = noon_to_noon_day, y = mean_n, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(limits = c(1,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab("Day") +
    ylab("Minimum Distance (m)") +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "none" )
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent")

# STATS  
df7 <- df5 %>% 
  group_by(trial, strain, sex, name, noon_to_noon_day) %>% 
  tally(sum(dist)) %>% 
  group_by(trial, strain, sex, name, noon_to_noon_day) %>% 
  dplyr::rename(sum_dist = n)

# output csv
write.csv(df7, paste(wd, "Figure2C_data.csv", sep = '/'), row.names = FALSE)

df7$trial <- as.factor(df7$trial)
df7$strain <- as.factor(df7$strain)
df7$sex <- as.factor(df7$sex)
df7$name <- as.factor(df7$name)

mod1 = lmer(sum_dist ~ strain*sex*log(noon_to_noon_day) + (1|trial) + (log(noon_to_noon_day)|name), data = df7) 
mod2 = lmer(sum_dist ~ strain*sex*log(noon_to_noon_day) + (log(noon_to_noon_day)|name), data = df7) 
mod3 = lmer(sum_dist ~ strain*sex*log(noon_to_noon_day) + (1|trial) + (1|name), data = df7) 
AIC(mod1,mod2, mod3)
summary(mod1)
anova(mod1)
emmeans(mod1, pairwise ~ strain*sex*log(noon_to_noon_day))
write.table(summary(mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# daily models: change for days 1 through 10 and report
d = lmer(sum_dist ~ strain*sex + (1|trial), data = subset(df7, noon_to_noon_day == 10)) 
em.d = emmeans(d, pairwise ~ strain*sex)
write.table(em.d, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# Fig2D_Percent time spent in most occupied zone -----------------------
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) %>% 
  group_by(name) %>% 
  slice_max(percent_time, n = 2) %>% 
  mutate(rank_order = rank(desc(percent_time))) %>% 
  complete(rank_order = 1:2, fill = list(n = 0, percent_time = 0)) %>% 
  mutate(rank_order = as.factor(rank_order))
df <- merge(df, meta, by  = "name") 

df1 <- df %>% 
  select(trial, name, strain, sex, rank_order, antenna, n, percent_time) %>% 
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  filter(rank_order == 1) %>% 
  dplyr::rename(duration_min = n)

# output csv
write.csv(df1, paste(wd, "Figure2D_data.csv", sep = '/'), row.names = FALSE)

table(round(df1$percent_time, 6), df1$strain_sex)
mode(df1$percent_time)

(p <- ggplot(df1, aes(x=strain_sex, y=percent_time, fill = strain_sex)) + 
    geom_violin(width=1, alpha = 1) +
    geom_boxplot(width=0.1, color = "white", alpha=0.1, size = 0.7) +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("red", "blue", "green", "darkgrey")) +
    xlab("") +
    ylab("Percent Total Time") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8, face = "bold"),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"),
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.position = "none")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5)

# STATS
summary(df1)
df1[sapply(df1, is.character)] <- lapply(df1[sapply(df1, is.character)],as.factor)

## there is a theoretical reason for going with mod2 as it is proportion data. logit transformations cannot be used if there are 0 or 1s
## in the case where there are 0s or 1s in proportion data, you need to do an arcsin transformatoin of the percentage data. 
## if you do the asin transformation, you need to UNDO it when you want to report the estimate values in the text
## to UNDO this, you need to do the sin(estimate)^2 
## i.e. sin(0.94)^2
m2 = lmer(asin(sqrt(percent_time)) ~ strain_sex + (1|trial), data = df1) # random effects
summary(m2)
anova(m2)
emmeans(m2, pairwise ~ strain_sex)
write.table(summary(m2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

#C57-F untransformed estimate
sin(0.94079)^2
#C57-M untransformed estimate
sin(0.26396+0.94079)^2
#OB-F untransformed estimate
sin(0.25090 +0.94079)^2
#OB-M untransformed estimate
sin(0.34030+0.94079)^2

# Fig2E-F_cumulative Proportion territory captures over trial period line plot -----------------------
# options(scipen=999)

#get total time spent in a zone per day per trial for each MALE or FEMALE ONLY
df <- data %>% 
  # filter(trial == "T001") %>%
  # filter(strain == "C57") %>%
  # filter(strain == "NYOB") %>%
  # filter(sex  == "M") %>%
  filter(sex  == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, noon_to_noon_day, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(trial_day_antenna = paste0(trial, "_",noon_to_noon_day, "_", antenna)) %>% 
  dplyr::rename(total_duration_min = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  select(trial_day_antenna, total_duration_min)

df1 <- data %>% 
  # filter(trial == "T001") %>%
  # filter(strain == "C57") %>%
  # filter(strain == "NYOB") %>%
  # filter(sex  == "M") %>%
  filter(sex  == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  mutate(trial_day_antenna = paste0(trial, "_", noon_to_noon_day, "_", antenna)) %>% 
  group_by(name, trial_day_antenna, noon_to_noon_day, antenna) %>% 
  tally(sum(duration_min)) %>%
  dplyr::rename(mus_duration_min = n)#specify dplyr due to conflict

df2 <- merge(df1, df, by = "trial_day_antenna", all = TRUE) # bring in rest of males that did not win any days. 

df3 <- merge(df2, meta, by = "name", all = FALSE) # bring in metadata

df4 <- df3 %>%
  select(trial, name, noon_to_noon_day, antenna, mus_duration_min, total_duration_min)

#testing summing daily adjusted scores and taking single number to avoid guessing of slice_max when scores vacillate in the negative range. 
df5 <- df4 %>% 
  mutate(mus_percent_capture = (mus_duration_min / total_duration_min)) %>% 
  group_by(trial, name, noon_to_noon_day) %>% 
  mutate(penalty = if(any(mus_percent_capture > 0.5)) 0 else -1) %>% # PENALTY #1: If on any day you dont capture greater than 50% for any zone, take off 1 whole point for each zone you failed to reach >50%. 
  group_by(trial, name) %>%
  complete(noon_to_noon_day = 1:10, fill = list(penalty = -1, mus_percent_capture = 0)) %>% # PENALTY #2: Not observed at all penalty. Doesnt effect anyone after filtering.  
  arrange(name, noon_to_noon_day) %>% 
  group_by(name, noon_to_noon_day) %>% 
  mutate(sum_daily_capture = sum(mus_percent_capture)) %>% # for each day sum percent capture pre-penalty application. 
  group_by(name,noon_to_noon_day) %>% 
  mutate(disc_col = paste0(name, "_", noon_to_noon_day, "_", sum_daily_capture)) %>% #Need to do this to drop repeated rows. 
  distinct(disc_col, .keep_all = TRUE) %>% # drop repeated daily sums on days with two zone rows. 
  mutate(sum_daily_capture_penalty = sum(sum_daily_capture+penalty)) %>% 
  group_by(name) %>% 
  mutate(csum_daily_capture_penalty = cumsum(sum_daily_capture_penalty))

## CLEANING
df6 <- df5 %>% 
  filter(!(name == "George")) %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & noon_to_noon_day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & noon_to_noon_day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & noon_to_noon_day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & noon_to_noon_day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & noon_to_noon_day >= 10))

df7 <- merge(df6, meta, by = "name", all = FALSE) # bring in metadata

df8 <- df7 %>%
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  select(trial.x, strain_sex, strain,sex, name, noon_to_noon_day, sum_daily_capture, penalty, sum_daily_capture_penalty, csum_daily_capture_penalty) %>% 
  dplyr::rename(trial = trial.x) %>% 
  mutate(label = if_else(noon_to_noon_day == max(noon_to_noon_day), as.character(name), NA_character_))

# output csv
# write.csv(df8, paste0(output_fp, "/", "Figure2E-F_male_data.csv"), row.names = FALSE)
write.csv(df8, paste0(output_fp, "/", "Figure2E-F_female_data.csv"), row.names = FALSE)

# line plot
(p <- ggplot(df8, aes(x=noon_to_noon_day, y=csum_daily_capture_penalty, group = name, color = name)) + #y=csum_adj_mus_percent_capture_score
    geom_line(size =1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(1,10,by=1), limits = c(1,10)) +
    xlab("Day") +
    ylab("Cumulative Territory Score") +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "none") 
  # geom_label_repel(aes(label = label),
  #                  point.padding = 0.2,
  #                  size = 4,
  #                  force = 3,
  #                  # arrow = arrow(length = unit(0.015, "npc")),
  #                  # color = "black",
  #                  nudge_x = 2,
  #                  nudge_y = 5,
  #                  na.rm = TRUE,
  #                  direction = "both",
  #                  hjust = 0,
  #                  segment.alpha = 0.5)
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent") #main figs

# density plot
(p <- ggplot(subset(df8, noon_to_noon_day == 10), aes(x=csum_daily_capture_penalty, group = strain_sex, fill = strain_sex)) + 
    geom_density(adjust = 0.25,alpha = 0.4) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("red", "blue", "green", "darkgrey")) +
    xlab("Territory Capture Score") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.75) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6), 
          legend.key.size = unit(0.25, 'cm'),
          legend.position = "top")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=1.5, bg = "transparent") #main figs


library(diptest)
df9 <- df8 %>% 
  filter(noon_to_noon_day == 10)   #, strain == "NYOB", )
dip.test(df9$csum_daily_capture_penalty)

dS <- (dip(df9$csum_daily_capture_penalty, full.result = TRUE))
plot(dS)




# Fig2G_PCA of space use and movement patterns  --------------------------------------------------------------
# total_zone_time_h
df <- data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(name) %>% 
  tally(sum(duration_hours)) 
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# mean_bout_visit_s
df <- data %>% 
  group_by(name) %>% 
  tally(mean(duration_s)) 
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_unique_zones_per_trial
df <- data %>% 
  group_by(name,antenna) %>% 
  tally() %>% 
  dplyr::count(name)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# mean_dist_zones_per_night
df <- data %>% # 
  group_by(name, noon_to_noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(name, noon_to_noon_day) %>%
  tally() %>% # get number of distinct zone visits
  select(name, n) %>% 
  summarize_all(funs(mean))
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# %_time_rank_1_zone
options(scipen=999)
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) %>% 
  group_by(name) %>% 
  slice(which.max(percent_time)) %>% 
  select(name, percent_time)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# %_time_rank_2_zone
options(scipen=999)
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) %>% 
  group_by(name) %>% 
  slice(-which.max(percent_time)) %>% 
  slice(which.max(percent_time)) %>% 
  select(name, percent_time)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_time_rank_1_zone
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  arrange(desc(n)) %>% 
  group_by(name) %>% 
  slice(which.max(n)) %>% 
  select(name, n)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_time_rank_2_zone
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  arrange(desc(n)) %>% 
  group_by(name) %>% 
  slice(-which.max(n)) %>% 
  slice(which.max(n)) %>% 
  select(name, n)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

## distance variables
df <- data %>% 
  select(trial, strain_sex,sex,strain,name, noon_to_noon_day, zone, zone_x, zone_y, antenna)

df$zone_x[df$zone_x =="A"] <- 3.75
df$zone_x[df$zone_x =="B"] <- 11.25
df$zone_y[df$zone_y =="A"] <- 7.6
df$zone_y[df$zone_y =="B"] <- 15.2
df$zone_y[df$zone_y =="C"] <- 22.8
df$zone_y[df$zone_y =="D"] <- 30.4

df$zone_x <- as.numeric(df$zone_x)
df$zone_y <- as.numeric(df$zone_y)
summary(df)

ids <- unique(df$name)
data_list <- list()
aa = ids[60]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) 
  # delete consecutive repeat antenna hits, only keep rows where antennas change. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$zone)]
  df2$dist <- NA
  n <- nrow(df2)
  
  if(n==1){
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  } else{
    df2$dist[2:n] <- sqrt((df2$zone_x[2:n] - df2$zone_x[1:n-1]) ^ 2 + (df2$zone_y[2:n] - df2$zone_y[1:n-1]) ^ 2)
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  }
}

df3 <- do.call("rbind", data_list)
summary(df3)
df3$noon_to_noon_day <- as.numeric(df3$noon_to_noon_day)

df4 <- df3 %>% 
  group_by(trial, strain_sex,strain, sex, name, noon_to_noon_day) %>% 
  tally(sum(dist)) %>% 
  complete(noon_to_noon_day = 1:10, fill = list(n = 0)) %>%  # fill in days where mouse doesnt appear with 0s.
  dplyr::rename(dist = n)

# data cleaning. 
df5 <- df4 %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & noon_to_noon_day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & noon_to_noon_day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & noon_to_noon_day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & noon_to_noon_day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & noon_to_noon_day >= 10))  

# avg_min_dist_traveled_per_night
df4 <- df3 %>% 
  group_by(name, noon_to_noon_day) %>% 
  tally(mean(dist_moved)) %>% 
  complete(noon_to_noon_day = 1:10, fill = list(n = 0)) # fill in days where mouse doesnt appear with 0s
df4 <- subset(df4, n > 0 | (name != "Anubis" & name != "Hare" & name != "Gilmore" & name != "Isis" & name != "Ray" & name != "Rose"))## if mouse has all 0s consecutively till the end, delete 0 rows for that animal. 
df5 <- df4 %>% 
  group_by(name) %>%
  tally(mean(n))
write.table(df5, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD 

# total min distance traveled
df4 <- df3 %>% 
  group_by(strain_sex, name, noon_to_noon_day) %>% 
  tally(sum(dist_moved)) %>% 
  group_by(name) %>%
  tally(sum(n))
write.table(df4, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD 

## GRAPH
df <- read.csv("Figure2G_data.csv", row.names = 1) # note, George removed from the analysis
pca <- prcomp(df, scale = TRUE) # always scale the data. 
plot(pca$x[,1], pca$x[,2])
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal Component", ylab = "Percent Variation")
biplot(pca) 
loading_scores <- pca$rotation[,1]
move_scores <- abs(loading_scores)
move_scores_ranked <- sort(move_scores, decreasing = TRUE)
top_moves <- names(move_scores_ranked)
pca$rotation[top_moves,1]


df2 <- data.frame(name = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
df3 <- left_join(df2, meta) 
df4 <- df3 %>% 
  mutate(strain_sex = paste(strain, sex, sep= "-")) %>% 
  select(strain_sex, strain,sex, name, X, Y)

(p <- ggplot(df4, aes(x = X, y = Y, color = strain_sex, shape = strain_sex)) +
    geom_point(aes(shape = factor(strain_sex))) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.key.size = unit(0.2, 'cm'), 
          legend.text = element_text(size=7),
          legend.position = "top")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent")


# FigS1B_Number of zone visits and time in zone correlation -----------------------
# number of zone visits
df <- data %>% 
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(trial, strain, sex, name, antenna) %>% 
  tally()

# duration of time spent in zone
df2 <- data %>% 
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(trial, strain, sex, name, antenna) %>% 
  tally(sum(duration_min))

df3 <- as.data.frame(cbind(df, df2$n))
df3$strain_sex <- paste0(df3$strain,"-",df3$sex)
colnames(df3) <- c("trial", "strain", "sex", "name", "antenna", "num_visits", "total_time_min", "strain_sex")
#write data frame
write.csv(df3, paste(wd, "FigureS1B_data.csv", sep = '/'), row.names = FALSE)

(p <- ggplot(data = df3, aes(x = num_visits, y = total_time_min, color = strain_sex)) +
    geom_point() +
    scale_color_brewer(palette = "PuOr") + 
    geom_smooth(method="lm") + 
    xlab("Number of zone visits") +
    ylab("Time spent in zone (min)") +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "black")) +
    stat_cor(aes(color = strain_sex), label.x = 4, method = "pearson", p.accuracy = 0.001) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
)
#save
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=5, height=3)

##stats
# m1 = lmer(total_time_min ~ num_visits*sex*strain + (1|trial) + (num_visits|name), data = df3)
# qqnorm(resid(m1))
# qqline(resid(m1))
# anova(m1)
# summary(m1)
# emmeans(m1, pairwise ~ num_visits*sex*strain)
# AIC(m1)



# FigS1C-F_Number of zone Visits per night heatmap ---------------------------
# GP heat map scale should be set to the max visits for any given category. Suri = 184 visits. 
df <- data %>%
  # filter(Trial == "T001", sex == "M") %>%
  # filter(Trial == "T001", sex == "F") %>%
  # filter(Trial == "T002", sex == "M") %>%
  # filter(Trial == "T002", sex == "F") %>%
  # filter(Trial == "T003", sex == "M") %>%
  # filter(Trial == "T003", sex == "F") %>%
  # filter(Trial == "T004", sex == "M") %>%
  filter(Trial == "T004", sex == "F") %>%
  # filter(Trial == "T005", sex == "M") %>%
  # filter(Trial == "T005", sex == "F") %>%
  # filter(Trial == "T006", sex == "M") %>%
  # filter(Trial == "T006", sex == "F") %>%
  # filter(Trial == "T007", sex == "M") %>%
  # filter(Trial == "T007", sex == "F") %>%
  group_by(name, noon_to_noon_day, antenna) %>%
  tally() #number of visits

#check maximum visit number! adjust GP heatmap settings accordingly. 
# Suri = 123 visits

ids <- unique(df$name)#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
idlist <- list() # NOTE THAT DATA FRAMED WILL BE IN ORDER OF THIS LIST
daylist<- list()
flag=1
# aa = ids[1]
for(aa in ids[1:length(ids)]){ # LOOP THROUGH EACH INDIVIDUAL AND PULL OUT NUMBER OF VISITS PER UNIQUE zone AND PUT INTO 2X4 GRID THAT LOCALIZES TO THE 
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 4, ncol = 2)) # ENCLOSURE SETUP. PUT EACH NIGHT OF ACTIVITY TO THE RIGHT FOR 10 NIGHTS. COPY AND PASTE DIRECLY INTO PRISM. 
  df1 <- subset(df, name == aa)
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY zone
    if(nrow(df2) == 0){
      stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
    } else {
      
      for(cc in 1:nrow(df2)){
        if(df2$antenna[cc] == 1){
          stats[4,1] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 2){
          stats[4,2] <-print(df2$n[cc])
        } else if(df2$antenna[cc] == 3){
          stats[3,1] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 4){
          stats[3,2] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 5){
          stats[2,1] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 6){
          stats[2,2] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 7){
          stats[1,1] <- print(df2$n[cc])
        } else if (df2$antenna[cc] == 8){
          stats[1,2] <- print(df2$n[cc])
        } else {print("fuck")}
        
      }
    }
    daylist[[bb]] <- stats
    stats <- data.frame(matrix(0, nrow = 4, ncol = 2)) #CHANGE FROM NA 'S TO 0 'S
  }
  master_class <- do.call("cbind",daylist) #THIS THROWS THE ERROR
  idlist[[flag]] <- master_class
  flag=flag+1
}
master_SASS <- do.call("rbind",idlist) # RBIND DATA FOR EACH INDIVIDUAL
# View(master_SASS)
head(master_SASS)
write.table(master_SASS, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  
ids # PRINT ORDER OF THE DATA SET, COPY INTO GRAPHPAD

# FigS2C_Cumulative estimated time spent in resource zones --------
df <- data %>% 
  select(strain, sex, name, noon_to_noon_day,duration_s) %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, noon_to_noon_day) %>% 
  tally(sum(duration_min)) %>% 
  mutate(csum_zone_min = cumsum(n)) %>%  #get cumulative # of novel mice met
  complete(name, noon_to_noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, noon_to_noon_day) %>% 
  fill(csum_zone_min) %>% ## fill cumulative sum data from last observed day
  mutate(group = paste0(strain, "-", sex)) %>% 
  select(group, strain, sex, name, noon_to_noon_day, csum_zone_min)

# data cleaning. 
df <- df %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & noon_to_noon_day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & noon_to_noon_day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & noon_to_noon_day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & noon_to_noon_day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & noon_to_noon_day >= 10))  

# output csv
write.csv(df, paste(wd, "FigureS2C_data.csv", sep = '/'), row.names = FALSE)

df1 <- df %>% 
  filter(group == "C57-M")


## Graph
df1 <- df %>% 
  group_by(group, noon_to_noon_day) %>% 
  summarise(mean_n = mean(csum_zone_min), sd_n = sd(csum_zone_min),count = n(), se_n = (sd_n/(sqrt(count))))

# PLOT
(p <- ggplot(df1, aes(x=noon_to_noon_day, y=mean_n, color = group)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(breaks = seq(1,11,by=1), limits = c(1,10.2)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "black")) +
    xlab("Day") +
    ylab("Cumulative time (min)") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          # legend.position = c(0.15,0.8))
          legend.position = "top") +
    guides(color=guide_legend(override.aes=list(fill=NA)))
)

ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=4, height=3) #supp figs. 

# Average visit bout duration by sex and strain (in progress) ---------------------------








# OLD CODE_DELETE PRIOR TO PUBLICATION ----------------------------------------------------------------
# [GP] Time in zones per night heatmap ------------------------------------------------------------
## SUM THE DURATIONS SPENT BY INDIVIDUAL PER DAY PER zone
# Havent used these because the range of durations on the first day screws up the scale. 
df <- data %>% 
  filter(Trial == "T001", sex == "M") %>%
  # filter(Trial == "T001", sex == "F") %>%
  # filter(Trial == "T002", sex == "M") %>%
  # filter(Trial == "T002", sex == "F") %>%
  # filter(Trial == "T003", sex == "M") %>%
  # filter(Trial == "T003", sex == "F") %>%
  # filter(Trial == "T004", sex == "M") %>%
  # filter(Trial == "T004", sex == "F") %>%
  # filter(Trial == "T005", sex == "M") %>%
  # filter(Trial == "T005", sex == "F") %>%
  # filter(Trial == "T006", sex == "M") %>%
  # filter(Trial == "T006", sex == "F") %>%
# filter(Trial == "T007", sex == "M") %>%
# filter(Trial == "T007", sex == "F") %>%
mutate(duration_min = duration_s / 60) %>% 
  group_by(name, noon_to_noon_day, antenna) %>% 
  # tally(sum(duration_s))
  tally(sum(duration_min))


#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
ids <- unique(df$name)

# CREATE EMPTY LISTS. 
idlist <- list()
daylist<- list()

# LOOP THROUGH EACH INDIVIDUAL AND PULL OUT NUMBER OF VISITS PER UNIQUE zone AND PUT INTO 2X4 GRID THAT LOCALIZES TO THE 
# ENCLOSURE SETUP. PUT EACH NIGHT OF ACTIVITY TO THE RIGHT FOR 10 NIGHTS. COPY AND PASTE DIRECLY INTO PRISM. 

flag=1
aa = ids[1]
for(aa in ids[1:length(ids)]){
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
  df1 <- subset(df, name == aa)
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY zone
    if(nrow(df2) == 0){
      stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
    } else {
      
      for(cc in 1:nrow(df2)){
        if(df2$antenna[cc] == 1){
          stats[4,1] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 2){
          stats[4,2] <-print(df2$n[cc])
        } else if(df2$antenna[cc] == 3){
          stats[3,1] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 4){
          stats[3,2] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 5){
          stats[2,1] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 6){
          stats[2,2] <- print(df2$n[cc])
        } else if(df2$antenna[cc] == 7){
          stats[1,1] <- print(df2$n[cc])
        } else if (df2$antenna[cc] == 8){
          stats[1,2] <- print(df2$n[cc])
        } else {print("ERROR")}
        
      }
    }
    daylist[[bb]] <- stats
    stats <- data.frame(matrix(0, nrow = 4, ncol = 2)) #CHANGE FROM NA 'S TO 0 'S
  }
  master_class <- do.call("cbind",daylist) #THIS THROWS THE ERROR
  idlist[[flag]] <- master_class
  flag=flag+1
}

# RBIND DATA FOR EACH INDIVIDUAL
master_SASS <- do.call("rbind",idlist)
head(master_SASS)
## IN THE FUTURE, CHANGE THE COLUMN nameS
View(master_SASS)

# COPY THE OUTPUT TO THE CLIPBOARD  
write.table(master_SASS, "clipboard", sep="\t", row.names=FALSE)

# PRINT ORDER OF THE DATA SET, COPY INTO GRAPHPAD
ids
write.table(ids, "clipboard", sep="\t", row.names=FALSE)


# [GP, line] Percent time in top-ranked zone per night -----------------------
# https://stackoverflow.com/questions/45341541/group-by-in-dplyr-and-calculating-percentages

options(scipen=999)
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, noon_to_noon_day, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) 


dplyr::count(name, noon_to_noon_day)

group_by(name, noon_to_noon_day, n) %>% 
  transmute(noon_to_noon_day, Percentage = )
mutate(daily_sum = sum(name, noon_to_noon_day, n))


ids <- unique(df$name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) %>% 
    mutate(percent_time = n / sum(n)) %>% 
    arrange(desc(percent_time))
  
  df2 <- data.frame(df1[,7]) #grab percent time
  df2[(nrow(df2)+1):8,] <- 0
  df3 <- df2[1:8,]
  
  data_list[[aa]] <- df3
}
df4 <- do.call("cbind", data_list)
# COPY OUTPUT TO CLIPBOARD. Port into excel, sort, and put into prism. 
write.table(df4, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE)




# [GP, line, IN PROGRESS] % Time in Top Ranked zone Per Night -----------------------
options(scipen=999)

df <- data %>% 
  select(strain, sex, name, noon_to_noon_day, zone) %>%
  group_by(strain, sex, name, zone) %>% 
  distinct(zone, .keep_all = TRUE) %>% # get distinct zones ever visited, keep associated day it was visited
  group_by(strain, sex, name, noon_to_noon_day) %>% #regroup by day
  tally() %>% #tally unique zones visited per day
  mutate(csum_novel_zones = cumsum(n)) %>%  #get cumulative # of novel mice met
  complete(name, noon_to_noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, noon_to_noon_day) %>% 
  fill(csum_novel_zones) %>% ## fill cumulative sum data from last observed day
  mutate(group = paste0(sex, "-", strain)) %>% 
  select(group, strain, sex, name, noon_to_noon_day, csum_novel_zones)

write.csv(df, file = "zone_accumulation_df.csv")
p <- ggplot(df, aes(x=noon_to_noon_day, y=csum_novel_zones, colour = group)) + 
  geom_smooth()
p

ggsave("ouput.png", plot = p, device='png', path = output_fp)


df <- data %>% 
  # filter(strain == "C57", sex == "M") %>%
  # filter(strain == "C57", sex == "F") %>%
  # filter(strain == "NYOB", sex == "M") %>%
  filter(strain == "NYOB", sex == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min))

ids <- unique(df$name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) %>% 
    mutate(percent_time = n / sum(n)) %>% 
    arrange(desc(percent_time))
  
  df2 <- data.frame(df1[,6])
  df2[(nrow(df2)+1):8,] <- 0
  df3 <- df2[1:8,]
  
  data_list[[aa]] <- df3
}
df4 <- do.call("cbind", data_list)
# COPY OUTPUT TO CLIPBOARD. Port into excel, sort, and put into prism. 
write.table(df4, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE)




# [GP, line, NOT USED] Ranked zone duration (min) -----------------------
## create table, copy and paste into Prism. 
df <- data %>% 
  # filter(strain == "C57", sex == "M") %>%
  # filter(strain == "C57", sex == "F") %>%
  # filter(strain == "NYOB", sex == "M") %>%
  # filter(strain == "NYOB", sex == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min))

ids <- unique(df$name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  df1 <- df %>% 
    filter(name == aa) %>% 
    arrange(desc(n))
  
  df2 <- data.frame(df1[,5])
  df2[(nrow(df2)+1):8,] <- 0
  df3 <- df2[1:8,]
  
  data_list[[aa]] <- df3
}

df4 <- do.call("cbind", data_list)

# COPY OUTPUT TO CLIPBOARD. Port into excel, sort, and put into prism. 
write.table(df4, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE)





# [GP, line, NOT USED] Ranked zone visits (freq) -----------------------
## create table, copy and paste into Prism. 
df <- data %>% 
  # filter(strain == "C57", sex == "M") %>%
  # filter(strain == "C57", sex == "F") %>%
  # filter(strain == "NYOB", sex == "M") %>%
  filter(strain == "NYOB", sex == "F") %>%
  group_by(strain, sex, name, antenna) %>% 
  tally()

ids <- unique(df$name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) %>% 
    arrange(desc(n))
  
  df2 <- data.frame(df1[,5])
  df2[(nrow(df2)+1):8,] <- 0
  df3 <- df2[1:8,]
  
  data_list[[aa]] <- df3
}

df4 <- do.call("cbind", data_list)

# COPY OUTPUT TO CLIPBOARD. Port into excel, sort, and put into prism. 
write.table(df4, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE)



# [R, ] strain AND sex zone VISIT DURATION (MIN) VIOLIN plot --------------------------------------
df <- data %>% 
  group_by(strain, sex, name) %>% 
  mutate(duration_hours = duration_s / 3600) %>% 
  tally(sum(duration_hours))

p <- ggplot(df) +
  aes(x = strain, y = n, fill = sex) +
  geom_violin(trim = TRUE) + 
  labs(x = "strain", y = "Total time in zones (h)") +
  theme_classic()
p

## ADD STATS
p + stat_compare_means(method = "t.test")

ggsave("FigS2Astrain by sex zone Visit Duration (min), Violin.svg", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='svg', 
       path = output_fp)






# [GGPLOT] INDIVIDUAL DURATIONS OBSERVED BY TRIAL BALL AND STICK -----------------------------------
trial <- unique(data$Trial)
strain <- c("C57", "C57", "C57", "OB", "OB", "C57", "OB")
## LOOP ACROSS TRIALS
i=trial[1]
flag = 1
for(i in trial[1:length(trial)]) {
  current_trial <- print(i)
  this_strain <-  strain[flag]
  
  df <- data %>% 
    mutate(name_sex = paste(name,sex, sep= "-")) %>% 
    filter(Trial == i) %>%
    group_by(name_sex) %>% 
    tally(sum(duration_s / 60)) %>% 
    arrange(desc(n))
  
  p <- df %>%
    filter(!is.na(n)) %>%
    arrange(n) %>%
    mutate(name_sex=factor(name_sex, name_sex)) %>%
    ggplot(aes(name_sex, n) ) +
    geom_segment( aes(x = name_sex, xend = name_sex, y = 0, yend = n), color="grey") +
    geom_point(size=3, color="#69b3a2") +
    coord_flip() +
    theme_ipsum() +
    theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position="none"
    ) +
    ylab("Time (min)") + 
    ylim(0,2500) +
    xlab("") +
    ggtitle(current_trial, this_strain)
  
  ggsave(filename=paste(current_trial, this_strain, "Summed Mouse zone Duration (min).png", sep = " "), 
         plot=p, width = 6.5, height = 5, dpi = 300, units = "in", device='png', path = output_fp)
  
  flag <- flag +1
}

p

# [PNG] TRIAL zone VISITS, VIOLIN ----------------------------------------------------
df <- data %>% 
  group_by(Trial, strain, name) %>%  
  tally()

p <- ggplot(df) +
  aes(x = Trial, y = n, fill = strain) +
  geom_violin() +
  labs(title = "Trial zone Visits",
       subtitle = "",
       x = "Trial",
       y = "zone Visits",
       caption = "") +
  theme_classic() 
p

p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

ggsave("Trial zone Visits.png", 
       # plot=p, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)


# [PNG] TRIAL zone DURATION (MIN), VIOLIN ----------------------------------------------------
df <- data %>% 
  group_by(Trial, strain, name) %>%  
  mutate(duration_min = duration_s / 60) %>% 
  tally(sum(duration_min))

p <- ggplot(df) +
  aes(x = Trial, y = n, fill = strain) +
  geom_violin() +
  labs(title = "Trial zone Duration (min)",
       subtitle = "",
       x = "Trial",
       y = "Time (min)",
       caption = "") +
  theme_classic() 
p

p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

ggsave("Trial zone Duration (Min), Violin.png", 
       # plot=p, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)



# [PNG] strain BY sex zone VISIT FREQUENCY, VIOLIN --------------------------------------
library(ggplot2)

df <- data %>% 
  group_by(strain, sex, name) %>% 
  tally()

p <- ggplot(df) +
  aes(x = strain, y = n, fill = sex) +
  geom_violin(trim = TRUE) + 
  labs(title = "zone Visits",
       subtitle = "",
       x = "strain",
       y = "Number of zone Visits",
       caption = "") +
  theme_classic() 
p

## ADD STATS
p + stat_compare_means(method = "t.test")

ggsave("strain by sex zone Visit Frequency, Violin.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)





# CREATE DF OF AVERAGE DURATION OF VISITS PER zone PER MALE ---------------

# GET TOTAL DURATION OF ALL VISITS ACROSS ALL TRIALS. 
sum(data$duration_s)

means <- data %>%
  filter(Trial == "T001", sex == "M") %>% 
  group_by(sex, name, noon_to_noon_day, zone) %>%
  dplyr::summarize(Mean = mean(duration_s, na.rm=TRUE)) #average duration of visit. 




# Q: IS THERE A DIFFERENCE IN THE NUMBER OF MALE AND FEMALE VISITS? --------------------------------------------

df <- data %>% 
  group_by(sex, Full_ID) %>% 
  tally()

# COPY OUTPUT TO CLIPBOARD. 
write.table(df, "clipboard", sep="\t", row.names=FALSE)

# PRODUCE SUMMARY STATS
df %>% 
  group_by(sex) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS

# NEXT, QUICKLY PLOT THE DATA TO VISUALIZE
plot <- ggboxplot(df, 
                  x="sex", 
                  y="n", 
                  color = "sex",
                  # palette = c("#00AFBB", "#E7B800"),
                  palette = "jco",
                  add = "jitter",
                  ylab = "RFID Reads",
                  xlab = "sex",
                  # title = "Difference in MF RFID Trials?",
                  # subtitle = "All Trials")
)
plot

# ADD TEST RESULTS TO GGPLOT
plot +
  stat_compare_means(method = "wilcox.test")

### CHECK ASSUMPTIONS
### ASSUMPTION 1: ARE THE TWO SAMPLES INDEPENDENT? POTENTIALLY NOT, AS THEY ARE IN THE SAME PADDOCK. 

### ASSUMPTION 2: DO THE DATA FROM EACH OF THE TWO GROUPS FOLLOWING A NORMAL DISTRIBUTION (PARAMETRIC)? 
# IF YES, USE T.TEST. IF NOT, USE WILCOXIN TEST

#GROUP 1
df1 <- df %>% filter(sex == "M")
ggdensity(df1$n)
ggqqplot(df1$n)
shapiro.test(df1$n) #p-val >0.05 implies data distribution not sig.dif. from a normal distribution. 

# GROUP 2
df1 <- df %>% filter(sex == "F")
ggdensity(df1$n)
ggqqplot(df1$n)
shapiro.test(df1$n) 

#### ASSUMPTION 3: DO THE TWO POPS HAVE THE SAME VARIANCE? USE THE F-TEST
# P-VALUE > 0.05, NO SIGNIFICIANT DIFFERENCE BETWEEN THE VARIANCES OF THE TWO DATA SETS.
var.test(n ~ sex, data = df)

# CHOOSE YOUR TEST. PARAMETRIC (T-TEST) OR NON-PARAMETRIC(WILCOX.TEST)
compare_means(n~sex, data = df, method = "wilcox.test", paired = FALSE)




# PLOT TOTAL # OF VISITS VISITS BY MALES AND FEMALES, MF ----------------------------------------------------

# SUBSAMPLE DATA
df <- data %>% 
  group_by(Trial,sex) %>% 
  tally()
df

plot <- ggbarplot(df, 
                  x="Trial",
                  y="n",
                  fill = "sex",
                  color = "sex",
                  palette = "Paired",
                  position = position_dodge(0.8),
                  ylab = "# OF zone VISITS",
                  xlab = "Trial",
                  title = "# OF zone VISITS")
plot



# DIFFERENCE IN MF RFID READS BY TRIAL? --------
df <- data %>% 
  group_by(Trial, sex, Full_ID) %>% 
  tally()
df

plot <- ggbarplot(df, 
                  x="Trial",
                  y="n",
                  fill = "sex",
                  color = "sex",
                  palette = "Paired",
                  position = position_dodge(0.8),
                  add = c("mean_sd", "jitter"),
                  ylab = "# OF VISITS",
                  xlab = "Trial",
                  # title = "Total RFID Reads Per Trial")
)
plot




# Q: DIFFERENCE IN TOTAL VISITS FOR NYOB AND C57 MALES? Y: -----------------------------

df <- data %>% 
  group_by(Trial,strain, sex, Full_ID) %>% 
  filter(sex == "M") %>% 
  tally()


# PRODUCE SUMMARY STATS
df %>% 
  group_by(strain) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS

# NEXT, QUICKLY PLOT THE DATA TO VISUALIZE
plot <- ggboxplot(df, 
                  x="strain", 
                  y="n", 
                  color = "strain",
                  palette = "jco",
                  add = "jitter",
                  ylab = "# OF zone VISITS",
                  xlab = "sex",
                  title = "Male"
                  # subtitle = ""
)
plot

# ADD TEST RESULTS TO GGPLOT
plot +
  stat_compare_means(method = "wilcox.test")



### ASSUMPTION 1: ARE THE TWO SAMPLES INDEPENDENT? POTENTIALLY NOT, AS THEY ARE IN THE SAME PADDOCK. 

### ASSUMPTION 2: DO THE DATA FROM EACH OF THE TWO GROUPS FOLLOWING A NORMAL DISTRIBUTION (PARAMETRIC)? 
# IF YES, USE T.TEST. IF NOT, USE WILCOXIN TEST
# SHAPIRO TEST: #p-val >0.05 implies data distribution not sig.dif. from a normal distribution.

#GROUP 1
df1 <- df %>% filter(strain == "C57")
ggdensity(df1$n)
ggqqplot(df1$n, 
         title = "C57 Males")
shapiro.test(df1$n)  

#GROUP 2
df1 <- df %>% filter(strain == "NYOB")
ggdensity(df1$n)
ggqqplot(df1$n,
         title = "NYOB Males")
shapiro.test(df1$n) 

#### ASSUMPTION 3: DO THE TWO POPS HAVE THE SAME VARIANCE? USE THE F-TEST
# P-VALUE > 0.05, NO SIGNIFICIANT DIFFERENCE BETWEEN THE VARIANCES OF THE TWO DATA SETS.
var.test(n ~ strain, data = df)

# CHOOSE YOUR TEST. PARAMETRIC (T-TEST) OR NON-PARAMETRIC(WILCOX.TEST)
compare_means(n~strain, data = df, method = "wilcox.test", paired = FALSE)

# DIFFERENCE IN # OF zone VISITS FOR NYOB AND C57 FEMALES? ----------------------

df <- data %>% 
  group_by(Trial,strain, sex, Full_ID) %>% 
  filter(sex == "F") %>% 
  tally()


# PRODUCE SUMMARY STATS
df %>% 
  group_by(strain) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS

# NEXT, QUICKLY PLOT THE DATA TO VISUALIZE
plot <- ggboxplot(df, 
                  x="strain", 
                  y="n", 
                  color = "strain",
                  palette = "jco",
                  add = "jitter",
                  ylab = "# OF VISITS",
                  xlab = "strain",
                  title = "Female",
                  # subtitle = "All Trials")
)
plot

plot +
  stat_compare_means(method = "wilcox.test")

# ADD TEST RESULTS TO GGPLOT
#CREATE COMPARISONS FOR PVALS AND SIG BARS.
comps <- list(c("C57", "NYOB"))

plot +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     label.y.npc = 0.8, 
                     label.x.npc = 0.5,
                     bracket.size = 2,
                     comparisons = comps)

# THERE ARE SOME POTENTIAL OUTLIERS HERE THAT COULD BE REMOVED. 
# https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/



### ASSUMPTION 1: ARE THE TWO SAMPLES INDEPENDENT? POTENTIALLY NOT, AS THEY ARE IN THE SAME PADDOCK. 

### ASSUMPTION 2: DO THE DATA FROM EACH OF THE TWO GROUPS FOLLOWING A NORMAL DISTRIBUTION (PARAMETRIC)? 
# IF YES, USE T.TEST. IF NOT, USE WILCOXIN TEST
# SHAPIRO TEST: #p-val >0.05 implies data distribution not sig.dif. from a normal distribution.

#GROUP 1
df1 <- df %>% filter(strain == "C57")
ggdensity(df1$n)
ggqqplot(df1$n, 
         title = "C57 Females")
shapiro.test(df1$n)  

#GROUP 2
df1 <- df %>% filter(strain == "NYOB")
ggdensity(df1$n)
ggqqplot(df1$n, 
         title = "NYOB Females")
shapiro.test(df1$n) 

#### ASSUMPTION 3: DO THE TWO POPS HAVE THE SAME VARIANCE? USE THE F-TEST
# P-VALUE > 0.05, NO SIGNIFICIANT DIFFERENCE BETWEEN THE VARIANCES OF THE TWO DATA SETS.
var.test(n ~ strain, data = df)

# CHOOSE YOUR TEST. PARAMETRIC (T-TEST) OR NON-PARAMETRIC(WILCOX.TEST)
compare_means(n~strain, data = df, method = "wilcox.test", paired = FALSE)



# DIFFERENCE IN TOTAL # OF VISITS BY PADDOCK? --------------------------------------
df <- data %>% 
  group_by(Paddock) %>% 
  tally()
df

# PRODUCE SUMMARY STATS OF NUMBER OF RFID READS PER PADDOCK
df %>% 
  group_by(Paddock) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS


plot <- ggbarplot(df, 
                  x="Paddock", 
                  y="n", 
                  # color = "Trial",
                  fill = "lightblue",
                  # palette = "Paired",
                  ylab = "Total RFID Reads",
                  xlab = "Paddock",
                  # title = "Total RFID Reads Per Trial",
)
plot


# DIFFERENCE IN TOTAL READS BY zone AND PADDOCK? --------------------------

df <- data %>% 
  group_by(Paddock, antenna) %>% 
  tally()
df

# PRODUCE SUMMARY STATS OF NUMBER OF RFID READS PER PADDOCK
df %>% 
  group_by(Paddock, antenna) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS



plot <- ggbarplot(df, 
                  x="antenna",
                  y="n",
                  fill = "Paddock",
                  color = "Paddock",
                  palette = "Paired",
                  position = position_dodge(0.8),
                  # add = c("mean_sd", "jitter"),
                  ylab = "RFID Reads",
                  xlab = "antenna",
                  # title = "Total RFID Reads Per Trial")
)
plot



# TOTAL MALE READS PER zone PER TRIAL? ------------------------------------

df <- data %>% 
  group_by(Trial, sex, Full_ID, zone) %>% 
  filter(sex == "M") %>% 
  filter(Trial == "T007") %>% 
  tally()
df

# PRODUCE SUMMARY STATS 
df %>% 
  group_by(Full_ID, zone) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS



plot <- ggbarplot(df, 
                  x="antenna",
                  y="n",
                  facet.by = "Full_ID", 
                  title = "T007 Males"
)
plot



# TOTAL FEMALE READS PER zone PER TRIAL? ------------------------------------

df <- data %>% 
  group_by(Trial, sex, Full_ID, zone) %>% 
  filter(sex == "F") %>% 
  filter(Trial == "T007") %>% 
  tally()
df

# PRODUCE SUMMARY STATS 
df %>% 
  group_by(Full_ID, zone) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS



plot <- ggbarplot(df, 
                  x="zone",
                  y="n",
                  facet.by = "Full_ID", 
                  title = "T007 Females"
)
plot





# TOTAL RFID READS AND MALE AGD SCATTER PLOT? ----------------------------------------------------

df <- data %>% 
  group_by(sex, strain, Full_ID) %>% 
  tally()

meta <- metadata %>% select(strain,
                            sex,
                            ID,
                            full_ids,
                            dec.tag.id,
                            field.age,
                            pre.mass, 
                            post.mass, 
                            body.mm, 
                            full_body,
                            tail,
                            agd,
                            testes,
                            test_perc,
                            embryo.total)

meta$n <- df$n[match(meta$full_ids,df$Full_ID)]
match(meta$full_ids,df$Full_ID)
meta$agd <- as.numeric(meta$agd)
meta$test_perc <- as.numeric(meta$test_perc)

z<- meta %>% 
  group_by(strain,sex,agd) %>% 
  filter(sex == "M") %>% 
  filter(strain == "NYOB")

ggplot(z, aes(x=test_perc, y=n, color=strain)) + 
  geom_point(size=3) +
  scale_x_continuous("AGD",breaks=scales::pretty_breaks(n=10))
# stat_summary(fun.data=mean_cl_normal) + 
# geom_smooth(method='lm', formula= y~x)
# 
# 
