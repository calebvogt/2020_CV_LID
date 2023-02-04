
# Figs_Stats_Analysis_Code.R
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# Load packages -----------------------------------------------------------
library(readxl)
library(tidyverse)
library(data.table)
library(gganimate)
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
library(jtools)
library(sjPlot)
library(effects)
library(ggeffects)
library(ggrepel)
# library(rgl) # for 3D PCA plot. 

# Set directories and load metadata ---------------------------------------------------------
wd <- setwd("C:/Users/Caleb/Box/0_CV_Shared_Analyses/7_LID_2020/RFID_analysis_v10/1_manuscript")
output_fp <- paste("C:/Users/caleb/Desktop")
meta <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
meta_short <- meta %>% 
  select(trial, paddock, strain, sex, name, code, family_group)

# Load rfid_data and clean ----------------------------------------------------------
rfid_data <- as.data.frame(fread("LID_2020_ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE))

# Clean data to 10 days and triage relevant subjects as identified by descriptive analyses. 
rfid_data <- rfid_data %>%
  filter(noon_to_noon_day >=1 & noon_to_noon_day <= 10) %>% 
  #T004: George only mouse to cross between trials on Day 3. triage. 
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

rfid_data$strain_sex <- paste0(rfid_data$strain, "-", rfid_data$sex)
rfid_data$field_time <- as.POSIXct(rfid_data$field_time, format="%Y-%m-%d %H:%M:%OS")

# Load move_data and clean ----------------------------------------------------------
move_data <- as.data.frame(fread("LID_2020_ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE))
move_data <- subset(move_data, select = -c(V1))
move_data$field_time <- as.POSIXct(move_data$field_time, format="%Y-%m-%d %H:%M:%OS") # CONVERT TIME FORMAT
move_data$field_time_STOP <- as.POSIXct(move_data$field_time_STOP, format="%Y-%m-%d %H:%M:%OS")
move_data$strain_sex <- paste0(move_data$strain, "-", move_data$sex)

#Data Cleaning: removed George from all analyses due to movement between trials
move_data <- move_data %>% 
  filter(!(name == "George")) # george crosses from T004 to T005 on Day 2. swaps back and forth before settling as floater in T005 From Day 3 to day 10. Triage completely

# Load social_data -----------------------------------------------------------
filenames <- list.files(wd, pattern = "*MOVEBOUT_GBI.csv")
social_data = lapply(filenames, fread) ## READ IN ALL FILES


# move_data data summary statements -----------------------------------------------

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
  group_by(name, noon_to_noon_day) %>% 
  tally(sum(duration_hours)) 
sum(df$n)
summary(df$n)
mean(df$n)
std.error(df$n)

## Mean resource zones visited per mouse per night 
df <- move_data %>% # 
  group_by(name, noon_to_noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(name, noon_to_noon_day) %>%
  tally() # get number of distinct zone visits
mean(df$n)
std.error(df$n)

# Figure 1D and Figure SX-SX: Number of zone visits per night (heatmaps) --------
# GP heat map scale should be set to the max visits for any given category. Suri = 184 visits. 
df <- move_data %>%
  # filter(trial == "T001", sex == "M") %>%
  # filter(trial == "T001", sex == "F") %>%
  # filter(trial == "T002", sex == "M") %>%
  # filter(trial == "T002", sex == "F") %>%
  # filter(trial == "T003", sex == "M") %>%
  # filter(trial == "T003", sex == "F") %>%
  # filter(trial == "T004", sex == "M") %>%
  filter(trial == "T004", sex == "F") %>%
  # filter(trial == "T005", sex == "M") %>%
  # filter(trial == "T005", sex == "F") %>%
  # filter(trial == "T006", sex == "M") %>%
  # filter(trial == "T006", sex == "F") %>%
  # filter(trial == "T007", sex == "M") %>%
  # filter(trial == "T007", sex == "F") %>%
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
        } else {print("barnacles")}
        
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



# Figure 2A and Figure S2A-B: Unique zones visited per day -------------------------------------------------------------------
df <- move_data %>% # 
  group_by(trial, strain_sex, sex, strain, name, noon_to_noon_day, antenna) %>% 
  # filter(trial == "T007") %>% #Used for per trial graphs
  tally() %>%  # get number of visits to each zone
  group_by(trial, strain_sex, sex, strain, name, noon_to_noon_day, .drop = FALSE) %>%
  tally() %>% # get number of unique zone visits
  complete(noon_to_noon_day = 1:10, fill = list(n = 0)) %>% # fill in days where mouse doesnt appear with 0s
  dplyr::rename(unique_zones_visited = n)

# data cleaning. 
df <- df %>% 
  #T004: George only mouse to cross between trials on Day 3. triage. 
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

# output csv
# write.csv(df, paste(wd, "Figure2A_data.csv", sep = '/'), row.names = FALSE)

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



# Figure 2B: Cumulative unique zones visited over trial period --------
df <- move_data %>% 
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
# write.csv(df, paste(wd, "Figure2B_data.csv", sep = '/'), row.names = FALSE)

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


# Figure 2C: Estimated minimum distance traveled -------------------------------------------------------------------
df <- move_data %>% 
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
# write.csv(df7, paste(wd, "Figure2C_data.csv", sep = '/'), row.names = FALSE)

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

# Figure 2D: Percent time spent in most occupied zone -----------------------
df <- move_data %>% 
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
# write.csv(df1, paste(wd, "Figure2D_data.csv", sep = '/'), row.names = FALSE)

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

# Figure 2E-F: cumulative resource zone Priority Access Scores -----------------------

#get total time spent in a zone per day per trial for each MALE or FEMALE ONLY
# options(scipen=999)
df <- move_data %>% 
  filter(sex  == "M") %>%
  # filter(sex  == "F") %>%
  # filter(trial == "T001") %>%
  # filter(strain == "C57") %>%
  # filter(strain == "NYOB") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, noon_to_noon_day, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(trial_day_antenna = paste0(trial, "_",noon_to_noon_day, "_", antenna)) %>% 
  dplyr::rename(total_duration_min = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  select(trial_day_antenna, total_duration_min)

df1 <- move_data %>% 
  filter(sex  == "M") %>%
  # filter(sex  == "F") %>%
  # filter(trial == "T001") %>%
  # filter(strain == "C57") %>%
  # filter(strain == "NYOB") %>%
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
# write.csv(df8, paste0(output_fp, "/", "Figure2E-F_female_data.csv"), row.names = FALSE)

# line plot
(p <- ggplot(df8, aes(x=noon_to_noon_day, y=csum_daily_capture_penalty, group = name, color = name)) + #y=csum_adj_mus_percent_capture_score
    geom_line(size =1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(1,10,by=1), limits = c(1,10)) +
    scale_y_continuous(breaks = seq(-10,20, by = 5), limits = c(-10,20)) +
    xlab("Day") +
    ylab("Cumulative PA Score") +
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
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent") #main figs

# density plot
(p <- ggplot(subset(df8, noon_to_noon_day == 10), aes(x=csum_daily_capture_penalty, group = strain_sex, fill = strain_sex)) + 
    geom_density(adjust = 0.25,alpha = 0.4) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("red", "blue", "green", "darkgrey")) +
    xlab("Priority Access Score") +
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

# Figure 2G and Figure S2X: PCA of space use and movement patterns  --------------------------------------------------------------
## To generate PCA variable data from raw data, see Section 2

## Section 1: Create PCA graphs
df <- read_excel("Figs_data_v3.xlsx", sheet = "Fig2G", col_names = TRUE) # note, only includes mice that made it to day 10. 
df1 <- df %>% 
  column_to_rownames("name")

pca <- prcomp(df1, scale = TRUE) # always scale the data. 
plot(pca$x[,1], pca$x[,2])
biplot(pca)
scores <- pca$x
loadings <- pca$rotation
sqrt(1/ncol(df1)) # cutoff for significant loadings
write.table(loadings, "clipboard", sep="\t", row.names=FALSE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per[1:3]
sum(pca.var.per[1:3])

# PC scree plot
barplot(pca.var.per, 
        xlab = "PC", 
        ylab = "Percent Variation", 
        names.arg = 1:length(pca.var.per), las = 1)
abline(h=1/ncol(df)*100, col = "red")

# Quick Biplot
biplot(scores[, 1:2], 
       loadings[, 1:2], 
       col = c("red", "blue"),
       cex=c(0.8, 0.8),
       xlabs = rep("o", nrow(scores)),
       xlim = c(-10, 10))


## ggplot
loading_scores <- pca$rotation[,1]
move_scores <- abs(loading_scores)
move_scores_ranked <- sort(move_scores, decreasing = TRUE)
top_moves <- names(move_scores_ranked)
pca$rotation[top_moves,1]
df2 <- data.frame(name = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
df3 <- left_join(df2, meta) 
df4 <- df3 %>% 
  mutate(strain_sex = paste(strain, sex, sep= "-")) %>% 
  select(trial, strain_sex, strain,sex, name, X, Y)

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

## PCA Stats
# 1. Take PCA scores for each mouse and associate with other (per mouse) data
# 2. Then do something like the follow. Note that Comp.1 here is equivalent to PC1. 
tracking60s.hometub=read.csv("/Users/michaelsheehan/Downloads/tracking60s.hometub.csv", header=T)

pca.mod1=lmer(Comp.1~night+(1|mouse), data=tracking60s.hometub)
anova(pca.mod1)

pca.mod2=lmer(Comp.2~night+(1|mouse), data=tracking60s.hometub)
anova(pca.mod2)



# Section 2: Get all data that you want to put into the PCA. 
# total_zone_time_h
df <- move_data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(name) %>% 
  tally(sum(duration_hours)) 
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# mean_bout_visit_s
df <- move_data %>% 
  group_by(name) %>% 
  tally(mean(duration_s)) 
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_unique_zones_per_trial
df <- move_data %>% 
  group_by(name,antenna) %>% 
  tally() %>% 
  dplyr::count(name)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# mean_dist_zones_per_night
df <- move_data %>% # 
  group_by(name, noon_to_noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(name, noon_to_noon_day) %>%
  tally() %>% # get number of distinct zone visits
  select(name, n) %>% 
  summarize_all(funs(mean))
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# %_time_rank_1_zone
options(scipen=999)
df <- move_data %>% 
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
df <- move_data %>% 
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
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  arrange(desc(n)) %>% 
  group_by(name) %>% 
  slice(which.max(n)) %>% 
  select(name, n)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_time_rank_2_zone
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  arrange(desc(n)) %>% 
  group_by(name) %>% 
  slice(-which.max(n)) %>% 
  slice(which.max(n)) %>% 
  select(name, n)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

## Create df to get avg and min daily distance travelled. 
df <- move_data %>% 
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

# total_min_distance_traveled
df4 <- df3 %>% 
  group_by(strain_sex, name, noon_to_noon_day) %>% 
  tally(sum(dist_moved)) %>% 
  group_by(name) %>%
  tally(sum(n))
write.table(df4, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD 

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

# Figure 2X: Discriminant Analysis --------------------------------------------
library(MASS)
library(dplyr)
# library(rattle)
library(heplots) # for BoxM test
library(caret)
library(klaR)



df <- read_excel("Figs_data_v3.xlsx", sheet = "Fig2G", col_names = TRUE) # note, only includes mice that made it to day 10. 
df1 <- merge(df, meta_short, by = "name")
df2 <- df1 %>% 
  mutate(strain_sex = paste(strain, sex, sep= "-")) %>% 
  relocate(strain_sex) %>% 
  dplyr::select(-c(name, paddock, sex, strain, family_group, trial, code))

#check assumptions for LDA
boxM(df2[,2:12], df2$strain_sex) # technically, a significant value suggests that LDA assumptions are not met, and therefore quadratic discriminant analysis should be used... 

model <- qda(strain_sex ~ ., data = df2)
model.p <- predict(model)
# model accuracy
table(model.p$class, df2$strain_sex, dnn = c('Predicted Group', 'Actual Group'))
mean(model.p$class == df2$strain_sex) # 75% accuracy. 

# Partition plots
library(klaR)
partimat(df2[2:12], strain_sex, data = df2)

## get the x,y coordinates for the LDA plot
data.qda.values <- predict(model)$class
## create a dataframe that has all the info we need to draw a graph
plot.data <- data.frame(X=data.qda.values$x[,1], Y=data.lda.values$x[,2], strain_sex = df2$strain_sex)
head(plot.data)



# LDA model
model <- lda(strain_sex ~ ., data = df2)
model.p <- predict(model)$class
# determine how good are the model predictions versus true categorization
table(model.p, df2$strain_sex)
# cross validate the model, leave one out approach. 
model2 <- lda(strain_sex ~ ., data = df2, CV = TRUE)
# Look at the assigned classes for the observation
table(model2$class, df2$strain_sex) # performs worse. 

## get the x,y coordinates for the LDA plot
data.lda.values <- predict(model)$class
## create a dataframe that has all the info we need to draw a graph
plot.data <- data.frame(X=data.lda.values$x[,1], Y=data.lda.values$x[,2], strain_sex = df2$strain_sex)
head(plot.data)

## draw a graph using ggplot2
(p <- ggplot(data=plot.data, aes(x=X, y=Y, color = strain_sex, shape = strain_sex)) +
  geom_point(aes(shape = factor(strain_sex))) +
  scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                     values=c("red", "blue", "green", "darkgrey")) +
  theme_classic() +
  # xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  # ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
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









# Figure 3: Add Mike's code here ------------------------------------------

# Figure 4: add social network code here ----------------------------------

# Table S1: Numbers and Descritpive Stats  -------------------------------------
# Mean RFID reads per trial
library(plotrix)
df <- rfid_data %>% 
  group_by(trial) %>%  
  tally()
mean(df$n)
std.error(df$n)

# Mean mouse RFID reads per day
df <- rfid_data %>% 
  # filter(trial == "T007") %>% #per trial measures. omit for all trials. 
  group_by(name, noon_to_noon_day) %>%  
  tally()
mean(df$n)
std.error(df$n)


# Figure S1A: Daily inter-zone travel times ---------------------------------------
df <- rfid_data
df$zone <- df$antenna
ids <- unique(df$name)
day_list <- list()
data_list <- list()
aa = ids[4]
for(aa in ids[1:length(ids)]){
  df1 <- df %>% 
    filter(name == aa) 
  #get daily inter-travel times. 
  for(bb in 1:10) {
    df2 <- df1 %>% 
      filter(noon_to_noon_day == bb)
    
    # delete consecutive repeat antenna hits, only keep rows where antennas change. 
    df3 <- as.data.table(df2)[, .SD[1], by = rleid(df2$antenna)]
    day_list[[bb]] <- df3 %>% 
      select(trial, strain, sex, name, zone, field_time) %>% 
      mutate(diff = field_time - lag(field_time), 
             diff_secs = as.numeric(diff, units = 'secs'), 
             diff_mins = as.numeric(diff, units = 'mins'), 
             diff_hours = as.numeric(diff, units = 'hours'))
  }
  data_list[[aa]] <- do.call("rbind", day_list)
}
df4 <- do.call("rbind", data_list)
df5 <- na.omit(df4)
min(df5$diff_secs)
mean(df5$diff_mins)
max(df5$diff_hours)

# Graph
# sdat <- summary(df5$diff_mins)
# summStr <- paste(names(sdat), format(sdat, digits = 0), collapse = "; ")
# op <- par(mar = c(7,4,4,2) + 0.1)
hist(df5$diff_mins, 
     xlim = c(0, 1000),
     breaks = 10000,
     # main = stuff,
     main = "",
     xlab = "Inter-zone travel time (min)"
)
# title(sub = summStr, line = 5.5)
# par(op)
# export as svg. 


# Figure S1B: Within-zone inter-RFID interval and time window capture thresholds --------

df <- rfid_data
df$zone <- df$antenna
ids <- unique(df$name)
big_data_list <- list()
data_list <- list()
flag <- 1
aa = ids[1]
for(aa in ids[1:length(ids)]){
  print(paste("Processing mouse ",flag, " out of ", length(ids), sep=''))
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) 
  
  days <- unique(df1$noon_to_noon_day)
  day_list <- list()
  cc=days[1]
  for(cc in days[1:length(days)]){
    df2 <- df1 %>% 
      filter(noon_to_noon_day == cc)
    
    zones <- unique(df2$zone)
    zone_list <- list()
    bb = zones[1]
    for(bb in zones[1:length(zones)]){
      zone_list[[bb]] <- df2 %>% 
        filter(zone == bb) %>% 
        select(trial, strain, sex, name, noon_to_noon_day, zone, field_time) %>% 
        mutate(diff = field_time - lag(field_time), 
               diff_secs = as.numeric(diff, units = 'secs')) %>% 
        #remove 0s
        filter(diff_secs > 0)
    }
    #list of a mouses data per day per zone
    data_list [[aa]] <- do.call("rbind", zone_list)
  }
  #list of all mouse data per day per zone
  big_data_list[[cc]] <- do.call("rbind", data_list)
  flag <- flag + 1
  
}

# after loop finishes
df3 <- do.call("rbind", big_data_list)

#remove NAs
df3 <- df3[!is.na(df3$diff_secs),]
sdat <- summary(df3$diff_secs)
sdat
## get time interval where 95% of intervals below that value. 
sort(df3$diff_secs)[0.95*length(df3$diff_secs)]  # 10 days, all time, 95% = 11 seconds
sort(df3$diff_secs)[0.99*length(df3$diff_secs)]  # 10 days, all time, 99% = 153 seconds <<<< We select the most conservative estimate for all mice at all times

#base plot
options(scipen=999)

## Graph
# svg(file = paste0(output_fp, "/", "output.svg"))
summStr <- paste(names(sdat), format(sdat, digits = 2), collapse = "; ")
summStr
op <- par(mar = c(7,4,4,2) + 0.1)
hist(df3$diff_secs, 
     xlim = c(0, 200),
     # log = "y",
     breaks = 30000,
     main = "",
     xlab = "Within Tub inter-read interval (s)"
)
abline(v=c(11,153), col=c("red","blue"), lty=c(1,2), lwd=c(3,3))
title(sub = summStr, line = 5.5)
par(op)
dev.copy(svg, file = paste0(output_fp, "/", "output.svg"))
dev.off()

# Figure S1C: Number of zone visits and time in zone correlation -----------------------
# number of zone visits
df <- move_data %>% 
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(trial, strain, sex, name, antenna) %>% 
  tally()

# duration of time spent in zone
df2 <- move_data %>% 
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(trial, strain, sex, name, antenna) %>% 
  tally(sum(duration_min))

df3 <- as.data.frame(cbind(df, df2$n))
df3$strain_sex <- paste0(df3$strain,"-",df3$sex)
colnames(df3) <- c("trial", "strain", "sex", "name", "antenna", "num_visits", "total_time_min", "strain_sex")

library(ggpubr)
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

# Figure S1D: Circadian activity by strain and sex ------------------------
df <- rfid_data %>% 
  select(strain_sex, strain, sex, name, noon_to_noon_day, time) %>% 
  filter(strain == "C57")
# filter(strain == "NYOB")

class(df$time)
df$time <- as.POSIXct(df$time, format="%H:%M:%OS")
df$hours <- as.numeric(format(df$time, format="%H"))
df$strain_sex <- as.factor(df$strain_sex)

## GRAPH
(p <- ggplot(df, aes(x = hours, fill = strain_sex)) +
    geom_histogram(binwidth = 1) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("red", "blue", "green", "black")) +
    theme_classic() +
    xlab("Hours") +
    ylab("RFID Reads") +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          # legend.position = c(0.2,0.8))
          legend.position = "top")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=3, height=3) #main figs





# Figure S2C: Cumulative estimated time spent in resource zones --------
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

# Figure S2X: Average visit bout duration by sex and strain (in progress) ---------------------------












