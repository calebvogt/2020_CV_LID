# 4_analyze_ALLTRIAL_MOVEBOUT.R
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# Install Packages --------------------------------------------------------

# install.packages("emmeans")   ## estimated marginal means and p-vals
# install.packages("sjstats")   ## partial eta squared and cohens f effect size
# install.packages("lme4")      ## estimated multi level model  
# install.packages("lmerTest")  ## comprehensive anova output with p-vals
# install.packages("MuMIn")     ## R2 for multi model inference


# LOAD PACKAGES AND DATA --------------------------------------------------
library(tidyverse)
library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
library(gganimate)
library(ggpubr)
library(hrbrthemes)
library(gifski)
library(av)
library(data.table)
library(readr)
library(reshape)
library(lubridate)
library(scales)
library(plot.matrix)
library(transformr)
library(rstatix)
library(emmeans)
library(sjstats)
library(lme4)
library(lmerTest)
library(MuMIn)

wd <- setwd("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_analysis_v5")
# wd <- setwd("C:/Users/caleb/Box/CV_Shared_Analyses/7_LID_2020/RFID_analysis_v5")
output_fp <- paste("C:/Users/Caleb Vogt/Desktop")
metadata <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
weather <- read_excel("LID_2020_metadata.xlsx", sheet = 2, skip = 0)
weather$Weather_Time <- as.POSIXct(weather$Weather_Time, format = "%m/%d/%Y %H:%M") 
class(weather$Weather_Time)
data <- as.data.frame(fread("LID_2020_ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE))
data <- subset(data, select = -c(V1))
data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS") # CONVERT TIME FORMAT
data$Field_Time_STOP <- as.POSIXct(data$Field_Time_STOP, format="%Y-%m-%d %H:%M:%OS")

######### TESTING
## do some checks on the data
testing <- data %>% 
  # filter(Trial == "T005") %>% 
  # filter(Name == "Barron") %>%
  # filter(noon_to_noon_day == c(1:10))## this is throwing away literally 10's of thousands of rows of RFID data. I do not understand how. 
  filter(Name == "Barron")
# filter(noon_to_noon_day == c(1:10))## this is throwing away literally 10's of thousands of rows of RFID data. I do not understand how. 
# filter(noon_to_noon_day %in% (1:10))
# filter(noon_to_noon_day >=1 & noon_to_noon_day <= 10)
######### TESTING


# [PRISM] Number of unique zones visited per night, graph and stats -----------------------------------
df <- data %>% 
  # filter(Strain == "C57", Sex == "M") %>%
  # filter(Strain == "C57", Sex == "F") %>%
  # filter(Strain == "NYOB", Sex == "M") %>%
  filter(Strain == "NYOB", Sex == "F") %>%
  group_by(Name, noon_to_noon_day, Antenna) %>% 
  tally() %>% 
  group_by(Name, noon_to_noon_day) %>%
  tally()
# View(df)

#quick graph
# ggplot(df, aes(x=noon_to_noon_day, y = n)) +
#   # geom_line() +
#   geom_point() 
#   # geom_errorbar(aes(ymin=n-sd, ymax=n+sd))
  

# Prism Graph
ids <- unique(df$Name)
# CREATE EMPTY LISTS. 
idlist <- list()
for(aa in ids[1:length(ids)]){
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 10, ncol = 1))
  df1 <- subset(df, Name == aa)
  bb=1
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY ZONE
    if(nrow(df2) == 0){
      stats[bb,1] <- 0
    } else {
      stats[bb,1] <- paste(df2$n)
    }
    idlist[[aa]] <- stats
  }
}

master_class <- do.call("cbind",idlist) #THIS THROWS THE ERROR
# View(master_class)

# COPY THE OUTPUT TO THE CLIPBOARD & PASTE TO PRISM
write.table(master_class, "clipboard", sep="\t", row.names = FALSE, col.names =  FALSE)


## Run stats


x <- merge(df,metadata, by.x = "Name", by.y = "name")

meta$n <- df$n[match(meta$full_ids,df$Full_ID)]
match(meta$full_ids,df$Full_ID)
meta$agd <- as.numeric(meta$agd)
meta$test_perc <- as.numeric(meta$test_perc)




############### lme model
library(lme4)
library(lmerTest)


#model 1
reg_results = lmer(n ~ sex + strain + noon_to_noon_day + #fixed effects
                     sex*strain + sex*noon_to_noon_day + strain*noon_to_noon_day + sex*strain*noon_to_noon_day + #interaction effects
                     (1|trial) + (1|Name), data = x) #random effects

anova(reg_results)

#model 2: favored model
library(emmeans)
model_2 = lmer(n ~ sex + strain + noon_to_noon_day + #fixed effects
                 sex*strain + sex*noon_to_noon_day + strain*noon_to_noon_day + #interaction effects
                 (1|trial) + (1|Name), data = x) #random effects
anova(reg_results)

#overall no effect of strain or sex (Fsub1,292.11) =o.2354 , p = 0.63)
#pattern over time, strain by sex interaction, sex:noon_to_noon_day, strain:noon_to_noon_day

#post-hoc test to get effect which is driving the significance. 

#yields means of the different effects

m_model_2 = emmeans(model_2, ~ sex + strain + sex*strain)

contrast(m_model_2, 'tukey')
##Estimate means that the F C57s are visiting on average 1.9 more tubs than M C57s/night. 






#model 3
reg_results = lmer(n ~ sex + strain + #fixed effects
                     sex*strain + #interaction effects
                     (1|trial) + (1|Name) + (1|noon_to_noon_day), data = x) #random effects
anova(reg_results)


#model 4

reg_results = lmer(n ~ sex + strain + #fixed effects
                     sex*strain + #interaction effects
                     (1|trial) + (1|Name) + (1|noon_to_noon_day), data = x) #random effects
anova(reg_results)



#factor is a grouping variable
#integer is performing a regression. likely makes more sense if there is a thought that 






# [PRISM] OB NUMBER OF UNIQUE ZONES VISITED PER NIGHT -----------------------------------
# CHANGE FILTER TO GET DATA FOR EACH SEX AND PASTE INTO GRAPHPAD. 
df <- data %>% 
  filter(Strain == "NYOB", Sex == "M") %>% 
  group_by(Name, noon_to_noon_day, Antenna) %>% 
  tally() %>% 
  group_by(Name, noon_to_noon_day) %>%
  tally()

#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
ids <- unique(df$Name)

# CREATE EMPTY LISTS. 
idlist <- list()
for(aa in ids[1:length(ids)]){
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 10, ncol = 1))
  df1 <- subset(df, Name == aa)
  bb=1
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY ZONE
    if(nrow(df2) == 0){
      stats[bb,1] <- 0
    } else {
      stats[bb,1] <- paste(df2$n)
    }
    idlist[[aa]] <- stats
  }
}

master_class <- do.call("cbind",idlist) #THIS THROWS THE ERROR
View(master_class)

# COPY THE OUTPUT TO THE CLIPBOARD & PASTE TO PRISM
write.table(master_class, "clipboard", sep="\t", row.names = FALSE, col.names =  FALSE)








# [PRISM] Heatmaps for durations spent in zones per day, individual/all trials ------------------------------------------------------------
## SUM THE DURATIONS SPENT BY INDIVIDUAL PER DAY PER ZONE
df <- data %>% 
  filter(Trial == "T001", Sex == "M") %>%
  # filter(Trial == "T001", Sex == "F") %>%
  # filter(Trial == "T002", Sex == "M") %>%
  # filter(Trial == "T002", Sex == "F") %>%
  # filter(Trial == "T003", Sex == "M") %>%
  # filter(Trial == "T003", Sex == "F") %>%
  # filter(Trial == "T004", Sex == "M") %>%
  # filter(Trial == "T004", Sex == "F") %>%
  # filter(Trial == "T005", Sex == "M") %>%
  # filter(Trial == "T005", Sex == "F") %>%
  # filter(Trial == "T006", Sex == "M") %>%
  # filter(Trial == "T006", Sex == "F") %>%
# filter(Trial == "T007", Sex == "M") %>%
# filter(Trial == "T007", Sex == "F") %>%
mutate(duration_min = duration_s / 60) %>% 
  group_by(Name, noon_to_noon_day, Antenna) %>% 
  # tally(sum(duration_s))
  tally(sum(duration_min))


#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
ids <- unique(df$Name)

# CREATE EMPTY LISTS. 
idlist <- list()
daylist<- list()

# LOOP THROUGH EACH INDIVIDUAL AND PULL OUT NUMBER OF VISITS PER UNIQUE ZONE AND PUT INTO 2X4 GRID THAT LOCALIZES TO THE 
# ENCLOSURE SETUP. PUT EACH NIGHT OF ACTIVITY TO THE RIGHT FOR 10 NIGHTS. COPY AND PASTE DIRECLY INTO PRISM. 

flag=1
aa = ids[1]
for(aa in ids[1:length(ids)]){
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
  df1 <- subset(df, Name == aa)
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY ZONE
    if(nrow(df2) == 0){
      stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
    } else {
      
      for(cc in 1:nrow(df2)){
        if(df2$Antenna[cc] == 1){
          stats[4,1] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 2){
          stats[4,2] <-print(df2$n[cc])
        } else if(df2$Antenna[cc] == 3){
          stats[3,1] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 4){
          stats[3,2] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 5){
          stats[2,1] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 6){
          stats[2,2] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 7){
          stats[1,1] <- print(df2$n[cc])
        } else if (df2$Antenna[cc] == 8){
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
## IN THE FUTURE, CHANGE THE COLUMN NAMES
View(master_SASS)

# COPY THE OUTPUT TO THE CLIPBOARD  
write.table(master_SASS, "clipboard", sep="\t", row.names=FALSE)

# PRINT ORDER OF THE DATA SET, COPY INTO GRAPHPAD
ids
write.table(ids, "clipboard", sep="\t", row.names=FALSE)



# [PRISM] Heatmaps for frequency of visits to zones per day, individual/all trials ------------------------------------------------------------
# Heat map scale should be set to max 75 visits. Carter has most visits on a given day with 72. 

df <- data %>% 
  # filter(Trial == "T001", Sex == "M") %>%
  filter(Trial == "T001", Sex == "F") %>%
  # filter(Trial == "T002", Sex == "M") %>%
  # filter(Trial == "T002", Sex == "F") %>%
  # filter(Trial == "T003", Sex == "M") %>%
  # filter(Trial == "T003", Sex == "F") %>%
  # filter(Trial == "T004", Sex == "M") %>%
  # filter(Trial == "T004", Sex == "F") %>%
  # filter(Trial == "T005", Sex == "M") %>%
  # filter(Trial == "T005", Sex == "F") %>%
  # filter(Trial == "T006", Sex == "M") %>%
  # filter(Trial == "T006", Sex == "F") %>%
  # filter(Trial == "T007", Sex == "M") %>%
# filter(Trial == "T007", Sex == "F") %>%
group_by(Name, noon_to_noon_day, Antenna) %>% 
  tally() #number of visits


#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
ids <- unique(df$Name)

# NOTE THAT DATA FRAMED WILL BE IN ORDER OF THIS LIST
# View(ids)

# CREATE EMPTY LISTS. 
idlist <- list()
daylist<- list()

# LOOP THROUGH EACH INDIVIDUAL AND PULL OUT NUMBER OF VISITS PER UNIQUE ZONE AND PUT INTO 2X4 GRID THAT LOCALIZES TO THE 
# ENCLOSURE SETUP. PUT EACH NIGHT OF ACTIVITY TO THE RIGHT FOR 10 NIGHTS. COPY AND PASTE DIRECLY INTO PRISM. 

flag=1
# aa = ids[1]
for(aa in ids[1:length(ids)]){
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
  df1 <- subset(df, Name == aa)
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY ZONE
    if(nrow(df2) == 0){
      stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
    } else {
      
      for(cc in 1:nrow(df2)){
        if(df2$Antenna[cc] == 1){
          stats[4,1] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 2){
          stats[4,2] <-print(df2$n[cc])
        } else if(df2$Antenna[cc] == 3){
          stats[3,1] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 4){
          stats[3,2] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 5){
          stats[2,1] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 6){
          stats[2,2] <- print(df2$n[cc])
        } else if(df2$Antenna[cc] == 7){
          stats[1,1] <- print(df2$n[cc])
        } else if (df2$Antenna[cc] == 8){
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

# RBIND DATA FOR EACH INDIVIDUAL
master_SASS <- do.call("rbind",idlist)
# View(master_SASS)
head(master_SASS)

# PRINT ORDER OF THE DATA SET, COPY INTO GRAPHPAD
ids
write.table(ids, "clipboard", sep="\t", row.names=FALSE)

# COPY THE OUTPUT TO THE CLIPBOARD  
write.table(master_SASS, "clipboard", sep="\t", row.names=FALSE)



# [Figure 1] STRAIN AND SEX ZONE VISIT DURATION (MIN) VIOLIN plot --------------------------------------
library(ggplot2)

df <- data %>% 
  group_by(Strain, Sex, Name) %>% 
  mutate(duration_hours = duration_s / 3600) %>% 
  tally(sum(duration_hours))


# COPY OUTPUT TO CLIPBOARD. 
write.table(df, "clipboard", sep="\t", row.names=FALSE)



p <- ggplot(df) +
  aes(x = Strain, y = n, fill = Sex) +
  geom_violin(trim = TRUE) + 
  labs(x = "Strain", y = "Time in Zones (h)") +
  theme_classic()
p

## ADD STATS
p + stat_compare_means(method = "t.test")

ggsave("Strain by Sex Zone Visit Duration (min), Violin.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)





# [Figure 1, PRISM] Ranked zone duration (% time) -----------------------
## create table, copy and paste into Prism. 
df <- data %>% 
  filter(Strain == "C57", Sex == "M") %>%
  # filter(Strain == "C57", Sex == "F") %>%
  # filter(Strain == "NYOB", Sex == "M") %>%
  # filter(Strain == "NYOB", Sex == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally(sum(duration_min))

ids <- unique(df$Name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(Name == aa) %>% 
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



# [Figure 1, PRISM] Ranked zone duration (min) -----------------------
## create table, copy and paste into Prism. 
df <- data %>% 
  # filter(Strain == "C57", Sex == "M") %>%
  # filter(Strain == "C57", Sex == "F") %>%
  # filter(Strain == "NYOB", Sex == "M") %>%
  filter(Strain == "NYOB", Sex == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally(sum(duration_min))

ids <- unique(df$Name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(Name == aa) %>% 
    arrange(desc(n))
  
  df2 <- data.frame(df1[,5])
  df2[(nrow(df2)+1):8,] <- 0
  df3 <- df2[1:8,]
  
  data_list[[aa]] <- df3
}

df4 <- do.call("cbind", data_list)

# COPY OUTPUT TO CLIPBOARD. Port into excel, sort, and put into prism. 
write.table(df4, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE)




# [Figure 1, PRISM] Ranked zone visit (frequency) -----------------------
## create table, copy and paste into Prism. 
df <- data %>% 
  # filter(Strain == "C57", Sex == "M") %>%
  # filter(Strain == "C57", Sex == "F") %>%
  # filter(Strain == "NYOB", Sex == "M") %>%
  filter(Strain == "NYOB", Sex == "F") %>%
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally()

ids <- unique(df$Name)
data_list <- list()
aa = ids[2]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(Name == aa) %>% 
    arrange(desc(n))
  
  df2 <- data.frame(df1[,5])
  df2[(nrow(df2)+1):8,] <- 0
  df3 <- df2[1:8,]
  
  data_list[[aa]] <- df3
}

df4 <- do.call("cbind", data_list)

# COPY OUTPUT TO CLIPBOARD. Port into excel, sort, and put into prism. 
write.table(df4, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE)


# [Figure 1] Correlation coefficient zone visit frequency and time -----------------------
df <- data %>% 
  # filter(bout_status == "START") %>%
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally()

df2 <- data %>% 
  # filter(bout_status == "START") %>%
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally(sum(duration_min))

df3 <- as.data.frame(cbind(df, df2$n))
# df3 <- as.data.frame(cbind(df$n, df2$n))
# colnames(df3)
# cor(df3)
# plot(df3)
# df3

p <- ggscatter(df3, x = "n", y = "...6", 
               # color = "Strain",
               # color = "Sex",
               # palette = "jco",
               add = "reg.line",
               add.params = list(color = "blue", fill = "red"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               xlab = "RZ Visit Frequency",
               ylab = "Time (min)"
)

p + stat_cor(method = "pearson", p.accuracy = 0.001)
# p + stat_cor(aes(color = Sex), label.x = 2, method = "pearson", p.accuracy = 0.001)


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
    mutate(Name_Sex = paste(Name,Sex, sep= "-")) %>% 
    filter(Trial == i) %>%
    group_by(Name_Sex) %>% 
    tally(sum(duration_s / 60)) %>% 
    arrange(desc(n))
  
  p <- df %>%
    filter(!is.na(n)) %>%
    arrange(n) %>%
    mutate(Name_Sex=factor(Name_Sex, Name_Sex)) %>%
    ggplot(aes(Name_Sex, n) ) +
    geom_segment( aes(x = Name_Sex, xend = Name_Sex, y = 0, yend = n), color="grey") +
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
  
  ggsave(filename=paste(current_trial, this_strain, "Summed Mouse Zone Duration (min).png", sep = " "), 
         plot=p, width = 6.5, height = 5, dpi = 300, units = "in", device='png', path = output_fp)
  
  flag <- flag +1
}

p

# [PNG] TRIAL ZONE VISITS, VIOLIN ----------------------------------------------------
df <- data %>% 
  group_by(Trial, Strain, Name) %>%  
  tally()

p <- ggplot(df) +
  aes(x = Trial, y = n, fill = Strain) +
  geom_violin() +
  labs(title = "Trial Zone Visits",
       subtitle = "",
       x = "Trial",
       y = "Zone Visits",
       caption = "") +
  theme_classic() 
p

p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

ggsave("Trial Zone Visits.png", 
       # plot=p, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)


# [PNG] TRIAL ZONE DURATION (MIN), VIOLIN ----------------------------------------------------
df <- data %>% 
  group_by(Trial, Strain, Name) %>%  
  mutate(duration_min = duration_s / 60) %>% 
  tally(sum(duration_min))

p <- ggplot(df) +
  aes(x = Trial, y = n, fill = Strain) +
  geom_violin() +
  labs(title = "Trial Zone Duration (min)",
       subtitle = "",
       x = "Trial",
       y = "Time (min)",
       caption = "") +
  theme_classic() 
p

p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

ggsave("Trial Zone Duration (Min), Violin.png", 
       # plot=p, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)





# [PNG] STRAIN BY SEX ZONE VISIT FREQUENCY, VIOLIN --------------------------------------
library(ggplot2)

df <- data %>% 
  group_by(Strain, Sex, Name) %>% 
  tally()

p <- ggplot(df) +
  aes(x = Strain, y = n, fill = Sex) +
  geom_violin(trim = TRUE) + 
  labs(title = "Zone Visits",
       subtitle = "",
       x = "Strain",
       y = "Number of Zone Visits",
       caption = "") +
  theme_classic() 
p

## ADD STATS
p + stat_compare_means(method = "t.test")

ggsave("Strain by Sex Zone Visit Frequency, Violin.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)




# keep adding things hereee.... -------------------------------------------





# CREATE DF OF AVERAGE DURATION OF VISITS PER ZONE PER MALE ---------------

# GET TOTAL DURATION OF ALL VISITS ACROSS ALL TRIALS. 
sum(data$duration_s)

means <- data %>%
  filter(Trial == "T001", Sex == "M") %>% 
  group_by(Sex, name, noon_to_noon_day, Zone) %>%
  dplyr::summarize(Mean = mean(duration_s, na.rm=TRUE)) #average duration of visit. 




# Q: IS THERE A DIFFERENCE IN THE NUMBER OF MALE AND FEMALE VISITS? --------------------------------------------

df <- data %>% 
  group_by(Sex, Full_ID) %>% 
  tally()

# COPY OUTPUT TO CLIPBOARD. 
write.table(df, "clipboard", sep="\t", row.names=FALSE)

# PRODUCE SUMMARY STATS
df %>% 
  group_by(Sex) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS

# NEXT, QUICKLY PLOT THE DATA TO VISUALIZE
plot <- ggboxplot(df, 
                  x="Sex", 
                  y="n", 
                  color = "Sex",
                  # palette = c("#00AFBB", "#E7B800"),
                  palette = "jco",
                  add = "jitter",
                  ylab = "RFID Reads",
                  xlab = "Sex",
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
df1 <- df %>% filter(Sex == "M")
ggdensity(df1$n)
ggqqplot(df1$n)
shapiro.test(df1$n) #p-val >0.05 implies data distribution not sig.dif. from a normal distribution. 

# GROUP 2
df1 <- df %>% filter(Sex == "F")
ggdensity(df1$n)
ggqqplot(df1$n)
shapiro.test(df1$n) 

#### ASSUMPTION 3: DO THE TWO POPS HAVE THE SAME VARIANCE? USE THE F-TEST
# P-VALUE > 0.05, NO SIGNIFICIANT DIFFERENCE BETWEEN THE VARIANCES OF THE TWO DATA SETS.
var.test(n ~ Sex, data = df)

# CHOOSE YOUR TEST. PARAMETRIC (T-TEST) OR NON-PARAMETRIC(WILCOX.TEST)
compare_means(n~Sex, data = df, method = "wilcox.test", paired = FALSE)




# PLOT TOTAL # OF VISITS VISITS BY MALES AND FEMALES, MF ----------------------------------------------------

# SUBSAMPLE DATA
df <- data %>% 
  group_by(Trial,Sex) %>% 
  tally()
df

plot <- ggbarplot(df, 
                  x="Trial",
                  y="n",
                  fill = "Sex",
                  color = "Sex",
                  palette = "Paired",
                  position = position_dodge(0.8),
                  ylab = "# OF ZONE VISITS",
                  xlab = "Trial",
                  title = "# OF ZONE VISITS")
plot



# DIFFERENCE IN MF RFID READS BY TRIAL? --------
df <- data %>% 
  group_by(Trial, Sex, Full_ID) %>% 
  tally()
df

plot <- ggbarplot(df, 
                  x="Trial",
                  y="n",
                  fill = "Sex",
                  color = "Sex",
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
  group_by(Trial,Strain, Sex, Full_ID) %>% 
  filter(Sex == "M") %>% 
  tally()


# PRODUCE SUMMARY STATS
df %>% 
  group_by(Strain) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS

# NEXT, QUICKLY PLOT THE DATA TO VISUALIZE
plot <- ggboxplot(df, 
                  x="Strain", 
                  y="n", 
                  color = "Strain",
                  palette = "jco",
                  add = "jitter",
                  ylab = "# OF ZONE VISITS",
                  xlab = "Sex",
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
df1 <- df %>% filter(Strain == "C57")
ggdensity(df1$n)
ggqqplot(df1$n, 
         title = "C57 Males")
shapiro.test(df1$n)  

#GROUP 2
df1 <- df %>% filter(Strain == "NYOB")
ggdensity(df1$n)
ggqqplot(df1$n,
         title = "NYOB Males")
shapiro.test(df1$n) 

#### ASSUMPTION 3: DO THE TWO POPS HAVE THE SAME VARIANCE? USE THE F-TEST
# P-VALUE > 0.05, NO SIGNIFICIANT DIFFERENCE BETWEEN THE VARIANCES OF THE TWO DATA SETS.
var.test(n ~ Strain, data = df)

# CHOOSE YOUR TEST. PARAMETRIC (T-TEST) OR NON-PARAMETRIC(WILCOX.TEST)
compare_means(n~Strain, data = df, method = "wilcox.test", paired = FALSE)

# DIFFERENCE IN # OF ZONE VISITS FOR NYOB AND C57 FEMALES? ----------------------

df <- data %>% 
  group_by(Trial,Strain, Sex, Full_ID) %>% 
  filter(Sex == "F") %>% 
  tally()


# PRODUCE SUMMARY STATS
df %>% 
  group_by(Strain) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS

# NEXT, QUICKLY PLOT THE DATA TO VISUALIZE
plot <- ggboxplot(df, 
                  x="Strain", 
                  y="n", 
                  color = "Strain",
                  palette = "jco",
                  add = "jitter",
                  ylab = "# OF VISITS",
                  xlab = "Strain",
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
df1 <- df %>% filter(Strain == "C57")
ggdensity(df1$n)
ggqqplot(df1$n, 
         title = "C57 Females")
shapiro.test(df1$n)  

#GROUP 2
df1 <- df %>% filter(Strain == "NYOB")
ggdensity(df1$n)
ggqqplot(df1$n, 
         title = "NYOB Females")
shapiro.test(df1$n) 

#### ASSUMPTION 3: DO THE TWO POPS HAVE THE SAME VARIANCE? USE THE F-TEST
# P-VALUE > 0.05, NO SIGNIFICIANT DIFFERENCE BETWEEN THE VARIANCES OF THE TWO DATA SETS.
var.test(n ~ Strain, data = df)

# CHOOSE YOUR TEST. PARAMETRIC (T-TEST) OR NON-PARAMETRIC(WILCOX.TEST)
compare_means(n~Strain, data = df, method = "wilcox.test", paired = FALSE)



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


# DIFFERENCE IN TOTAL READS BY ZONE AND PADDOCK? --------------------------

df <- data %>% 
  group_by(Paddock, Antenna) %>% 
  tally()
df

# PRODUCE SUMMARY STATS OF NUMBER OF RFID READS PER PADDOCK
df %>% 
  group_by(Paddock, Antenna) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS



plot <- ggbarplot(df, 
                  x="Antenna",
                  y="n",
                  fill = "Paddock",
                  color = "Paddock",
                  palette = "Paired",
                  position = position_dodge(0.8),
                  # add = c("mean_sd", "jitter"),
                  ylab = "RFID Reads",
                  xlab = "Antenna",
                  # title = "Total RFID Reads Per Trial")
)
plot



# TOTAL MALE READS PER ZONE PER TRIAL? ------------------------------------

df <- data %>% 
  group_by(Trial, Sex, Full_ID, Zone) %>% 
  filter(Sex == "M") %>% 
  filter(Trial == "T007") %>% 
  tally()
df

# PRODUCE SUMMARY STATS 
df %>% 
  group_by(Full_ID, Zone) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS



plot <- ggbarplot(df, 
                  x="Antenna",
                  y="n",
                  facet.by = "Full_ID", 
                  title = "T007 Males"
)
plot



# TOTAL FEMALE READS PER ZONE PER TRIAL? ------------------------------------

df <- data %>% 
  group_by(Trial, Sex, Full_ID, Zone) %>% 
  filter(Sex == "F") %>% 
  filter(Trial == "T007") %>% 
  tally()
df

# PRODUCE SUMMARY STATS 
df %>% 
  group_by(Full_ID, Zone) %>%  #what are the groups
  get_summary_stats(n, type = "mean_sd")  #PRINT SUMMARY STATISTICS



plot <- ggbarplot(df, 
                  x="Zone",
                  y="n",
                  facet.by = "Full_ID", 
                  title = "T007 Females"
)
plot





# TOTAL RFID READS AND MALE AGD SCATTER PLOT? ----------------------------------------------------

df <- data %>% 
  group_by(Sex, Strain, Full_ID) %>% 
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
