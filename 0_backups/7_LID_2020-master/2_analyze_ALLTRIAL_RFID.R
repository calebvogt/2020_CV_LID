# 2_analyze_ALLTRIAL_RFID
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# LOAD PACKAGES & DATA -----------------------------------------------------------
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

wd <- setwd("C:/Users/Caleb/Box/0_CV_Shared_Analyses/7_LID_2020/RFID_analysis_v10")
output_fp <- paste("C:/Users/caleb/Desktop")
meta <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
# weather <- read_excel("LID_2020_metadata.xlsx", sheet = 2, skip = 0)
# weather$Weather_Time <- as.POSIXct(weather$Weather_Time, format = "%m/%d/%Y %H:%M")
data <- as.data.frame(fread("LID_2020_ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE))
data$strain_sex <- paste0(data$strain, "-", data$sex)
## CLEAN DATA DOWN TO 10 DAYS + remove george from the analysis
data <- data %>%
  filter(noon_to_noon_day >=1 & noon_to_noon_day <= 10) 
#removed george from all analyses due to movement between trials
data <- data %>% 
  filter(name != "George")
data$field_time <- as.POSIXct(data$field_time, format="%Y-%m-%d %H:%M:%OS")


# Summary RFID Statistics ---------------------------------------------
# mean rfid reads per trial and coefficient of variation
df <- data %>% 
  group_by(trial) %>%  
  tally()
summary(df$n)
mean(df$n)
std <- function(x) sd(x)/sqrt(length(x))
std(df$n)
cv <- sd(df$n) / mean(df$n) * 100 #coeffcient of variation
cv


# mean mouse reads per night 
df <- data %>% 
  group_by(name, noon_to_noon_day) %>%  
  tally()
summary(df$n)
mean(df$n)
std(df$n)
cv <- sd(df$n) / mean(df$n) * 100 #coeffcient of variation
cv


# FigS1A_Within-zone RFID read interval and time window capture thresholds -----------------------------------------------
df <- data
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
# pdf(file = paste0(output_fp, "/", "output.pdf"))
summStr <- paste(names(sdat), format(sdat, digits = 4), collapse = "; ")
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

dev.copy(pdf, file = paste0(output_fp, "/", "output.pdf"))
dev.off()


library(pastecs)
stat.desc(df3$diff_secs)
summary(df3$diff_secs)

## ggplot of log frequency counts of various within tub-read intervals
df4 <- df3 %>% 
  group_by(diff_secs) %>% 
  tally()
options(scipen=999)

ggplot(df4) + 
  geom_histogram(aes(x=diff_secs, y =..density.., weight = log(n)))


# FigS2D_Circadian activity by strain and sex ------------------------------------
df <- data %>% 
  select(strain_sex, strain, sex, name, noon_to_noon_day, time) %>% 
  filter(strain == "C57")
  # filter(strain == "NYOB")

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

class(df$time)
# df$time <- as.POSIXct(df$time, format="%H:%M:%OS")

df$time <- strptime(df$time, format="%H:%M:%S")
df$hours <-  as.numeric(format(df$time, format="%H"))
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





# [Figure S3] Between tub individual travel time from RFID hits ---------------------------

df$Zone <- df$Antenna

ids <- unique(df$Name)
data_list <- list()
aa = ids[4]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(Name == aa) 
  
  # delete consecutive repeat antenna hits, only keep rows where antennas change. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$Antenna)]
  
  data_list[[aa]] <- df2 %>% 
    select(Trial, Strain, Sex, Name, Zone, field_time) %>% 
    # mutate(diff = field_time - lag(field_time), diff_secs = as.numeric(diff, units = 'secs'))
    mutate(diff = field_time - lag(field_time), diff_mins = as.numeric(diff, units = 'mins'))
}

df4 <- do.call("rbind", data_list)


write.csv(df4, file = "temp.csv")


df5 <- df4 %>% 
  filter(Strain == "C57", Sex == "M")
# filter(Strain == "C57", Sex == "F")
# filter(Strain == "NYOB", Sex == "M")
# filter(Strain == "NYOB", Sex == "F")

stuff <- paste(unique(df5$Strain), unique(df5$Sex))
sdat <- summary(df4$diff_mins)
summStr <- paste(names(sdat), format(sdat, digits = 4), collapse = "; ")
op <- par(mar = c(7,4,4,2) + 0.1)
# hist(df5$diff_secs)
hist(df4$diff_mins, 
     xlim = c(0, 60),
     breaks = 10000,
     main = stuff,
     # main = "",
     xlab = "Between tub travel time (min)"
)
title(sub = summStr, line = 5.5)
par(op)




# [R, Violin] # of RFID reads per mouse per trial -------------------------
df <- data %>% 
  group_by(Trial, Strain, Name) %>%  
  tally()

# plot
ggplot(df) +
  aes(x = Trial, y = n, fill = Strain) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  labs(x = "Trial",  y = "Total RFID Reads Per Mouse") +
  theme_test() 
#resize in window
ggsave("output.png", dpi = 300,device='png',path = output_fp)

## stats
summary(df) 
df$Trial <- as.factor(df$Trial)
df$Strain <- as.factor(df$Strain)
df$Name <- as.factor(df$Name)

#anova
one.way <- aov(n ~ Trial, data = df)
summary(one.way)

mean(df$n)



##
model = lmer(n ~ sex + strain + noon_to_noon_day + 
               sex*strain + sex*noon_to_noon_day + strain*noon_to_noon_day +
               (1|trial) + (1|Name), data = stats) 
summary(model)
anova(model)
eta_sq(model, partial = TRUE) # partial eta sq
r.squaredGLMM(model) # adjust R2 for the model as an alternative

##Post-hoc for main effects
emmeans(model, pairwise ~ noon_to_noon_day) #throws error

# post-hoc for significant interactions
emmeans(model, pairwise ~ sex | strain)
emmeans(model, pairwise ~ sex | noon_to_noon_day) #noon_to_noon day is an integer, so its performing a regression
emmeans(model, pairwise ~ strain | noon_to_noon_day)


#quick assumptions check
plot(model)

qqnorm(resid(model))
hist(resid(model))


# [R, Violin, DONE] # of RFIDs per mouse by sex and strain --------------------------------------
df <- data %>% 
  group_by(Strain, Sex, Name) %>% 
  tally()

## stats

## plot
p <- ggplot(df) +
  aes(x = Sex, y = n, fill = Sex) +
  geom_violin(trim = TRUE) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  labs(x = "Sex",
       y = "Total RFID Reads Per Mouse",
       caption = "") +
  theme_test() +
  facet_wrap(~ Strain)
p

## ADD STATS
# p + stat_compare_means(method = "anova")

#size in rstudio window
ggsave("plot.png", 
       dpi = 300, 
       device='png', 
       path = output_fp)


# GIF: 12 HOUR INDIVIDUAL PADDOCK ACTIVITY, ALL MICE -----------------------
library(gganimate)
library(gifski)
library(av)
df <- clean
ids <- unique(df$Name)
ids
i=ids[1]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  
  # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  move_df <- subset(df, Name == current_mouse)
  
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(move_df$Trial), unique(move_df$Strain), unique(move_df$Sex), sep = " ")
  
  ## SET TIME VARIABLE
  
  p3 <- ggplot(move_df, aes(zone_x, zone_y)) +
    ggtitle(paste(current_mouse_info, current_mouse, "Movement", sep = " ")) +
    geom_point(show.legend = FALSE, alpha = 0.7, size = 2) +
    xlim("A","B") +
    ylim("A","B","C","D") +
    geom_jitter(width = 0.1, height = 0.1) +
    scale_color_viridis_d() +
    labs(x = "none", y = "none") +
    theme(plot.background = element_blank(),
          panel.grid.major = element_line(colour="black", size = 0.5),
          panel.grid.minor = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 20))
  
  # RENDER THE ANIMATION
  p3a <- p3 +
    # geom_line(aes(group=seq_along(Field.Time))) + #remove?
    transition_time(field_time) +
    
    ## FEATURE: REPLACE DATE WITH NIGHT # + HOURLY TIME TRANSITION
    labs(subtitle = "Time: {frame_time}")
  
  
  plot(p3)
  # SAVE THE ANIMATION
  animate(p3a, duration = 60, fps = 24, width = 300, height = 300, renderer = gifski_renderer())
  anim_save(paste(current_mouse_info, current_mouse, "Paddock Activity.gif", sep ="_"), path = output_fp)
}




# PNG: 12 HOUR MOUSE ACTIVITY, ALL MICE LOOP -----------------------
df <- clean
# df <- subset(clean, Trial == "T001")
ids <- unique(df$Name)
ids
i=ids[1]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  
  # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  move_df <- subset(df, Name == current_mouse)
  
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(move_df$Trial), unique(move_df$Strain), unique(move_df$Sex), sep = " ")
  
  ## PNG: STATIC INDIVIDUAL MOVEMENT ACROSS TUBS FOR ENTIRE TIME PERIOD.
  p <- ggplot(move_df) +
    aes(x = field_time, y = Antenna) +
    geom_point(na.rm=TRUE, size=1, color = "black") +
    ggtitle(paste(current_mouse_info, current_mouse, "Activity", sep = " ")) +
    xlab("Date") + 
    ylab("Zone") +
    scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
  
  ggsave(filename=paste(current_mouse_info, current_mouse, "Activity.png", sep = " "), 
         plot=p, 
         width = 5, 
         height = 4, 
         dpi = 300, 
         units = "in", 
         device='png', 
         path = output_fp)
} 


# PNG: 12 HOUR TRIAL READS, VIOLIN ----------------------------------------------------
df1 <- data %>% 
  group_by(Trial, Strain, Name) %>%  
  tally()

p <- ggplot(df1) +
  aes(x = Trial, y = n, fill = Strain) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  labs(x = "Trial",
       y = "RFID Reads per mouse",
       caption = "") +
  theme_classic() 
p

ggsave("12h Trial Reads.png", 
       # plot=p, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)


# PNG: 12 HOUR C57 NIGHTLY READS, BAR -------------------------------------
library(ggsignif)

## CREATE DATA FRAME
df <- data %>% 
  filter(Strain == "C57") %>% 
  group_by(Sex, Name, noon_to_noon_day) %>%  
  tally()

## GET SUMMARY STATISTICS
df.summary <- df %>% 
  group_by(noon_to_noon_day, Sex) %>% 
  summarise(mean_grp = mean(n), # MEAN
            sd_grp = sd(n, na.rm = TRUE), # STANDARD DEVIATION
            n_grp = n(), # SAMPLE SIZE
            se_grp = sd(n)/sqrt(n())) ## STANDARD ERROR

## PLOT SUMMARY STATISTICS DATAFRAME WITH SEM ERROR BARS
p <- ggplot(df.summary) +
  aes(x = noon_to_noon_day, y = mean_grp, fill = Sex) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "C57 RFID Reads Per Night",
       subtitle = "",
       x = "Night",
       y = "RFID Reads",
       caption = "Error bars indicate s.e.m.") +
  scale_x_discrete(limits = c(1:10)) +
  geom_errorbar(aes(ymin = mean_grp - se_grp, ymax = mean_grp + se_grp),
                width=.2,
                position=position_dodge(.9)) +
  theme_classic()
p

ggsave("12h C57 nightly reads.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)

# PNG: 12 HOUR NYOB NIGHTLY READS, BAR -------------------------------------
## CREATE DATA FRAME
df <- clean %>% 
  filter(Strain == "NYOB") %>% 
  group_by(Sex, Name, noon_to_noon_day) %>%  
  tally()

## GET SUMMARY STATISTICS
df.summary <- df %>% 
  group_by(noon_to_noon_day, Sex) %>% 
  summarise(mean_grp = mean(n), # MEAN
            sd_grp = sd(n, na.rm = TRUE), # STANDARD DEVIATION
            n_grp = n(), # SAMPLE SIZE
            se_grp = sd(n)/sqrt(n())) ## STANDARD ERROR

## PLOT SUMMARY STATISTICS DATAFRAME WITH SEM ERROR BARS
p <- ggplot(df.summary) +
  aes(x = noon_to_noon_day, y = mean_grp, fill = Sex) +
  geom_bar(position = "dodge", stat = "identity") +
  labs(title = "NYOB RFID Reads Per Night",
       subtitle = "",
       x = "Night",
       y = "RFID Reads",
       caption = "Error bars indicate s.e.m.") +
  scale_x_discrete(limits = c(1:10)) +
  geom_errorbar(aes(ymin = mean_grp - se_grp, ymax = mean_grp + se_grp),
                width=.2,
                position=position_dodge(.9)) +
  theme_classic()
p

ggsave("12h NYOB nightly reads.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)




# GIF: 12 HOUR ALL TRIAL MALE PADDOCK ACTIVITY -----------------------
trial <- unique(clean$Trial)

## FILTER BY SEX, 20 GRIDS IS TOO MUCH. 
df <- clean %>% 
  filter(Sex == "M")

## LOOP ACROSS TRIALS
i=trial[4]
for(i in trial[1:length(trial)]) {
  current_trial <- print(i)
  
  move_df <- subset(df, Trial == current_trial)
  
  sex <- unique(df$Sex)
  
  # PLOT 4: ALL INDIVIDUAL MOVMENT ANIMATION ANIMATION
  p4 <- ggplot(move_df, aes(zone_x, zone_y, color = Name)) +
    ggtitle("Group Movement") +
    geom_point(show.legend = TRUE, alpha = 0.7, size = 2) +
    xlim("A","B") +
    ylim("A","B","C","D") +
    geom_jitter(width = 0.1, height = 0.1) +
    scale_color_viridis_d() +
    labs(x = "none", y = "none") +
    facet_wrap(~Name) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_line(colour="black", size = 0.5),
          panel.grid.minor = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 20))
  
  # RENDER THE ANIMATION
  p4a <- p4 +
    transition_time(field_time) +
    labs(subtitle = "Time: {frame_time}")
  
  # SAVE THE ANIMATION
  animate(p4a, duration = 60, fps = 24, width = 500, height = 500, renderer = gifski_renderer())
  anim_save(paste(current_trial,sex, "Group_Movement_FULL.gif", sep ="_"), path = output_fp)
}


# GIF: 12 HOUR MOUSE ACTIVITY, ALL MICE  -----------------------
df <- clean
ids <- unique(df$Name)
ids
i=ids[1]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  
  # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  move_df <- subset(df, Name == current_mouse)
  
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(move_df$Trial), unique(move_df$Strain), unique(move_df$Sex), sep = " ")
  
  p2 <- ggplot(move_df) +
    aes(x = field_time, 
        y = Antenna, 
        color = factor(Name)) +
    geom_line(na.rm=TRUE, color="red", size=1) +
    ggtitle(paste(current_mouse_info, current_mouse, "Movement", sep = " ")) +
    scale_color_viridis_d() +
    labs(x = "Date", y = "Zone") +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
  
  # RENDER THE ANIMATION
  p2a <- p2 +
    geom_point(aes(group=seq_along(field_time))) +
    transition_reveal(field_time)
  
  # SAVE THE ANIMATION
  animate(p2a, duration = 10, fps = 20, width = 300, height = 300, renderer = gifski_renderer())
  anim_save(paste(current_mouse_info, current_mouse, "Activity.gif", sep =" "), path = output_fp)
  
}




# ?? STUFF FOR sna NETWORKS, WHAT DO I DO WITH THIS??? ------------------------


# Network Measures
# GET GRAPH CENTRALITY MEASURES AT NODE AND NETWORK LEVEL
degree.cent <- centr_degree(net_graph, mode = "all")

#NODE LEVEL CENTRALITY MEASURES
degree.cent$res

# NETWORK GRAPH LEVEL CENTRALITY MEASURE (COMPARE ACROSS DAYS?)
degree.cent$centralization
degree(g_undir, mode='all')
degree(g_undir, mode='in')

# FOR NEW, CONVERT TIMES COLUMN FROM SECONDS BACK INTO FIELD TIME. 
new1<- data %>% 
  filter(Trial == "T007") %>% 
  filter(Sex == "M") %>% 
  select(Time_sec,
         field_time)

test <- merge(new, new1, by.x = "times", by.y = "Time_sec")
test <- test %>% 
  relocate(times, .after = field_time)
test$sum_rows <- rowSums(test[,1:10]) #NUMBER OF COLUMNS CHANGES EACH TIME... SO ADJUST OR FIND A WORKAROUND. 

test1<- unique(test)

test2 <- test1 %>% 
  filter(sum_rows > 1)


test <- df$n[match(meta$full_ids,df$Full_ID)]
match(meta$full_ids,df$Full_ID)
# test <- new %>% 
#   mutate(sum_rows = sum(new[,1:10]))


# WORKING WITH GMM RDATA FILES 
# READ IN GMM.RDATA FILES. 

# EXTRACT OUTPUT FROM INDIVIDUAL GMM FILES. 
gbi <- gmm_data$gbi
events <- gmm_data$metadata
observations_per_event <- gmm_data$B

# Split up location and date data
tmp <- strsplit(events$Location,"_")
tmp <- do.call("rbind",tmp)
events$Date <- tmp[,1]
events$Location <- tmp[,2]

#Turning into network
r_network <- get_network(association_data = gmm_data$gbi, data_format = "GBI")

r_network_ALL <- get_network(association_data = gmm_data$gbi, data_format = "GBI")
r_net_ALL <- graph.adjacency(r_network_ALL, mode = "undirected", weighted = TRUE, diag = FALSE)

#Making network MF
r_network_MF <- r_network_ALL
r_network_MF[which(global_ids =="*-M-*"),which(global_ids =="*-M-*")] <- 0
r_network_MF[which(global_ids =="*-F-*"),which(global_ids =="*-F-*")] <- 0
save(r_network, file = "All RFID Network")
save(r_network_MF, file = "MF RFID Network")
write.csv(r_network_MF, file = "MF RFID Network.csv")

#Unweighted
r_network_MF_uw <- r_network_MF
r_network_MF_uw[r_network_MF_uw > 0] <- 1
degree_rfidMF <- rowSums(r_network_MF_uw)
degree_rfidMF
save(r_network_MF_uw, file = "MF RFID Network UW") #?? as csv???


# [in prog] bad plotRFID READS BY INDIVIDUAL, ZONE, ALL TRIAL PLOTS -----------------------------------
df <- data %>% 
  group_by(Trial, Sex, Name, Antenna) %>% 
  filter(Trial == "T001", Sex == "M") %>%    # CHANGE
  tally()

ggbarplot(df, 
          x="Antenna",
          y="n",
          facet.by = c("Name"), ncol = 5,
          ylim = c(0, 500000),
          xlab = "Zone",
          ylab = "RFID Reads",
          title = "T001")  # CHANGE



