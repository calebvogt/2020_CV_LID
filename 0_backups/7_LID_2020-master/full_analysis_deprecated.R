# 4_full_analysis.R
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# LOAD PACKAGES -----------------------------------------------------------
## Load packages
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


# [#########SECTION 1#########] ALLTRIAL_RFID ANALYSES & STATS --------------------------------------------------
#set wd and output file path
wd <- setwd("G:/My Drive/PROJECTS/7_LID_2020/RFID_analysis_v8")
output_fp <- paste("C:/Users/caleb/Desktop")

# load data
meta <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
# weather <- read_excel("LID_2020_metadata.xlsx", sheet = 2, skip = 0)
# weather$Weather_Time <- as.POSIXct(weather$Weather_Time, format = "%m/%d/%Y %H:%M") 

# load all trial RFID data
data <- as.data.frame(fread("LID_2020_ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE))
data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS")

## CLEAN DATA DOWN TO 10 DAYS
data <- data %>%
  filter(noon_to_noon_day >=1 & noon_to_noon_day <= 10) 

## quick stats
# average reads per mouse nights ()
df <- data %>% 
  # filter(Trial == "T007" & Sex == "M") %>% 
  group_by(Trial) %>%  
  tally()
summary(df$n)
mean(df$n)

# [R, S1A, DONE] Within-tub Read Interval & Time window capture thresholds, histogram -----------------------------------------------
## use this to determine your MOVEBOUT time threshold value. Critical! Takes inter-read intervals within a tub for each mouse for each day
# and removes all 0 values (minimum difference is 1 second, as we remove miliseconds from our data) and finds 95% and 99% capture threshold
df <- data
df$Zone <- df$Antenna
ids <- unique(df$Name)
big_data_list <- list()
data_list <- list()
flag <- 1
aa = ids[1]
for(aa in ids[1:length(ids)]){
  print(paste("Processing mouse ",flag, " out of ", length(ids), sep=''))
  ## create df of percent time
  df1 <- df %>% 
    filter(Name == aa) 
  
  days <- unique(df1$noon_to_noon_day)
  day_list <- list()
  cc=days[1]
  for(cc in days[1:length(days)]){
    df2 <- df1 %>% 
      filter(noon_to_noon_day == cc)
    
    zones <- unique(df2$Zone)
    zone_list <- list()
    bb = zones[1]
    for(bb in zones[1:length(zones)]){
      zone_list[[bb]] <- df2 %>% 
        filter(Zone == bb) %>% 
        select(Trial, Strain, Sex, Name,noon_to_noon_day, Zone, Field_Time) %>% 
        mutate(diff = Field_Time - lag(Field_Time), 
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



# [R, FigS2] 24hour Strain circadian activity (RFID reads per hour) ------------------------------------

## ISSUE 1: REMOVE INDIVIDUALS KNOWN TO BE NESTING UNDER THE RFID READERS. ESPECIALLY TRUE FOR OBS? 

df <- data %>% 
  select(Strain, Sex, noon_to_noon_day, Time) %>% 
  filter(Strain == "C57")
# filter(Strain == "NYOB")
df$Time <- as.POSIXct(df$Time, format="%H:%M:%OS")

# base r plot
pdf(file = paste0(output_fp, "/", "output.pdf"))
hist(df$Time, 
     breaks = "hours", 
     col = "gray", 
     xlab = "Hours",
     ylab = "RFID Reads",
     # ylim = c(0,25000),
     main = "C57",
     # main = "OB",
     freq = TRUE)
dev.off()


## ggplot attempt
ggplot(df, aes(x = Time)) +
  geom_histogram() +
  xlab("Hours") +
  ylab("RFID Reads") +
  stat_bin(bins = 24) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
# stat_count(width = 0.5)


scale_x_continuous(breaks = seq(1,10, by =1), limits = c(1,10)) +
  scale_y_continuous(breaks = seq(1,8, by =1), limits =c(1,8)) +
  
  
  
  
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
    select(Trial, Strain, Sex, Name, Zone, Field_Time) %>% 
    # mutate(diff = Field_Time - lag(Field_Time), diff_secs = as.numeric(diff, units = 'secs'))
    mutate(diff = Field_Time - lag(Field_Time), diff_mins = as.numeric(diff, units = 'mins'))
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
    transition_time(Field_Time) +
    
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
    aes(x = Field_Time, y = Antenna) +
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
    transition_time(Field_Time) +
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
    aes(x = Field_Time, 
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
    geom_point(aes(group=seq_along(Field_Time))) +
    transition_reveal(Field_Time)
  
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
         Field_Time)

test <- merge(new, new1, by.x = "times", by.y = "Time_sec")
test <- test %>% 
  relocate(times, .after = Field_Time)
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









# [#########SECTION 2#########] ALLTRIAL_MOVEBOUT ANALYSES AND STATS ------------------------------------

# Load Data 
wd <- setwd("G:/My Drive/PROJECTS/7_LID_2020/RFID_analysis_v8")
output_fp <- paste("C:/Users/caleb/Desktop")

#load data
meta <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
# weather <- read_excel("LID_2020_metadata.xlsx", sheet = 2, skip = 0)
# weather$Weather_Time <- as.POSIXct(weather$Weather_Time, format = "%m/%d/%Y %H:%M") 

# Load movebout data: rows correspond to a per mouse visit bout to specific zones, with estimated start and stop times. includes other relevant metadata information. 
data <- as.data.frame(fread("LID_2020_ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE))
data <- subset(data, select = -c(V1))
data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS") # CONVERT TIME FORMAT
data$Field_Time_STOP <- as.POSIXct(data$Field_Time_STOP, format="%Y-%m-%d %H:%M:%OS")

# Movebout quick data exploration -----------------------------------------------


## Explore the data
df <- data %>% 
  # filter(Strain == "C57") %>%
  # filter(Strain == "NYOB") %>%
  # filter(Sex == "M") %>% 
  # filter(Sex == "F") %>% 
  # filter(Trial == "T001") %>% 
  # filter(Trial == "T002") %>% 
  # filter(Trial == "T003") %>% 
  # filter(Trial == "T004") %>% 
  # filter(Trial == "T005") %>% 
  # filter(Trial == "T006") %>% 
  # filter(Trial == "T007") %>% 
mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(Trial, Strain, Sex) %>% 
  # tally() #number of visits
  tally(sum(duration_hours)) 

sum(df$n)


#quickly graph the quantiative variables
qplot(Antenna, n, data = df, colour = Sex) +
  stat_summary(fun.data = mean_cl_normal) +
  geom_smooth(method="lm")


# [GP, FigS1, DONE] Number of Zone Visits per night heatmap ---------------------------
# GP heat map scale should be set to the max visits for any given category. Suri = 184 visits. 
df <- data %>%
  # filter(Trial == "T001", Sex == "M") %>%
  # filter(Trial == "T001", Sex == "F") %>%
  # filter(Trial == "T002", Sex == "M") %>%
  # filter(Trial == "T002", Sex == "F") %>%
  # filter(Trial == "T003", Sex == "M") %>%
  # filter(Trial == "T003", Sex == "F") %>%
  # filter(Trial == "T004", Sex == "M") %>%
  filter(Trial == "T004", Sex == "F") %>%
  # filter(Trial == "T005", Sex == "M") %>%
  # filter(Trial == "T005", Sex == "F") %>%
  # filter(Trial == "T006", Sex == "M") %>%
  # filter(Trial == "T006", Sex == "F") %>%
  # filter(Trial == "T007", Sex == "M") %>%
  # filter(Trial == "T007", Sex == "F") %>%
  group_by(Name, noon_to_noon_day, Antenna) %>%
  tally() #number of visits

#check maximum visit number! adjust GP heatmap settings accordingly. 
# Suri = 123 visits

ids <- unique(df$Name)#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
idlist <- list() # NOTE THAT DATA FRAMED WILL BE IN ORDER OF THIS LIST
daylist<- list()
flag=1
# aa = ids[1]
for(aa in ids[1:length(ids)]){ # LOOP THROUGH EACH INDIVIDUAL AND PULL OUT NUMBER OF VISITS PER UNIQUE ZONE AND PUT INTO 2X4 GRID THAT LOCALIZES TO THE 
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 4, ncol = 2)) # ENCLOSURE SETUP. PUT EACH NIGHT OF ACTIVITY TO THE RIGHT FOR 10 NIGHTS. COPY AND PASTE DIRECLY INTO PRISM. 
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
master_SASS <- do.call("rbind",idlist) # RBIND DATA FOR EACH INDIVIDUAL
# View(master_SASS)
head(master_SASS)
write.table(master_SASS, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  
ids # PRINT ORDER OF THE DATA SET, COPY INTO GRAPHPAD

# [GP, DONE] Time in zones per night heatmap ------------------------------------------------------------
## SUM THE DURATIONS SPENT BY INDIVIDUAL PER DAY PER ZONE
# Havent used these because the range of durations on the first day screws up the scale. 
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


# [R, FigS1, DONE] Correlation coefficient zone visit frequency and time, scatterplot -----------------------
df <- data %>% 
  # filter(bout_status == "START") %>%
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally()
# df$n =  Number of VISITS to a zone over the entire trial

df2 <- data %>% 
  # filter(bout_status == "START") %>%
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally(sum(duration_min))
#df$n = total time spent in a particular zone per mouse

df3 <- as.data.frame(cbind(df, df2$n))
df3$strain_sex <- paste0(df3$Strain,"-",df3$Sex)

p <- ggscatter(df3, x = "n", y = "...6", 
               # color = "strain_sex",
               color = "Sex",
               # color = "Strain", 
               # palette = "jco",
               add = "reg.line",
               add.params = list(color = "blue", fill = "red"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               xlab = "Number of zone visits",
               ylab = "Time spent in zone (min)"
)

p + stat_cor(method = "pearson", p.accuracy = 0.001)
# p + stat_cor(aes(color = Sex), label.x = 2, method = "pearson", p.accuracy = 0.001)

ggsave("ouput.pdf", plot = p, device='pdf', path = output_fp)

# plot 2
df3 %>% 
  ggplot(aes(x = n, 
             y = ...6, 
             color = strain_sex)) +
  geom_point() +
  scale_color_brewer(palette = "PuOr") + 
  geom_smooth(method="lm") + 
  xlab("Number of zone visits") +
  ylab("Time spent in zone (min)") +
  stat_cor(aes(color = strain_sex), label.x = 4, method = "pearson", p.accuracy = 0.001) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

#save
ggsave("output.pdf",device='pdf', path = output_fp)


# [GP, Fig2, DONE] Distinct zones visited per night, line graph--------
######### STATS
library(emmeans)
library(sjstats)
library(lme4)
library(lmerTest)
library(MuMIn)


options(scipen = 0) # turn of scientific notation

# Effect of sex/strain/time on unique number of zone visits? 
df <- data %>% # 
  group_by(Name, noon_to_noon_day, Antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(Name, noon_to_noon_day) %>%
  tally() # get number of distinct zone visits

stats <- merge(df,meta, by.x = "Name", by.y = "name")

stats <- stats %>% 
  select(trial, Name, sex, strain, noon_to_noon_day, n) # where n is the number of unique zone visits per night

stats$trial <- as.factor(stats$trial)
stats$sex <- as.factor(stats$sex)
stats$strain <- as.factor(stats$strain)
summary(stats) #noon_to_noon day left as numeric instead of factor. regression? when there is a categorical variable it 

# model. Stick with model one, despite AIC, simpler with pairwise interaction effects. 
model1 = lmer(n ~ sex + strain + noon_to_noon_day + # fixed effects 
               sex*strain + sex*noon_to_noon_day + strain*noon_to_noon_day + # interaction effects
               (1|trial) + (1|Name), data = stats) # random effects

# model assumptions check
plot(model1)
qqnorm(resid(model1))
hist(resid(model1))

# model stats
# summary(model) ### throws different p vals for the sex effect... unclear why. 
anova(model1)
# eta_sq(model, partial = TRUE) # partial eta sq. not reported
# r.squaredGLMM(model) # adjust R2 for the model as an alternative. not reported. 


model1emm <- emmeans(model1, pairwise ~ sex | noon_to_noon_day) #noon_to_noon day is an integer, so its performing a regression
emmeans(model, pairwise ~ noon_to_noon_day | strain)
m.model1=emmeans(model, ~ sex*strain)
contrast(m.model1, 'tukey', simple = 'each')


# 
# 
# 
# 
# 
# 
# 
# 
# model2 = lmer(n ~ sex*strain*noon_to_noon_day + (1|trial) + (1|Name), data = stats) # random effects
# # model assumptions check
# plot(model2)
# qqnorm(resid(model2))
# hist(resid(model2))
# 
# # model stats
# # summary(model) ### throws different p vals for the sex effect... unclear why. 
# anova(model2)
# eta_sq(model, partial = TRUE) # partial eta sq
# r.squaredGLMM(model) # adjust R2 for the model as an alternative
# 
# #check performanace of each model. 
# AIC(model1)
# AIC(model2) #smaller values in this case means a "better" model
# 
# # post-hoc for significant interactions
# model2emm <- emmeans(model2, ~ sex*strain*noon_to_noon_day)
# contrast(model2emm, method = "tukey")
# 
# 
# 
# 
# 
# emmeans(model2, pairwise ~ sex | strain, adjust = "tukey")
# 
# 


#### GRAPH 

df <- data %>% 
  # filter(Strain == "C57", Sex == "M") %>%
  # filter(Strain == "C57", Sex == "F") %>%
  # filter(Strain == "NYOB", Sex == "M") %>%
  # filter(Strain == "NYOB", Sex == "F") %>%
  group_by(Name, noon_to_noon_day, Antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(Name, noon_to_noon_day) %>%
  tally() # get number of unique zone visits per night

ids <- unique(df$Name)
idlist <- list()
for(aa in ids[1:length(ids)]){
  stats <- data.frame(matrix(0, nrow = 10, ncol = 1)) # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  df1 <- subset(df, Name == aa)
  bb=1
  for(bb in 1:10){
    df2 <- subset(df1, noon_to_noon_day == bb)
    if(nrow(df2) == 0){  # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY ZONE
      stats[bb,1] <- 0
    } else {
      stats[bb,1] <- paste(df2$n)
    }
    idlist[[aa]] <- stats
  }
}
master_class <- do.call("cbind",idlist)
write.table(master_class, "clipboard", sep="\t", row.names = FALSE, col.names =  FALSE)



# [GP, Fig2, DONE] Ranked zone use (% time), line graph -----------------------

options(scipen=999)
df <- data %>% 
  # filter(Strain == "C57", Sex == "M") %>%
  # filter(Strain == "C57", Sex == "F") %>%
  # filter(Strain == "NYOB", Sex == "M") %>%
  # filter(Strain == "NYOB", Sex == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(Strain, Sex, Name, Antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time))

# model. Stick with model one, despite AIC, simpler with pairwise interaction effects. 
model1 = lmer(n ~ sex + strain + noon_to_noon_day + # fixed effects 
                sex*strain + sex*noon_to_noon_day + strain*noon_to_noon_day + # interaction effects
                (1|trial) + (1|Name), data = stats) # random effects

# model assumptions check
plot(model1)
qqnorm(resid(model1))
hist(resid(model1))

# model stats
# summary(model) ### throws different p vals for the sex effect... unclear why. 
anova(model1)





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


# [GP, line] Percent time in top-ranked zone per night -----------------------
# https://stackoverflow.com/questions/45341541/group-by-in-dplyr-and-calculating-percentages

options(scipen=999)
df <- data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(Strain, Sex, Name, noon_to_noon_day, Antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time))


dplyr::count(Name, noon_to_noon_day)

group_by(Name, noon_to_noon_day, n) %>% 
  transmute(noon_to_noon_day, Percentage = )
mutate(daily_sum = sum(Name, noon_to_noon_day, n))





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


# [R, Fig2, DONE] Cumulative unique zones visited per night, line graph --------

######### STATS
library(emmeans)
library(sjstats)
library(lme4)
library(lmerTest)
library(MuMIn)

df <- data %>% 
  select(Trial, Strain, Sex, Name, noon_to_noon_day, Zone) %>%
  group_by(Trial, Strain, Sex, Name, Zone) %>% 
  distinct(Zone, .keep_all = TRUE) %>% # get distinct zones ever visited, keep associated day it was visited
  group_by(Trial, Strain, Sex, Name, noon_to_noon_day) %>% #regroup by day
  tally() %>% #tally unique zones visited per day
  mutate(csum_novel_zones = cumsum(n)) %>%  #get cumulative # of novel mice met
  complete(Name, noon_to_noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(Name, noon_to_noon_day) %>% 
  fill(csum_novel_zones) %>% ## fill cumulative sum data from last observed day
  mutate(group = paste0(Strain, "-", Sex)) %>% 
  select(Trial, group, Strain, Sex, Name, noon_to_noon_day, csum_novel_zones)


## model
# model. Stick with model one, despite AIC, simpler with pairwise interaction effects. 
model1 = lmer(csum_novel_zones ~ Sex + Strain + noon_to_noon_day + # fixed effects 
                Sex*Strain + Sex*noon_to_noon_day + Strain*noon_to_noon_day + # interaction effects
                (1|Trial) + (1|Name), data = df) # random effects

# model assumptions check
plot(model1)
qqnorm(resid(model1))
hist(resid(model1))
anova(model1)

m.model1=emmeans(model1, ~ Sex*Strain)
contrast(m.model1, 'tukey', simple = 'each')

#C57 stand out once again. 


# write.csv(df, file = "zone_accumulation_df.csv")
p <- ggplot(df, aes(x=noon_to_noon_day, y=csum_novel_zones, colour = group)) + 
  geom_smooth() +
  xlab("Night") +
  ylab("Cumulative unique zones visited") +
  scale_x_continuous(breaks = seq(1,10, by =1), limits = c(1,10)) +
  scale_y_continuous(breaks = seq(1,8, by =1), limits =c(1,8)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p

ggsave("ouput.pdf", plot = p, device='pdf', path = output_fp)

# [R, FigS2, DONE] Cumulative time spent in zones over trial, line graph --------
df <- data %>% 
  select(Strain, Sex, Name, noon_to_noon_day,duration_s) %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(Strain, Sex, Name, noon_to_noon_day) %>% 
  tally(sum(duration_min)) %>% 
  mutate(csum_zone_min = cumsum(n)) %>%  #get cumulative # of novel mice met
  complete(Name, noon_to_noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(Name, noon_to_noon_day) %>% 
  fill(csum_zone_min) %>% ## fill cumulative sum data from last observed day
  mutate(group = paste0(Strain, "-", Sex)) %>% 
  select(group, Strain, Sex, Name, noon_to_noon_day, csum_zone_min)

# write.csv(df, file = "zone_accumulation_df.csv")
p <- ggplot(df, aes(x=noon_to_noon_day, y=csum_zone_min, colour = group)) + 
  geom_smooth() +
  xlab("Night") +
  ylab("Cumulative time spent in zones (min)") +
  scale_x_continuous(breaks = seq(1,10, by = 1), limits = c(1,10)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
p

ggsave("ouput.pdf", plot = p, device='pdf', path = output_fp)

# [R, PCA, IN PROGRESS] PCA of zone usage???  -----------------------------
# https://www.youtube.com/watch?v=0Jp4gsfOLMs


# [GP, line, IN PROGRESS] % Time in Top Ranked Zone Per Night -----------------------
options(scipen=999)

df <- data %>% 
  select(Strain, Sex, Name, noon_to_noon_day, Zone) %>%
  group_by(Strain, Sex, Name, Zone) %>% 
  distinct(Zone, .keep_all = TRUE) %>% # get distinct zones ever visited, keep associated day it was visited
  group_by(Strain, Sex, Name, noon_to_noon_day) %>% #regroup by day
  tally() %>% #tally unique zones visited per day
  mutate(csum_novel_zones = cumsum(n)) %>%  #get cumulative # of novel mice met
  complete(Name, noon_to_noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(Name, noon_to_noon_day) %>% 
  fill(csum_novel_zones) %>% ## fill cumulative sum data from last observed day
  mutate(group = paste0(Sex, "-", Strain)) %>% 
  select(group, Strain, Sex, Name, noon_to_noon_day, csum_novel_zones)

write.csv(df, file = "zone_accumulation_df.csv")
p <- ggplot(df, aes(x=noon_to_noon_day, y=csum_novel_zones, colour = group)) + 
  geom_smooth()
p

ggsave("ouput.png", plot = p, device='png', path = output_fp)


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
  # filter(Strain == "C57", Sex == "M") %>%
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





# [GP, line, NOT USED] Ranked zone visits (freq) -----------------------
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



# [R, ] STRAIN AND SEX ZONE VISIT DURATION (MIN) VIOLIN plot --------------------------------------
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

# [#########SECTION 3#########] TRIAL_MOVEBOUT_GBI ANALYSIS AND STATS -----------------------------------

#set fp
# wd <- setwd("C:/Users/Caleb Vogt/Box/0_CV_Shared_Analyses/7_LID_2020/RFID_analysis_v6")
wd <- setwd("C:/Users/caleb/Box/0_CV_Shared_Analyses/7_LID_2020/RFID_analysis_v6")
# output_fp <- paste("C:/Users/Caleb Vogt/Desktop")
output_fp <- paste("C:/Users/caleb/Desktop")

#load data
metadata <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
meta_short <- metadata %>% 
  select(trial, paddock, strain, sex, name, code, family_group)

## READ IN DATA
filenames <- list.files(wd, pattern = "*MOVEBOUT_GBI.csv")

## READ IN ALL FILES
myfiles = lapply(filenames, fread)

# [R, Fig4, DONE] M, F, and M+F SNA plots per day per trial -----------------------------------------

# choose data
# data <- as.data.frame(fread("T001_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
# data <- as.data.frame(fread("T002_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
data <- as.data.frame(fread("T003_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
# data <- as.data.frame(fread("T004_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
# data <- as.data.frame(fread("T005_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
# data <- as.data.frame(fread("T006_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
# data <- as.data.frame(fread("T007_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
data <- subset(data, select = -c(V1))

i=1
for(i in 1:10) { # loop through days
  gbi <- data %>% 
    filter(Day == i) %>%
    select(matches("*C57*"))
  # select(matches("*NYOB*"))
  # select(matches("*-M-"))
  # select(matches("*-F-"))
  
  #remove strain-sex info from colnames  
  colnames(gbi)<-gsub("C57-M-","",colnames(gbi))
  colnames(gbi)<-gsub("C57-F-","",colnames(gbi))
  # colnames(gbi)<-gsub("NYOB-M-","",colnames(gbi))
  # colnames(gbi)<-gsub("NYOB-F-","",colnames(gbi))
  
  ids <- colnames(gbi)
  undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION, # TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
                              data_format = "GBI",
                              #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
                              association_index = "SRI")
  
  net <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
  net <- simplify(net) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
  
  # set node attributes
  V(net)$size=degree(net)*2 #because 1 is a small size for a node, I'm just multiplying it by 5
  V(net)$name
  V(net)$sex = as.character(meta_short$sex[match(V(net)$name,meta_short$name)])
  V(net)$color = V(net)$sex #assign the "Sex" attribute as the vertex color
  V(net)$color = gsub("F","red",V(net)$color) #Females will be red
  V(net)$color = gsub("M","blue",V(net)$color) #Males will be blue
  
  ## M+F plot
  pdf(file=paste0(output_fp,"/","Night ", i, ".pdf"))
  par(mar = c(0.4, 0.1, 2, 0.1))
  plot.igraph(net,
              vertex.label=NA,
              layout = layout.fruchterman.reingold
              # layout=layout_nicely
  )
  title(paste0("Night ", i), cex.main=3)
  dev.off()
  
  
  # 
  # 
  # ### plot 2
  # # x <- get.adjacency(net)
  # # View(x)
  # # plot(net)
  # graph.strength(net) #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
  # V(net)$label <- V(net)$name
  # V(net)$degree <- degree(net)
  # 
  # # CREATE YOUR PLOT
  # png(file=paste0(output_fp,"/","Night ", i, ".png"))
  # par(mar = c(0.4, 0.1, 2, 0.1))
  # plot.igraph(net,
  #             vertex.color = "lightblue", #change
  #             # vertex.color = "red",
  #             vertex.size = 50, #20
  #             # vertex.size = igraph::degree(net)*5, #SET NODE SIZE AS A FUNCTION OF DEGREE CENTRALITY MULTIPLIED BY A SCALAR
  #             vertex.label.color = "black",
  #             vertex.label.font = 4, #changes font type
  #             vertex.label.cex = 1.5, #0.75
  #             edge.width = E(net)$weight*100, #maintain original weights
  #             edge.color = 'black',
  #             edge.curved = 0.5,
  #             layout = layout_in_circle(net, order = ids) # SORT ALPHABETICALLY FOR REPEATED GRAPHS ACROSS DAYS
  # )
  # title(paste0("Night ", i), cex.main=3)
  # dev.off()    # SAVE THE PNG FILE
  # 
  
}






# [R, Fig4, DONE] Cumulative sum of unique mice encountered, line ------------
## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
males <- meta_short %>% 
  filter(sex == "M") %>% 
  select(name) %>% 
  filter(!is.na(name))
male_list <- dplyr::pull(males, name)

## female names list
females <- meta_short %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  select(name) %>% 
  filter(!is.na(name))
female_list <- dplyr::pull(females, name)

trial_stats <- list()
aa = 1
for(aa in 1:length(myfiles)){
  ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- myfiles[[aa]]
  df <- subset(df, select = -c(V1))
  
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum)
  
  colnames(df2)<-gsub("C57-M-","",colnames(df2))
  colnames(df2)<-gsub("C57-F-","",colnames(df2))
  colnames(df2)<-gsub("NYOB-M-","",colnames(df2))
  colnames(df2)<-gsub("NYOB-F-","",colnames(df2))
  
  ## get mouse column names starting at col 10
  col_ids <- colnames(df2[,10:ncol(df2)])
  bb = col_ids[1]
  all_mouse_list <- list()
  first_flag = 1
  for(bb in col_ids[1:length(col_ids)]) {
    df3 <- df2 %>% 
      filter((!!as.symbol(bb)) == 1) %>% 
      mutate(Name = bb) %>% 
      relocate(Name)
    
    # remove current mouse from the next loop to compare to other animals
    non_self_ids <- col_ids[!col_ids %in% bb]
    non_self_ids[1]
    novel_mouse_rows <- list()
    second_flag = 1
    for(i in non_self_ids[1:length(non_self_ids)]) {
      df4 <- df3 %>% 
        filter((!!as.symbol(i)) == 1) %>% 
        mutate(novel_mouse_met = i)
      #save first observed meeting of the focal and novel mouse to list
      novel_mouse_rows[[second_flag]] <- df4[1,]
      second_flag = second_flag + 1
    }
    all_mouse_list[[first_flag]] <- do.call("rbind", novel_mouse_rows)
    first_flag <- first_flag + 1
  }
  df5 <- do.call("rbind", all_mouse_list)
  
  #remove na rows which are introduced when a mouse does not ever meet a particular other mouse. 
  df5[rowSums(is.na(df5)) > 0,]
  df6 <- df5[complete.cases(df5), ]
  
  ## ADD RELEVANT METADATA INFORMATION. 
  df7 <- merge(df6, meta_short, by.x = "Name", by.y = "name")
  df8 <- df7 %>%
    select(trial, paddock, strain, sex, Name, code, family_group, novel_mouse_met, Day, Zone, Field_Time_Start, Field_Time_Stop, m_sum, f_sum, mf_sum, duration_s) %>% 
    relocate(trial, paddock, strain, sex, Name, code, family_group, novel_mouse_met, Day, Zone,  Field_Time_Start, Field_Time_Stop, m_sum, f_sum, mf_sum, duration_s)
  
  trial_stats[[aa]] <- df8
}
df9 <- do.call("rbind", trial_stats)
df10 <- df9[with(df9, order(Name, Day)),]

# df10 has a dataframe ordered by first meeting time with each mouse in the paddock .

df11 <- df10 %>% 
  mutate(group = paste0(strain, "-", sex)) %>% 
  # filter(group == "M-NYOB") %>%
  # filter(sex == "M") %>% 
  group_by(group, Name, Day) %>% 
  tally() %>% 
  mutate(csum = cumsum(n)) %>% #get cumulative # of novel mice met
  complete(Name, Day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(Name, Day) %>% 
  fill(csum) ## fill cumulative sum data from last observed day

# write.csv(df11, file = "output.csv")

p <- ggplot(df11, aes(x=Day, y=csum, colour = group)) + 
  scale_x_continuous(breaks = seq(1,10, by = 1)) +
  scale_y_continuous(breaks = seq(1,20, by = 1)) +
  geom_smooth() +
  xlab("Night") +
  ylab("Cumulative unique mice encountered") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

p

ggsave("output.pdf", plot = p, device='pdf', path = output_fp)



# [R, SNA] All mouse, all days SNA plots per trial -----------------------------------------
#yank data
data <- myfiles[[4]]

gbi <- data %>% 
  filter(Day == 10) %>% ## adjust trials
  # select(matches("*C57*"))
  select(matches("*NYOB*")) 

colnames(gbi)<-gsub("C57-M-","",colnames(gbi))
colnames(gbi)<-gsub("C57-F-","",colnames(gbi))
colnames(gbi)<-gsub("NYOB-M-","",colnames(gbi))
colnames(gbi)<-gsub("NYOB-F-","",colnames(gbi))

# CREATE IDS LIST FROM COLUMN NAMES
ids <- colnames(gbi)
undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION
                            data_format = "GBI",
                            #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
                            association_index = "SRI")

# CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
net <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE)
net <- simplify(net) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 

# set node attributes
V(net)$name
V(net)$sex = as.character(meta_short$sex[match(V(net)$name,meta_short$name)])
V(net)$color = V(net)$sex #assign the "Sex" attribute as the vertex color
V(net)$color = gsub("F","red",V(net)$color) #Females will be red
V(net)$color = gsub("M","blue",V(net)$color) #Males will be blue

png(file=paste0(output_fp,"/","output.png")) #OPEN PNG FOR SAVING
plot.igraph(net,
            vertex.label=NA,
            layout=layout.fruchterman.reingold)
dev.off()




#measures
centr_degree(net, mode = "all") # GET GRAPH CENTRALITY MEASURES AT NODE AND NETWORK LEVEL# RES = NODE CENTRALITY, #CENTRALIZATION = GRAPH CENTRALITY, THEORETICAL MAX = ALL POSSIBLE EDGES. 
degree(net, mode="all") #VERTEX DEGREE
edge_density(net, loops = FALSE) # EDGE DENSITY = RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF ALL POSSIBLE EDGES. Single Number
closeness(net, mode = "all", weights = NA) # CLOSENESS = MEASURES HOW MANY STEPS IT TAKES TO REACH EVERY OTHER VERTEX FROM A GIVEN VERTEX. # RUN ON GRAPHS THAT ARE WELL CONNECTED. 
betweenness(net, directed = F, weights = NA) # BETWEENNESS = VERTEX AND EDGE BETWEENESS ARE DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A VERTEX OR AN EDGE
edge_betweenness(net, directed = FALSE, weights = NA) # EDGE BETWEENESS = NUMBER OF TIMES A NODE LIES ON THE SHORTEST PATH BETWEEN OTHER NODES. OFTEN USED TO WEIGHT EDGES
graph.strength(net) #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)

E(net)$weight <- edge_betweenness(net) # SET THE EDGE  
# UNCLEAR IF THIS IS VERY USEFUL. 
V(net)$label <- V(net)$name
V(net)$degree <- degree(net)

hist(V(net)$degree,
     col = 'green',
     xlim = c(0,20),
     main = 'Histogram of Vertex Degree',
     ylab = 'Frequency',
     xlab = 'Vertex Degree')

# ADD COMMUNITY CLUSTERING
cluster <- cluster_spinglass(net)
cluster <- cluster_edge_betweenness(net)
cluster
str(cluster)

# CREATE YOUR PLOT
# png(file=paste0("Night ", i, ".png"), width=500, height=500) #OPEN PNG FOR SAVING

lo <- layout_with_fr(net)
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)
par(mfrow=c(1,2), mar=c(0,0,0,0))

plot.igraph(net, edge.arrow.width = .25,
            edge.arrow.size = .25,
            vertex.label = NA,
            vertex.size = 5, 
            rescale=FALSE, 
            layout=lo*0.25)
plot.igraph(net, edge.arrow.width = .25,
            edge.arrow.size = .25,
            vertex.label = NA,
            vertex.size = 5, 
            rescale=FALSE, 
            layout=lo*1)


plot(net,
     vertex.size = igraph::degree(net)*2, #SET NODE SIZE AS A FUNCTION OF DEGREE CENTRALITY MULTIPLIED BY A SCALAR
     vertex.label.color = "black",
     # vertex.label.font = 1,
     # vertex.label.cex = 1.5,
     # vertex.label = NA,
     
     ### EDGE SETTINGS
     # edge.width = 2,
     edge.width = 100*(edge_attr(net)$weight)/100,
     # edge.width = E(net)$weight,
     edge.color = 'black',
     # edge.color = 'gray50',
     # edge.curved = 0.5,
     # edge.label.font = 1,
     # edge.label.cex = 1,
     
     ### COMMUNITY CLUSTERING
     mark.groups = cluster,
     mark.border = NA,
     
     ### LAYOUTS
     layout = layout_nicely(net)
     
)
# dev.off()    # SAVE THE PNG FILE






# [R, SNA, IN PROGRESS] Global network properties -------------------------

# http://pablobarbera.com/big-data-upf/html/02b-networks-descriptive-analysis.html
