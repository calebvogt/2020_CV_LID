## analyze_ALLTRIAL_RFID_DATA

library(tidyverse)
library(data.table)
library(plyr)
library(dplyr)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(Hmisc)
library(plotly)
library(ggiraph)
library(readr)
library(reshape)
library(lubridate)
library(scales)
library(plot.matrix)
library(transformr)
library(rstatix)


# SET WD, LOAD DATA, FORMAT TIME SERIES
wd <- setwd("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_DATA/RFID_analysis_v4")
output_fp <- paste("C:/Users/Caleb Vogt/Desktop")

data <- as.data.frame(fread("LID_2020_ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE))
data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS")
df <- data

## Clean data if desired
## CLEAN DATA DOWN TO 10 DAYS OF OBSERVATION BETWEEN 6PM AND 6AM
df <- data %>%
  filter(noon_to_noon_day == c(1:10))
  # filter(noon_to_noon_day == c(1:10), Time <= "06:00:00" | Time >= "18:00:00")

# [Figure S2] Within mouse, within tub RFID hit interval ---------------------------

## use this to determine your WINDOW BOUT time threshold value. Critical!
df$Zone <- df$Antenna
ids <- unique(df$Name)
data_list <- list()
aa = ids[4]
aa = "Cole"
for(aa in ids[1:length(ids)]){
  
  print(paste("Processing mouse ",aa, " out of ", length(ids), sep=''))
  
  ## create df of percent time
  df1 <- df %>% 
    filter(Name == aa) 
  
  zones <- unique(df1$Zone)
  bb = zones[1]
  zone_list <- list()
  for(bb in zones[1:length(zones)]){
    zone_list[[bb]] <- df1 %>% 
      filter(Zone == bb) %>% 
      select(Trial, Strain, Sex, Name, Zone, Field_Time) %>% 
      mutate(diff = Field_Time - lag(Field_Time), diff_secs = as.numeric(diff, units = 'secs'))
  }
  
  
  data_list [[aa]] <- do.call("rbind", zone_list)
}

df2 <- do.call("rbind", data_list)

sdat <- summary(df2$diff_secs)
sdat
summStr <- paste(names(sdat), format(sdat, digits = 4), collapse = "; ")
summStr
op <- par(mar = c(7,4,4,2) + 0.1)
# hist(df5$diff_secs)
hist(df2$diff_secs, 
     xlim = c(0, 400),
     breaks = 300000,
     # main = stuff,
     main = "",
     xlab = "Within Tub inter-read interval (s)"
)
title(sub = summStr, line = 5.5)
par(op)

## get time interval where 95% of intervals below that value. 
sort(df2$diff_secs)[0.95*length(df2$diff_secs)]

# 10 days, all time = 258 seconds
# 10 days, 6pm-6AM = 345 seconds. We use this as the more conservative estimate of nocturnal transition times. 


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




# PNG: ERROR: 24H B6 MALE CIRCADIAN ------------------------------------
## ISSUE 1: REMOVE INDIVIDUALS KNOWN TO BE NESTING UNDER THE RFID READERS. ESPECIALLY TRUE FOR OBS? 
# df$Time <- as.POSIXct(df$Time, format="%H:%M:%OS")
df <- data %>% 
  filter(noon_to_noon_day == c(1:10),
         Strain == "C57", 
         Sex == "M")

df <- data %>% 
  select(Strain, Sex, noon_to_noon_day, Time) %>% 
  filter(Strain == "C57", Sex == "M") %>% 
  filter(noon_to_noon_day == 2)

df$Time <- as.numeric(as.character(df$Time))

ggplot(df) +
  aes(x = Time) +
  geom_bar()
# stat_count(width = 0.5)
p


filter(noon_to_noon_day == 1) %>% 
  filter(Strain)
Strain == "C57", 
Sex == "M")




## I fucked this up somehow. 
png(file = paste0(output_fp, "/", "B6_M_circadian.png"))
hist(df$Time, 
     breaks = "hours", 
     col = "red", 
     xlab = "Hours",
     ylab = "RFID Reads",
     # ylim = c(0,25000),
     main = "OB Males",
     freq = TRUE)
## SAVES TO CURRENT WORKING DIRECTORY. MOVE FILES TO 0_OUTPUT_PLOTS 
dev.off()


# PNG: TW CIRCULAR SOCIAL NETWORKS, WEIGHTED --------------------------------------------
library(sna)
library(timeordered)
library(ape)
library(asnipe)
library(igraph)

sex_list <- c("M","F")
trial_list <- unique(clean$Trial)
bb=trial_list[1]
for(bb in trial_list[1:length(trial_list)]) {
  df <- clean %>% 
    filter(Trial == bb)
  strain <- unique(df$Strain)
  aa=sex_list[1]
  for(aa in sex_list[1:length(sex_list)]) {
    # PULL RFID READS BY TRIAL ETC. 
    df1 <- df %>% 
      filter(Sex == aa) %>% 
      select(noon_to_noon_day, 
             Time_sec,
             Name, 
             Antenna)
    
    ## RENAME COLUMNS
    colnames(df1) = c("Date",
                      "Time", 
                      "ID",
                      "Location")
    
    df1$Location <- as.character(df1$Location)
    df1 <- unique(df1)
    ids <- sort(unique(df1$ID))
    ids
    i=8
    for(i in 1:10){
      group_by_individual <- get_associations_points_tw(df1, 
                                                        time_window = 60,
                                                        which_days = i,
                                                        which_locations = NULL) #
      
      gbi <- group_by_individual[[1]]
      # times <- group_by_individual[[2]]
      # days <- group_by_individual[[3]]
      # locations <- group_by_individual[[4]]
      # 
      # CREATE UNDIRECTED WEIGHTED ASSOCIATION MATRIX
      association_data <- gbi
      undir_matrix <- get_network(association_data,
                                  data_format = "GBI",
                                  #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
                                  association_index = "SRI")  
      
      # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX
      net_graph <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE)
      
      # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES.
      net_graph <- simplify(net_graph)
      
      #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
      strengths <- graph.strength(net_graph)
      strengths
      V(net_graph)$label <- V(net_graph)$name
      E(net_graph)$weight <- edge.betweenness(net_graph)
      
      ## CHANGE SEXY_COLOR
      if(aa == "M"){
        sexy_color <- c("lightblue")
      } else {
        sexy_color <- c("red")
      }
      
      ## GET IDS OF INDIVIDUALS WHO WERE READ ON THE CORRECT NIGHT
      df2 <- df1 %>% 
        filter(Date == i) %>% 
        select(ID)
      
      ids <- sort(unique(df2$ID))
      
      tryCatch({
        ## ERROR: T002_C57_F_NIGHT_2 ONLY SPITS 9X9 MATRIX. 
        ## REASON: NETGRAPH ONLY SHOWS NODES WHERE THERE WERE ACTUAL RFID READS. 
        ## UNCONNECTED NODES HAVE RFID READS, BUT DO NOT HAVE INTERACTIONS ACCORDING TO THE TW. 
        png(file=paste0(output_fp, "/", bb, "_", strain, "_", aa, "_Night_", i, "_TW_Social_Network.png" ), width=500, height=500) #OPEN PNG FOR SAVING
        plot(net_graph,
             layout = layout_in_circle(net_graph, order = ids), #sort alphabetically?
             # layout = layout_in_circle(net_graph, order = order(V(net_graph))),
             vertex.color = sexy_color,
             vertex.size = 50,
             vertex.label.color = "black",
             vertex.label.font=1,
             vertex.label.cex = 1.5,
             edge.color = 'black',
             edge.width=E(net_graph)$weight,
             # edge.label.font =1,
             main = paste0(strain, " Night ", i),     )  # PLOT THE NETWORK
        dev.off()    # SAVE THE PNG FILE
        print(paste0(bb, "_", strain, "_", aa, "_Night_", i, "_TW_Social_Network.png" ))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) ## COULD ALSO APPLY THE ALTERNATIVE NET GRAPH FOR THIS ISSUE. 
      
    }
    
  }
  
}

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
  labs(title = "reads per mouse",
       subtitle = "",
       x = "Trial",
       y = "Reads per mouse",
       caption = "") +
  theme_classic() 
p

p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5)

ggsave("12h Trial Reads.png", 
       # plot=p, 
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)

# PNG: 12 HOUR STRAIN/SEX READS, VIOLIN --------------------------------------
library(ggplot2)

df <- clean %>% 
  group_by(Strain, Sex, Name) %>% 
  tally()

p <- ggplot(df) +
  aes(x = Strain, y = n, fill = Sex) +
  geom_violin(trim = TRUE) + 
  labs(title = "Total RFID Reads",
       subtitle = "",
       x = "Strain",
       y = "Total Reads",
       caption = "") +
  theme_classic() 
p

## ADD STATS
p + stat_compare_means(method = "t.test")

ggsave("12h Sex_strain reads.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)

# PNG: 12 HOUR C57 READS, VIOLIN --------------------------------------
df <- clean %>% 
  filter(Strain == "C57") %>% 
  group_by(Sex, Name) %>% 
  tally()

p <- ggplot(df) +
  aes(x = Sex, y = n, fill = Sex) +
  geom_violin(trim = TRUE) + 
  geom_boxplot(width=0.05) +
  labs(title = "C57 Male and Female Reads",
       subtitle = "12 hour data",
       x = "Sex",
       y = "Total Reads",
       caption = "**P<0.005, T-test") +
  theme_classic() 
p

## ADD STATS
p + stat_compare_means(method = "t.test")
p + geom_signif(comparisons = list(c("F", "M")), 
                map_signif_level=TRUE, 
                test = "t.test")

ggsave("12h C57 Male & Female Reads.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)

# PNG: 12 HOUR NYOB READS, VIOLIN --------------------------------------
df <- clean %>% 
  filter(Strain == "NYOB") %>% 
  group_by(Sex, Name) %>% 
  tally()

p <- ggplot(df) +
  aes(x = Sex, y = n, fill = Sex) +
  geom_violin(trim = TRUE) + 
  geom_boxplot(width=0.05) +
  labs(title = "NYOB Male and Female Reads",
       subtitle = "12 hour data",
       x = "Sex",
       y = "Total Reads",
       caption = "T-test") +
  theme_classic() 
p

## ADD STATS
p + stat_compare_means(method = "t.test")
p + geom_signif(comparisons = list(c("F", "M")), 
                map_signif_level=TRUE, 
                test = "t.test")


ggsave("12h NYOB Male & Female Reads.png", 
       # plot = p,
       width = 5, 
       height = 4, 
       dpi = 300, 
       units = "in", 
       device='png', 
       path = output_fp)




# PNG: 12 HOUR C57 NIGHTLY READS, BAR -------------------------------------
library(ggsignif)

## CREATE DATA FRAME
df <- clean %>% 
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


