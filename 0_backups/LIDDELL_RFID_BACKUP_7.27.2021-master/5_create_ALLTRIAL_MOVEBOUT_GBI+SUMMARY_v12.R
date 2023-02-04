# 5_create_ALLTRIAL_MOVEBOUT_GBI+SUMMARY.R
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# load packages
library(tidyverse)
library(readxl)
library(plyr)
library(dplyr)
library(data.table)
library(readr)
library(lubridate)

# SET WD, LOAD DATA, FORMAT TIME SERIES
# wd <- setwd("C:/Users/caleb/Box/CV_Shared_Analyses/7_LID_2020/RFID_analysis_v5")
wd <- setwd("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_analysis_v5")
data <- as.data.frame(fread("LID_2020_ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE))
metadata <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)


data$Field_Time <- as.POSIXct(data$Field_Time, format="%Y-%m-%d %H:%M:%OS")
data$Field_Time_STOP <- as.POSIXct(data$Field_Time_STOP, format="%Y-%m-%d %H:%M:%OS")

data$Field_Time <- with_tz(data$Field_Time, "UTC")
data$Field_Time_STOP <- with_tz(data$Field_Time_STOP, "UTC")

data$Zone <- paste(data$Paddock, data$Antenna, sep = "_")
trials <- unique(data$Trial)


# Create MOVEBOUT_GBI -----------------------------------------------------

trial_list <- list()
aa = trials[4]
for(aa in trials[1:length(trials)]){
  # select data
  df <- data %>%
    filter(Trial == aa) %>% 
    mutate(Mouse = paste(Strain,Sex,Name, sep = "-")) %>% 
    select(Trial, Antenna, Zone, Strain, Sex, Full_ID, Name, Mouse,  Code, noon_to_noon_day, Field_Time, Field_Time_STOP, duration_s)
  
  mice <- unique(df$Mouse)  
  
  #loop 2
  zone_list <- list()
  zones <- unique(df$Zone)
  bb = zones[1]
  for(bb in zones[1:length(zones)]){
    df1 <- df %>% 
      filter(Zone == bb)
    df1 <- df1[order(df1$Field_Time, na.last=FALSE), ]
    df1$Field_Time <- as.POSIXct(df1$Field_Time, format="%Y-%m-%d %H:%M:%OS")
    df1$Field_Time_STOP <- as.POSIXct(df1$Field_Time_STOP, format="%Y-%m-%d %H:%M:%OS")
    
    #grab start and stops
    START <- as.data.frame(df1$Field_Time)
    colnames(START) <- "x"
    STOP <- as.data.frame(df1$Field_Time_STOP)
    colnames(STOP) <- "x"
    
    # create long row of starts and stops and order them
    start_stop <- as.data.frame(rbind(START, STOP))
    
    colnames(start_stop) <- "x"
    start_stop <- as.data.frame(start_stop[order(start_stop$x, na.last=FALSE),])
    colnames(start_stop) <- "x"
    
    ## loop 3
    list <- list()
    cc=1
    for(cc in 1:(nrow(start_stop)-1)) {
      Field_Time_Start <- as.character(start_stop$x[cc], format="%Y-%m-%d %H:%M:%OS")
      Field_Time_Stop <- as.character(start_stop$x[cc+1], format="%Y-%m-%d %H:%M:%OS")
      list[[cc]] <- cbind(Field_Time_Start, Field_Time_Stop)
    }
    
    # Create GBI data frame
    gbi <- as.data.frame(do.call("rbind", list))
    
    #convert to posixct
    gbi$Field_Time_Start <- as.POSIXct(gbi$Field_Time_Start, format="%Y-%m-%d %H:%M:%OS")
    gbi$Field_Time_Stop <- as.POSIXct(gbi$Field_Time_Stop, format="%Y-%m-%d %H:%M:%OS")
    #force convert to UTC without changing time
    gbi$Field_Time_Start <- force_tz(gbi$Field_Time_Start, tzone = "UTC")
    gbi$Field_Time_Stop <- force_tz(gbi$Field_Time_Stop, tzone = "UTC")
    
    #create duration
    gbi$duration_s <- gbi$Field_Time_Stop - gbi$Field_Time_Start
    gbi_cols <- c("Trial", "Day", "Zone")
    gbi[gbi_cols] <- NA
    # create center_time which will be compared later to determine participation in grouping event
    gbi$center_time <- gbi$Field_Time_Start + (gbi$duration_s/2)
    
    ## add columns for all trial mice
    gbi[mice] <- NA
    
    #loop 4
    dd = 2
    for(dd in 1:nrow(gbi)) {
      center <- as.POSIXct(gbi$center_time[dd])
      print(paste("Processing row ",dd," out of ",nrow(gbi), " for zone ",bb," in Trial ",aa, sep=''))
      ff =1
      for(ff in 1:nrow(df1)){
        int <- interval(df1$Field_Time[ff], df1$Field_Time_STOP[ff])
        int
        if(center %within% int) {
          cool_mouse <- df1$Mouse[ff]
          # add a 1 to the mouse column present in the interaction
          gbi[[cool_mouse]][[dd]] <- 1 
          gbi$Day[dd] <- df1$noon_to_noon_day[ff]
        }
      }
    }
    
    # now, remove grouping bouts with no detected animals. 
    gbi2 <- gbi[!is.na(gbi$Day),]
    
    # add other details
    gbi2$Trial <- paste(aa)
    gbi2$Zone <- paste(bb)
    
    #replace na's with 0s
    gbi2[is.na(gbi2)] <- 0
    
    #drop center time
    zone_list[[bb]] <- subset(gbi2, select = -c(center_time))
    
  }
  
  trial_list[[aa]] <- do.call("rbind", zone_list) 
  
}


meta_short <- metadata %>% 
  mutate(Mouse = paste(strain, sex, name, sep = "-")) %>% 
  select(trial, paddock, strain, sex, name, code, Mouse,family_group)

## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
males <- meta_short %>% 
  filter(sex == "M") %>% 
  select(Mouse) %>% 
  filter(!is.na(Mouse))
male_list <- dplyr::pull(males, Mouse)

## female names list
females <- meta_short %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  select(Mouse) %>% 
  filter(!is.na(Mouse))
female_list <- dplyr::pull(females, Mouse)


## add in columns with group composition information
## this probably should be just added directly into the loop, but helpful if you have already run the loop which takes awhile. 
trial_stats <- list()
aa = 1
for(aa in 1:length(trial_list)){
  ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- trial_list[[aa]]
  current_trial <- unique(df$Trial)
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum)
  
  df3 <- df2 %>% 
    relocate(Trial, Day, Zone, Field_Time_Start, Field_Time_Stop, duration_s, m_sum, f_sum, mf_sum)
  
  write.csv(df3, paste(current_trial, "_MOVEBOUT_GBI.csv", sep=""))
}


# CREATE ALLTRIAL_MOVEBOUT_GBI_SUMMARY --------------------------------------------------
## READ IN DATA
filenames <- list.files(wd, pattern = "*MOVEBOUT_GBI.csv")

## READ IN ALL FILES
myfiles = lapply(filenames, fread)

trial_stats <- list()
aa = 1
for(aa in 1:length(trial_list)){
  ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- trial_list[[aa]]
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum)
  
  ## get mouse column names starting at col 11 (Check this and confirm)
  col_ids <- colnames(df2[,10:ncol(df2)]) ## CHANGE
  col_ids
  
  stats <- list()
  bb = col_ids[1]
  for(bb in col_ids[1:length(col_ids)]) {
    df3 <- df2 %>% 
      # select(Day, Zone, Start, End, Duration, Field_Time_START,Field_Time_STOP, (bb)) %>% 
      ## HOLY SHIT THIS TOOK AWHILE. TO LOOP THROUGH THE COLUMNS, NEED TO CHANGE STRING AS VARIABLE TO SYMBOL.
      filter((!!as.symbol(bb)) == 1) %>% 
      mutate(Name = bb) %>% 
      select(Day, Field_Time_Start, Field_Time_Stop, Zone, duration_s, Name,  m_sum, f_sum, mf_sum)
    
    ## ADD RELEVANT METADATA INFORMATION. 
    df4 <- merge(df3, meta_short, by.x = "Name", by.y = "Mouse")
    df5 <- df4 %>% 
      relocate(trial, paddock, Day, Zone, 
               strain, sex, name, code, family_group, Name, 
               Field_Time_Start, Field_Time_Stop, duration_s, m_sum, f_sum, mf_sum)
    
    stats[[bb]] <- df5
  }
  trial_stats[[aa]] <- do.call("rbind", stats)
  
}

master_stats <- do.call("rbind", trial_stats)
write.csv(master_stats, "LID_2020_ALLTRIAL_MOVEBOUT_GBI_SUMMARY.csv")

# [DEPRECATED] Create adjacency matrices -------------------------------------------------
# library(igraph)

# # select data
# df <- data %>%
#   select(Trial, Antenna, Zone, Strain, Sex, Full_ID, Name, Code, noon_to_noon_day, Field_Time, Field_Time_STOP, duration_s) %>%
#   filter(Trial == "T005")
# 
# # Remove any rows with NA's.
# df <- na.omit(df)
# df1 <- df[order(df$Field_Time, na.last=FALSE), ]
# 
# Zone <- sort(unique(df1$Zone))
# flag <- 0
# aa = 1
# for (aa in 1:length(Zone)) {
#   zone_num <- Zone[aa]
#   zonewise <- df1[which(df1$Zone==zone_num),]
#   print(paste("Processing zone ",zone_num," out of ",length(Zone),sep=''))
#   bb = 1
#   for(bb in 1:(nrow(zonewise)-1)){
#     c1 <- zonewise[bb,,drop=FALSE] #pull a focal row that will be compared to all other rows
#     for(cc in (bb+1):(nrow(zonewise))){
#       c2 <- zonewise[cc, , drop=FALSE] # pull out the rows to compare to focal row
#       if(c1$Name != c2$Name){
#         ##if mouse 2 arrives during mouse 1's bout, pull mouse 1 and mouse 2 names and zone where it occurred
#         if(c2$Field_Time < c1$Field_Time_STOP & c2$Field_Time > c1$Field_Time){
#           ## add code to pull out estimate overlap_duration_s time.
#           xtemp <- matrix(c(as.character(c1$Name), as.character(c2$Name), c1$Zone), nrow=1)
#           colnames(xtemp)<-c("ID1","ID2","zone")
#           if(flag<1){
#             socialinteractions <- xtemp
#           } else {
#             socialinteractions <- rbind(socialinteractions,xtemp)
#           }
#           flag <- flag+1
#         }
#       }
#     }
#   }
# }
# 
# # This will create a directed data frame as combinations are repeated. no time association
# socialinteractions <- as.data.frame(socialinteractions)
# socialinteractions$ID2 <- factor(socialinteractions$ID2)
# socialinteractions$ID1 <- factor(socialinteractions$ID1)
# socialinteractions$zone <- factor(socialinteractions$zone)
# summary(socialinteractions)
# # Create directed adjacency matrix
# g.unit <- (table(socialinteractions))
# # note that this will result in combinations that dont happen (have frequency of 0).
# caca <- as.data.frame(g.unit)
# caca[which(max(caca$Freq)==caca$Freq),]
# ids <- list()
# ids[[1]]<-sort(unique(caca$ID1))
# ids[[2]]<-sort(unique(caca$ID2))
# #Create the totsmagoats.
# z=Zone[1]
# flag <- 0
# for(z in Zone[1:length(Zone)]){
#   flag <- flag+1
#   diszone <- caca[which(caca$zone==z),]
#   ##change number to reflect all subjects in ids.
#   present <- matrix(diszone$Freq,nrow=20,ncol=20,dimnames = ids) ### Change this.
#   if(flag<2){
#     totsmagoats<-present
#   } else {
#     totsmagoats<-totsmagoats+present
#   }
# }
# sum(totsmagoats)
# # This is a directed adjacency matrix, where ID2 is entering the zone after ID1 is there. Rows = ID1 (first mouse). and columns equal ID2 (second mouse).
# # Weighted towards which individual comes into the tub/approaches?
# # write.csv(totsmagoats, "directed_adjacency_matrix.csv")
# # Create Undirected adjacency matrix
# summary(socialinteractions)
# df <- subset(socialinteractions, select = -c(3))
# # coerces the data into a two-column matrix format that igraph likes
# el=as.matrix(df)
# el[,1]=as.character(el[,1])
# el[,2]=as.character(el[,2])
# # turns the edgelist into a 'graph object'
# g = graph.edgelist(el,directed=FALSE)
# #create undirected adjacency matrix from edgelist
# g <- get.adjacency(g, sparse = FALSE)
# # write.csv(g, "undirected_adjacency_matrix.csv")
