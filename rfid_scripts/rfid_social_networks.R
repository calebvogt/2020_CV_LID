## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(asnipe)
library(sna)
library(timeordered)
library(igraph)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")
rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
rfid_data <- rfid_data %>% 
  filter(noon_day %in% 1:12) %>% # traps dropped on 5/19 or day 13. so data up until day 12 is good. 
  mutate(sex_phase = paste(sex, phase, sep = "_"))


# Input data must be of the following form: Date is an integer for day (usually starting at 1 on the first
# day). Time are the number of seconds elapsed from the start (continuous across all dates). ID is a
# unique character string for each individual. Location is a unique character string for each location.
# ggsave(p, filename = "output/rfid_data_distance_travelled.png", device = "png", bg = "white")

df <- rfid_data %>%  
  filter(trial == "T001") %>% 
  select(noon_day, time_sec, name, zone) %>% 
  dplyr::rename(Date = noon_day, Time = time_sec, ID = name, Location = zone)
  
df$Location <- as.character(df$Location)
df <- unique(df) # remove RFID reads that happen in less than 1 sec of each other (thus repeated, since we drop ms)
ids <- sort(unique(df$ID))
ids
    
i=2
for(i in 1:max(df$Date)){ #daily networks. can also do hourly, but would need to adjust code above. See circadian code. 
  print(paste("Processing time", i))
  group_by_individual <- get_associations_points_tw(df, time_window = 50, which_days = i, which_locations = NULL) #use the 95% selection threshold
  gbi <- group_by_individual[[1]]
  times <- group_by_individual[[2]]
  days <- group_by_individual[[3]]
  locations <- group_by_individual[[4]]
  
  # CREATE UNDIRECTED WEIGHTED ASSOCIATION MATRIX
  undir_matrix <- get_network(gbi, data_format = "GBI", association_index = "SRI")  #simple ratio index, halfweight index (HWI)
  g <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX
  g <- simplify(g) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES.
  
  #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
  strengths <- graph.strength(g)
  V(g)$label <- V(g)$name
  # V(g)$color <- ifelse(V(g)$name %in% male_names, "lightblue", "red")
  E(g)$weight <- edge.betweenness(g)
  
  # ADD COMMUNITY CLUSTERING
  cc <- cluster_optimal(g)
  
  ids <- sort(unique(df$ID))
  # tryCatch({ #error, netgraph only shows nodes with actual reads, unconnected nodes have reads, but no TW interactions
    png(filename = paste0("output/rfid_data", "_social_net","_T001", "_Night_", i,".png"), width=1000, height=1000)
    plot(g,
         # layout = layout_in_circle(g),
         layout = layout.fruchterman.reingold(g),
         # layout = layout_in_circle(g, order = order(V(g))),
         ### COMMUNITY CLUSTERING
         mark.groups = cc,
         mark.border = NA,
         # vertex.color = sexy_color,
         # vertex.size = 50,
         # vertex.label.color = "black",
         # vertex.label.font=1,
         # vertex.label.cex = 1.5,
         # edge.color = 'black',
         edge.width=E(g)$weight,
         # edge.label.font =1,
         main = paste0("Night ", i)
         )  
    dev.off()    
    print(paste("Complete"))
  # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}



# PNG: REINGOLD DAILY SOCIAL NETWORK GRAPHS, ALL MICE ------------------------------
library(sna)
library(timeordered)
library(asnipe)
library(igraph)
library(dplyr)

trial_list <- unique(data$Trial)
bb=trial_list[1]
for(bb in trial_list[1:length(trial_list)]) {
  df <- data %>%  ## change this to clean if needed. 
    filter(Trial == bb)
  
  ## Set male names
  male_names <- df %>% 
    filter(Sex == "M") %>% 
    .$Name %>% 
    unique()
  
  ## Set male names
  female_names <- df %>% 
    filter(Sex == "F") %>% 
    .$Name %>% 
    unique()
  
  strain <- unique(df$Strain)
  # PULL RFID READS BY TRIAL ETC. 
  df1 <- df %>% 
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
  days <- unique(df1$Date)
  ids <- sort(unique(df1$ID))
  ids
  
  i=1
  for(i in 1:length(days)){
    group_by_individual <- get_associations_points_tw(df1, 
                                                      time_window = 60,
                                                      which_days = i,
                                                      which_locations = NULL) #
    
    ## old code ^^^
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
    
    
    
    # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
    g <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE)
    
    # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
    g <- simplify(g)
    
    #VIEW VERTICES (NODES) OF GRAPH
    V(g)
    
    # VIEW VERTICES/NODE NAMES
    V(g)$name
    
    V(g)$color
    
    # VIEW EDGES OF GRAPH
    E(g)
    
    # VIEW EDGE WEIGHT VALUES. HOW DOES ANSIPE WEIGHT THESE EDGES? 
    E(g)$weight
    
    
    
    # UNCLEAR IF THIS IS VERY USEFUL. 
    V(g)$label <- V(g)$name
    V(g)$degree <- degree(g)
    
    
    #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
    graph.strength(g)
    
   
    
    
    # CREATE YOUR PLOT
    png(file=paste0(output_fp, "/", "All_Mouse_Social_Network_", bb, "_Night_", i, ".png"), width=1000, height=1000) #OPEN PNG FOR SAVING
    plot(g,
         ### VERTEX SETTINGS
         # vertex.color = rainbow(16), # ADJUST FOR NUMBER OF NODES
         vertex.color = V(g)$color,
         vertex.size = igraph::degree(g)*3, #SET NODE SIZE AS A FUNCTION OF DEGREE CENTRALITY MULTIPLIED BY A SCALAR
         vertex.label.color = "black",
         
         ### EDGE SETTINGS
         ## divided by scalar if too large
         edge.width = E(g)$weight*3,
         
         
         
         ### LAYOUTS
         layout = layout_nicely(g),
         
         ### PLOT PROPERTIES
         main = paste0(bb, " Night ", i)
    )
    
    # SAVE THE PNG FILE
    dev.off() 
  }
}



