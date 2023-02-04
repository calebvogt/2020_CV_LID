## analyze_TW_GBI_CSV
## Caleb C. Vogt, Cornell University

library(tidyverse)
library(asnipe)
library(igraph)
library(data.table)
library(readxl)

# SET WD, OUTPUT PATH, LOAD DATA, FORMAT TIME SERIES
wd <- setwd("C:/Users/Caleb/Box/CV_Shared_Analyses/7_LID_2020/RFID_DATA/RFID_analysis_v3")
output_fp <- paste("C:/Users/Caleb Vogt/Box/CV_Shared_Analyses/7_LID_2020/RFID_DATA/RFID_analysis_v3/0_output_plots")
metadata <- read_excel("Liddell_2020_metadata.xlsx", sheet = 1, skip = 1)


# CREATE GMM_MOUSE_STATS --------------------------------------------------

meta_short <- metadata %>% 
  select(trial, paddock, strain, sex, name, code, family_group)

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

## READ IN DATA
filenames <- list.files(wd, pattern = "*TW_GBI_DATA.csv")

## READ IN ALL FILES
myfiles = lapply(filenames, fread)

## LOOP #1: 
trial_stats <- list()
aa = 1
for(aa in 1:length(myfiles)){
  ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- myfiles[[aa]]
  df2 <- df %>% 
    mutate(m_sum = rowSums(select(., contains(male_list)))) %>% 
    mutate(f_sum = rowSums(select(., contains(female_list)))) %>% 
    mutate(mf_sum = rowSums(select(., contains(c(male_list,female_list))))) %>% 
    relocate(m_sum,f_sum,mf_sum)
  
  ## get mouse column names starting at col 10
  col_ids <- colnames(df2[,10:ncol(df2)])
  col_ids
  
  stats <- list()
  bb = col_ids[1]
  for(bb in col_ids[1:length(col_ids)]) {
    df3 <- df2 %>% 
      # select(Day, Zone, Start, End, Duration, Field_Time_START,Field_Time_STOP, (bb)) %>% 
      ## HOLY SHIT THIS TOOK AWHILE. TO LOOP THROUGH THE COLUMNS, NEED TO CHANGE STRING AS VARIABLE TO SYMBOL.
      filter((!!as.symbol(bb)) == 1) %>% 
      mutate(Name = bb) %>% 
      select(days, times, Field_Time, Zone, locations,Name,  m_sum, f_sum, mf_sum)
    
    ## ADD RELEVANT METADATA INFORMATION. 
    df4 <- merge(df3, meta_short, by.x = "Name", by.y = "name")
    df5 <- df4 %>% 
      relocate(trial, paddock, days, Zone, locations, 
               strain, sex, Name, code, family_group, 
               times, Field_Time, m_sum, f_sum, mf_sum)
    
    stats[[bb]] <- df5
  }
  trial_stats[[aa]] <- do.call("rbind", stats)
  
}

master_stats <- do.call("rbind", trial_stats)
write.csv(master_stats, "LID_2020_TW_GBI_SUMMARY.csv")


# REMAKE PLOTS PLOTS WITH GMM_MOUSE_STATS FOR FIELD DURATION ANALYSIS ETC ----------------------------------------------

## LOOP #1: CREATE INDIVIDUAL PARTICIPATION IN SOCIAL EVENT DATA FRAME FOR MOVEMENT ANALYSIS. 
trial_stats <- list()
aa = 1
for(aa in 1:length(myfiles)){
  ## PULL OUT TRIAL DATAFRAME
  df <- myfiles[[aa]]
  col_ids <- colnames(df[,10:ncol(df)])
  col_ids <- col_ids[!grepl('dud_mouse', col_ids)]
  
  stats <- list()
  bb = col_ids[1]
  for(bb in col_ids[1:length(col_ids)]) {
    
    ## CREATE SUMMED DURATION PER NIGHT 
    df2 <- df %>% 
      # select(Day, Zone, Start, End, Duration, Field_Time_START,Field_Time_STOP, (bb)) %>% 
      ## HOLY SHIT THIS TOOK AWHILE. TO LOOP THROUGH THE COLUMNS, NEED TO CHANGE STRING AS VARIABLE TO SYMBOL.
      filter((!!as.symbol(bb)) == 1) %>% 
      mutate(Name = bb) %>% 
      group_by(Name, Day) %>%
      tally(sum(Duration))
    
    stats[[bb]] <- df2
  }
  trial_stats[[aa]] <- do.call("rbind", stats)
  
}







## LOOP #1: 
aa = myfiles[[1]]
for(aa in myfiles){
  df <- aa
  col_ids <- colnames(df[,9:ncol(df)])
  col_ids <- col_ids[!grepl('dud_mouse', col_ids)]
  
  
  # SUBSET DATA AND CREATE NETWORK OBJECT -------------------------------------------------------------
  
  # SUBSET DESIRED DF INTO GBI FORMAT FOR NETWORK ANALYSES... NOTE THAT THIS GBI IS PURELY 1'S AND 0'S. DOES NOT INCLUDE 
  # THE OBSERVATIONS_PER_EVENT GBI FROM THE GMM OUTPUT. WOULD NEED TO REMAKE GMM_GBI_DATA.CSVS TO GET THIS. 
  # The third item has the same structure as the group by individual matrix, but instead of being binary
  # (0 or 1), it contains the number of observations of each individual in each event.
  
  
  for(i in 1:10) {
    
    gbi <- data %>% 
      filter(Day == i) %>%
      select(all_of(col_ids))
    
    # CREATE IDS LIST FROM COLUMN NAMES
    ids <- colnames(gbi)
    
    # TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
    # HOW DOES THIS UNDIRECTED MATRIX WEIGHT THE RELATIONSHIPS? 
    undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION
                                data_format = "GBI",
                                #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
                                association_index = "SRI")
    
    # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
    net <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE)
    
    # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
    net <- simplify(net)
    
    #VIEW VERTICES (NODES) OF GRAPH
    V(net)
    
    # VIEW VERTICES/NODE NAMES
    V(net)$name
    
    # VIEW EDGES OF GRAPH
    E(net)
    
    # VIEW EDGE WEIGHT VALUES. HOW DOES ANSIPE WEIGHT THESE EDGES? 
    E(net)$weight
    
    # PLOT THE CURRENT NETWORK
    # plot(net)
    
    # PLOT THE CURRENT NETWORK WITH STANDARD WEIGHTS. WHAT IS THE APPROPRIATE EDGE WEIGHTING METRIC?
    # plot(net,
    #      vertex.label = V(net)$name, 
    #      layout = layout.fruchterman.reingold,
    #      edge.color = "black",
    #      edge.width = E(net)$weight # NOTE: THE EDGES HERE ARE SMALL DUE TO THE ASNIPE ALGO. 
    # )
    
    
    # GET NETWORK MEASURES ----------------------------------------------------
    
    # SUMMARY AND GRAPH INFORMATION
    summary(net)
    print_all(net)
    
    # GET GRAPH CENTRALITY MEASURES AT NODE AND NETWORK LEVEL
    # RES = NODE CENTRALITY, #CENTRALIZATION = GRAPH CENTRALITY, THEORETICAL MAX = ALL POSSIBLE EDGES. 
    centr_degree(net, mode = "all")
    
    #VERTEX DEGREE
    degree(net, mode="all") 
    
    # EDGE DENSITY = RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF ALL POSSIBLE EDGES. Single Number
    edge_density(net, loops = FALSE)
    
    # CLOSENESS = MEASURES HOW MANY STEPS IT TAKES TO REACH EVERY OTHER VERTEX FROM A GIVEN VERTEX. 
    # RUN ON GRAPHS THAT ARE WELL CONNECTED. 
    closeness(net, mode = "all", weights = NA)
    
    # BETWEENNESS = VERTEX AND EDGE BETWEENESS ARE DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A VERTEX OR AN EDGE
    betweenness(net, directed = F, weights = NA)
    
    # EDGE BETWEENESS = NUMBER OF TIMES A NODE LIES ON THE SHORTEST PATH BETWEEN OTHER NODES. OFTEN USED TO WEIGHT EDGES
    edge_betweenness(net, directed = FALSE, weights = NA)
    E(net)$weight <- edge_betweenness(net)
    
    # MEASURES TO ADD: NETWORK DIAMETER, 
    
    
    
    
    # UNCLEAR IF THIS IS VERY USEFUL. 
    V(net)$label <- V(net)$name
    V(net)$degree <- degree(net)
    
    
    #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
    graph.strength(net)
    
    # GRAPHING THE NETWORK  ---------------------------------------------------------
    
    # DEGREE CENTRALITY HISTOGRAM 
    # hist(V(net)$degree,
    #      col = 'green',
    #      main = 'Histogram of Vertex Degree',
    #      ylab = 'Frequency',
    #      xlab = 'Vertex Degree')
    
    # ADD COMMUNITY CLUSTERING
    # cluster <- cluster_spinglass(net)
    # cluster <- cluster_edge_betweenness(net)
    # cluster
    # str(cluster)
    
    # CREATE YOUR PLOT
    png(file=paste0("Night ", i, ".png"), width=500, height=500) #OPEN PNG FOR SAVING
    plot(net,
         
         ### VERTEX SETTINGS
         # vertex.color = vertex_attr<- (net)$cor, ### GAH see this webpage on how to assign vertex attributes.
         # https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwig5ru0rODuAhUAF1kFHZlxDWsQFjAFegQIAxAC&url=https%3A%2F%2Fmariliagaiarsa.weebly.com%2Fuploads%2F3%2F8%2F6%2F2%2F38628397%2Figraphtutorialeng.html&usg=AOvVaw0tsqEmkx30sLVNf7jsrTK5
         vertex.color = "lightblue",
         # vertex.color = "red",
         # vertex.color = rainbow(10), # ADJUST FOR NUMBER OF NODES
         # vertex.size = 20,
         vertex.size = 50,
         # vertex.size = igraph::degree(net)*5, #SET NODE SIZE AS A FUNCTION OF DEGREE CENTRALITY MULTIPLIED BY A SCALAR
         vertex.label.color = "black",
         vertex.label.font = 1,
         vertex.label.cex = 1.5,
         # vertex.label = NA,
         
         ### EDGE SETTINGS
         # edge.width = 2,
         # edge.width = 100*(edge_attr(net)$weight)/100,
         edge.width = E(net)$weight,
         edge.color = 'black',
         # edge.color = 'gray50',
         edge.curved = 0.5,
         # edge.label.font = 1,
         # edge.label.cex = 1,
         
         ### COMMUNITY CLUSTERING
         # mark.groups = cluster, 
         # mark.border = NA,
         
         ### LAYOUTS
         layout = layout_in_circle(net, order = ids), # SORT ALPHABETICALLY FOR REPEATED GRAPHS ACROSS DAYS
         # layout = layout_nicely(net),
         # layout = layout.graphopt(net),
         # layout = layout.fruchterman.reingold(net),
         # layout = layout.kamada.kawai(net),
         # layout = layout_as_bipartite(net),
         # layout = layout_as_star(net),
         # layout = layout_as_tree(net),
         
         ### PLOT PROPERTIES
         main = paste0("Night ", i)
    )
    dev.off()    # SAVE THE PNG FILE
    
    
  }
  
  ## note: need to add vertex color option to make males light blue, females red. 
  # 
  # # COMMUNITY DETECTION GRAPH
  # cnet <- cluster_edge_betweenness(net)
  
  # PLOT COMMUNITY DETECTION GRAPH. 
  # plot(cnet,
  #      net)
  # 
  
}


# MALE AND FEMALE FULL NETWORK FULL TRIAL --------------------------------------

# SET WD AND LOAD DATA
wd <- setwd("G:/My Drive/PROJECTS & WRITING/Liddell 2020  - Mus musculus social structure/Liddell_2020_RFID/Liddell_2020_RFID_full_10day_data")

metadata <- read_excel("Liddell_2020_metadata.xlsx", sheet = 1, skip = 1)
metadata <- metadata[1:140,]
metadata <- metadata %>%
  select(trial, name, code, sex, field_age, pre_mass)

# CHANGE DATA/TRIAL SET YOU WANT TO WORK WITH. 
data <- as.data.frame(fread("T001_GMM_GBI_DATA.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
data <- subset(data, select = -c(V1))


attributes <- metadata %>% 
  filter(trial == "T001")  ## CHANGE THIS


col_ids <- unique(attributes$name)

# SUBSET DATA AND CREATE NETWORK OBJECT -------------------------------------------------------------


gbi <- data %>% 
  select(all_of(col_ids))

# CREATE IDS LIST FROM COLUMN NAMES
ids <- colnames(gbi)

# TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
# HOW DOES THIS UNDIRECTED MATRIX WEIGHT THE RELATIONSHIPS? 
undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION
                            data_format = "GBI",
                            #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
                            association_index = "SRI")

# CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
net <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE)

# SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
net <- simplify(net)

#VIEW VERTICES (NODES) OF GRAPH
V(net)

# VIEW VERTICES/NODE NAMES
V(net)$name

# VIEW EDGES OF GRAPH
E(net)

# VIEW EDGE WEIGHT VALUES. HOW DOES ANSIPE WEIGHT THESE EDGES? 
E(net)$weight

# PLOT THE CURRENT NETWORK
plot(net)

# PLOT THE CURRENT NETWORK WITH STANDARD WEIGHTS. WHAT IS THE APPROPRIATE EDGE WEIGHTING METRIC?
plot(net,
     vertex.label = V(net)$name, 
     layout = layout.fruchterman.reingold,
     edge.color = "black",
     edge.width = E(net)$weight # NOTE: THE EDGES HERE ARE SMALL DUE TO THE ASNIPE ALGO. 
)


# GET NETWORK MEASURES ----------------------------------------------------

# SUMMARY AND GRAPH INFORMATION
summary(net)
print_all(net)

# GET GRAPH CENTRALITY MEASURES AT NODE AND NETWORK LEVEL
# RES = NODE CENTRALITY, #CENTRALIZATION = GRAPH CENTRALITY, THEORETICAL MAX = ALL POSSIBLE EDGES. 
centr_degree(net, mode = "all")

#VERTEX DEGREE
degree(net, mode="all") 

# EDGE DENSITY = RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF ALL POSSIBLE EDGES. Single Number
edge_density(net, loops = FALSE)

# CLOSENESS = MEASURES HOW MANY STEPS IT TAKES TO REACH EVERY OTHER VERTEX FROM A GIVEN VERTEX. 
# RUN ON GRAPHS THAT ARE WELL CONNECTED. 
closeness(net, mode = "all", weights = NA)

# BETWEENNESS = VERTEX AND EDGE BETWEENESS ARE DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A VERTEX OR AN EDGE
betweenness(net, directed = F, weights = NA)

# EDGE BETWEENESS = NUMBER OF TIMES A NODE LIES ON THE SHORTEST PATH BETWEEN OTHER NODES. OFTEN USED TO WEIGHT EDGES
edge_betweenness(net, directed = FALSE, weights = NA)

# SET THE EDGE  
E(net)$weight <- edge_betweenness(net)

# MEASURES TO ADD: NETWORK DIAMETER, 


# UNCLEAR IF THIS IS VERY USEFUL. 
V(net)$label <- V(net)$name
V(net)$degree <- degree(net)


#GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
graph.strength(net)

# GRAPHING THE NETWORK  ---------------------------------------------------------

# DEGREE CENTRALITY HISTOGRAM 
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
plot(net,
     
     ### VERTEX SETTINGS
     # vertex.color = vertex_attr <- (net)$cor, ### GAH see this webpage on how to assign vertex attributes.
     # https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwig5ru0rODuAhUAF1kFHZlxDWsQFjAFegQIAxAC&url=https%3A%2F%2Fmariliagaiarsa.weebly.com%2Fuploads%2F3%2F8%2F6%2F2%2F38628397%2Figraphtutorialeng.html&usg=AOvVaw0tsqEmkx30sLVNf7jsrTK5
     # vertex.color = "lightblue",
     # vertex.color = "red",
     vertex.color = rainbow(20), # ADJUST FOR NUMBER OF NODES
     # vertex.size = 20,
     # vertex.size = 50,
     vertex.size = igraph::degree(net), #SET NODE SIZE AS A FUNCTION OF DEGREE CENTRALITY MULTIPLIED BY A SCALAR
     vertex.label.color = "black",
     # vertex.label.font = 1,
     # vertex.label.cex = 1.5,
     # vertex.label = NA,
     
     ### EDGE SETTINGS
     edge.width = 2,
     # edge.width = 100*(edge_attr(net)$weight)/100,
     # edge.width = E(net)$weight,
     # edge.color = 'black',
     # edge.color = 'gray50',
     # edge.curved = 0.5,
     # edge.label.font = 1,
     # edge.label.cex = 1,
     
     ### COMMUNITY CLUSTERING
     mark.groups = cluster,
     mark.border = NA,
     
     ### LAYOUTS
     # layout = layout_in_circle(net, order = ids), # SORT ALPHABETICALLY FOR REPEATED GRAPHS ACROSS DAYS
     layout = layout_nicely(net),
     # layout = layout.graphopt(net),
     # layout = layout.fruchterman.reingold(net),
     # layout = layout.kamada.kawai(net),
     # layout = layout_as_bipartite(net),
     # layout = layout_as_star(net),
     # layout = layout_as_tree(net),
     
     ### PLOT PROPERTIES
     # main = paste0("Night ", i)
)
# dev.off()    # SAVE THE PNG FILE

