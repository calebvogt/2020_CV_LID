## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)
library(ggraph)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

el <- read.csv("data/ALLTRIAL_MOVEBOUT_edgelist.csv")
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)
m_pas <- read.csv("data/priority_access_scores_males.csv")
f_pas <- read.csv("data/priority_access_scores_females.csv")
pas <- rbind(m_pas,f_pas)

df <- el %>% 
  filter(trial == 'T007')

trial_info <- print(unique(df$trial))

net_stats_list <- list()
vertex_stats_list <- list()
node_stats_list <- list()

i=1
# for(i in 1:10) { # loop through days
  
  df2 <- df %>%  # select mouse columns
    # filter(day == i) %>% 
    select(ID1, ID2, sum_duration_s)
  
  ## scale values
  # df2$sum_duration_s_scaled <- scales::rescale(df2$sum_duration_s, to = c(1,10))
  
  ## create directed adjaceny matrix with duration weights
  # directed_graph <- graph.data.frame(df2, directed = TRUE)
  
  # create directed adjaceny matrix with duration weights
  undirected_graph <- graph.data.frame(df2, directed = F)
  
  ## add scaled sum_duration_s_scaled 
  # dir.adj.mat <- as_adjacency_matrix(directed_graph, type = "both", names = TRUE, attr = 'sum_duration_s_scaled', sparse = FALSE)
  # View(dir.adj.mat)
  
  ## add scaled sum_duration_s
  dir.adj.mat <- as_adjacency_matrix(undirected_graph, type = "both", names = TRUE, attr = 'sum_duration_s', sparse = FALSE)
  write.csv(dir.adj.mat, "edge_sn_duration_s_T007_full.csv")
  
  
  ## create undirected adjaceny matrix
  undirected_graph<- as.undirected(directed_graph, mode = "collapse", edge.attr.comb = "sum")
  undir.adj.mat <- as_adjacency_matrix(undirected_graph, type = "both", names = TRUE, attr = 'sum_duration_s_scaled', sparse = FALSE)
  View(undir.adj.mat)
  
  g <- graph.adjacency(undir.adj.mat, mode="undirected", weighted = T, diag = F) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
  g <- simplify(g) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
  
  # Create the network 
  social_net <- ggraph(g, layout = "stress") +                                                                                                         
    geom_node_point(size = 2) +                                         
    geom_node_text(aes(label = name), nudge_y = 0.05, nudge_x = 0.2)+ 
    geom_edge_link() +
    theme_void()
  
  # Render the network 
  show(social_net)
  
  ##
  V(g)$day <- i 
  V(g)$sex <- as.character(metadata$sex[match(V(g)$name,metadata$name)])
  V(g)$label <- V(g)$name # Add names to graph graphs (full 10 day plot)
  # V(g)$family_group = as.character(meta_short$family_group[match(V(g)$name,meta_short$name)])
  V(g)$color = V(g)$sex #assign the "Sex" attribute as the vertex color
  
  #C57 trials, Change
  V(g)$color = gsub("F","sienna1",V(g)$color) #Females will be red
  V(g)$color = gsub("M","sienna",V(g)$color) #Males will be blue
  
  #Outbred trials, Change
  # V(g)$color = gsub("F","skyblue",V(g)$color) #Females will be red
  # V(g)$color = gsub("M","skyblue4",V(g)$color) #Males will be blue
  
  # get days priority access score
  pas1 <- pas %>% 
    filter(noon_day == i) %>% 
    dplyr::select(name, csum_daily_capture_penalty)
  V(g)$node_priority_access_score <- as.character(pas1$csum_daily_capture_penalty[match(V(g)$name,pas1$name)])
  
  ## Vertex Measures
  V(g)$node_degree_centrality <- degree(g, mode="all") #node degree centrality, raw score. 
  V(g)$node_eigen_centrality <- eigen_centrality(g)$vector # measure of being connected to the well connected. not clear if this should be used or centr_eigen$vector
  V(g)$node_closeness_centrality <- closeness(g, mode = "all", normalized = TRUE, weights = NA) # measure of # steps required to access every other node. 
  V(g)$node_betweeness_centrality <- betweenness(g, directed = FALSE, weights = NA) # node BETWEENNESS = DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A VERTEX
  V(g)$node_edge_strength <- graph.strength(g) # node edge strength = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
  V(g)$node_page_rank <- page_rank(g)$vector
  V(g)$node_authority_score <- authority_score(g)$vector
  
  # create vertex stats
  vertex_stats <- igraph::as_data_frame(g, "vertices")
  vertex_stats_list[[i]] <- vertex_stats
  
  ## Edge measures
  # sort(edge_betweenness(g, directed = FALSE, weights = NA)) # EDGE BETWEENESS = DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A particular edge
  
  # Set up Network stats
  net_stats <- data.frame(matrix(ncol = 1,nrow = 1)) # create empty dataframe
  
  ## Network  Measures
  graph_centrality = centr_degree(g, mode = "all", loops = T, normalized = T) # https://igraph.org/r/doc/centr_degree.html
  net_stats$net_centrality <- graph_centrality$centralization #net level centrality score based on node level centrality. https://igraph.org/r/doc/centralize.html
  net_stats$net_eigen_centrality <- centr_eigen(g, directed =F, scale = T, normalized = T)$value #https://igraph.org/r/doc/centr_eigen.html
  net_stats$net_mean_dist <- mean_distance(g, directed = FALSE) # mean distance = avg. # of edges between any two nodes in network. 
  net_stats$net_edge_density <- edge_density(g, loops = FALSE) # EDGE DENSITY = RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF ALL POSSIBLE EDGES. Single Number
  net_stats$net_transitivity <- transitivity(g) # aka clustering coefficient, probability that adjacent nodes of a network are connected. 
  net_stats$net_components_num <- components(g)$no # number of network sub components including isolates
  net_stats$net_modularity_infomap <- modularity(cluster_infomap(g)) #infomap method attempts to map the flow of information in a network, and the different clusters in which information may get remain for longer periods. Similar to walktrap, but not necessarily maximizing modularity, but rather the so-called "map equation".
  net_stats$net_modularity_infomap_group_n <- length(cluster_infomap(g))
  net_stats$net_modularity_fast_greedy <- modularity(cluster_fast_greedy(g))
  net_stats$net_modularity_fast_greedy_group_n <- length(cluster_fast_greedy(g))
  net_stats$net_assortativity_priority_access_score <- assortativity(g, V(g)$node_priority_access_score)
  # net_stats$net_assortativity_family_group <- assortativity(g, as.factor(V(g)$family_group))
  
  # clustering dendrogram 
  # https://www.rdocumentation.org/packages/igraph/versions/1.2.6/topics/plot_dendrogram
  # height corresponds with strength of relationship 
  # plot_dendrogram(cluster_fast_greedy(g), use.modularity = F) 
  
  #clean nets_stats
  net_stats[1] <- NULL
  net_stats_list[[i]] <- net_stats
  
  
  svg(file=paste0(output_fp,"/", trial_info,"_SN_Day_", i, ".svg"), bg = "transparent") ## M+F plot
  par(mar = c(0.4, 0.1, 2, 0.1))
  plot(g,
       layout = layout.fruchterman.reingold,
       # layout = layout_nicely,
       vertex.size = 10,
       # vertex.size = V(g)$node_degree_centrality*10,
       # vertex.size = scales::rescale(V(g)$node_edge_strength, to = c(5,25)), #rescale vertex strength to reasonable min/max size
       # vertex.label= V(g)$name, # include labels
       vertex.label= NA, # Remove vertex labels
       # vertex.label.font = 2,
       # vertex.label.color = "black",
       # vertex.label.cex = 1,
       # vertex.label.degree = 2,
       # edge.width = scales::rescale(E(g)$weight, to = c(0.5,30)), # rescale edge weight to reasonable min/max size
       edge.color = "darkgray",
       edge.curved = 0.2,
       
       asp = 0.2
       ### COMMUNITY CLUSTERING.
       # Decide if you actually want this. Maybe for full 10 day plot only.
       # mark.groups = cluster_fast_greedy(g),
       # mark.border = "darkgray"
  )
  title(paste0("Day ", i), cex.main=3) #
  dev.off()
  # tkplot(g) # gui adjustments. cool!
  
  
  ############## PLOT 2: Circular Network ##################
  # x <- get.adjacency(g)
  # plot(g)
  # graph.strength(g) #GET NODE STRENGTHS FOR THIS NETWORK. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
  # V(g)$label <- V(g)$name
  # V(g)$degree <- degree(g)
  # 
  # svg(file=paste0(output_fp,"/","P2_Day_", i, ".svg", bg = "transparent"))
  # par(mar = c(0.4, 0.1, 2, 0.1))
  # plot.igraph(g,
  #             vertex.color = "lightblue", #change
  #             # vertex.color = "red",
  #             vertex.size = 50, #20
  #             # vertex.size = igraph::degree(g)*5, #SET NODE SIZE AS A FUNCTION OF DEGREE CENTRALITY MULTIPLIED BY A SCALAR
  #             vertex.label.color = "black",
  #             vertex.label.font = 4, #changes font type
  #             vertex.label.cex = 1.5, #0.75
  #             edge.width = E(g)$weight*100, #maintain original weights
  #             edge.color = 'black',
  #             edge.curved = 0.5,
  #             layout = layout_in_circle(g, order = ids) # SORT ALPHABETICALLY FOR REPEATED GRAPHS ACROSS DAYS
  # )
  # title(paste0("Day ", i), cex.main=3)
  # dev.off() 
  
}
# Get net Stats and copy the output directly 
net_stats <- do.call("rbind", net_stats_list)
# options(scipen = 999)
write.table(net_stats, "clipboard-16384", sep="\t", row.names=F, col.names = F) 

## Get vertex stats # COPY THE OUTPUT TO THE CLIPBOARD, paste directly into excel. 
vertex_stats <- do.call("rbind", vertex_stats_list)
vertex_stats <- merge(vertex_stats, meta_short, by.x = "name", by.y = "name")
vertex_stats <- vertex_stats %>% 
  dplyr::rename(sex = sex.x, family_group = family_group.x) %>% 
  dplyr::select(trial, strain, sex, name, code, family_group, day, 
                node_priority_access_score, node_degree_centrality, node_eigen_centrality, 
                node_closeness_centrality, node_betweeness_centrality, node_edge_strength, 
                node_page_rank, node_authority_score) %>% 
  filter(!(name == "George")) %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & day >= 10)) %>% 
  arrange(name, day)
write.table(vertex_stats, "clipboard-16384", sep="\t", row.names=F, col.names = F) 
