## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## change this whole thing into a loop per trial, will be much much easier and doesnt take too long to run. 
library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

dd <- paste(getwd(), "data", sep = "/")
filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
social_data = lapply(filenames, fread) ## READ IN ALL FILES
meta <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

# clean social data for triaged mice from social interaction bouts. 
# note that this cleaning step merely deletes columns and adds 0s where appropriate. Does not adjust GBI summary information (m_sum, f_sum, mf_sum)
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  ## code here. 
  df$strain <- ifelse(df$trial == "T001", "C57", 
                      ifelse(df$trial == "T002", "C57", 
                             ifelse(df$trial == "T003", "C57",
                                    ifelse(df$trial == "T004", "Outbred",
                                           ifelse(df$trial == "T005", "Outbred",
                                                  ifelse(df$trial == "T006", "C57",
                                                         ifelse(df$trial == "T007", "Outbred", NA)))))))
  
  df2 <- df %>% 
    #T004: George only mouse to cross between trials on Day 3. triage completely (drop column)
    dplyr::select(-one_of(c( "V1", "NYOB-M-George"))) %>% 
    #T003: Anubis visually confirmed dead by seizure on day 5.  
    mutate_at(vars(one_of(c("C57-M-Anubis"))), ~ ifelse(day >= 5, 0, .)) %>%
    #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
    mutate_at(vars(one_of(c("C57-F-Rae"))), ~ ifelse(day >= 2, 0, .)) %>%
    #T004: Hare only appears day 1. Not recovered, presumed dead. 
    mutate_at(vars(one_of(c("NYOB-M-Hare"))), ~ ifelse(day >= 2, 0, .)) %>%
    #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
    mutate_at(vars(one_of(c("NYOB-F-Isis"))), ~ ifelse(day >= 3, 0, .)) %>%
    #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
    mutate_at(vars(one_of(c("C57-F-Rose"))), ~ ifelse(day >= 10, 0, .))
  
  social_data[[aa]] <- df2
}

trial_net_stats_list <- list()
trial_vertex_stats_list <- list()
aa=1
for(aa in 1:7) {
  df <- social_data[[aa]] 
  geno <- print(unique(df$strain))
  trial <- print(unique(df$trial))
  
  # Get priority access scores for node attribute. 
  m_pas <- read_csv("data/priority_access_scores_males.csv")
  f_pas <- read_csv("data/priority_access_scores_females.csv")
  pas <- rbind(m_pas,f_pas)
  
  net_stats_list <- list()
  vertex_stats_list <- list()
  i=1
  for(i in 1:10) { # loop through days
    ## get PAS score for each day, not final day 10 PAS score. 
    pas1 <- pas %>% 
      filter(noon_day == i) %>% 
      dplyr::select(name, csum_daily_capture_penalty) %>% 
      dplyr::rename(priority_access_score = csum_daily_capture_penalty)
    
    gbi <- df %>%  # select mouse columns
      filter(day == i) %>% #comment out if you want the full 10 day network and run from here below
      dplyr::select(matches(c("*C57*", "*NYOB*")))   # choose your strain.
    
    colnames(gbi) <- gsub("C57-M-","",colnames(gbi))
    colnames(gbi) <- gsub("C57-F-","",colnames(gbi))
    colnames(gbi)<-gsub("NYOB-M-","",colnames(gbi))
    colnames(gbi)<-gsub("NYOB-F-","",colnames(gbi))
    ids <- colnames(gbi)
    
    ## create network
    undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION, # TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
                                data_format = "GBI",
                                association_index = "SRI") #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
    # write.csv(undir_matrix, "gbi_sn_SRI_T007_full.csv")
    
    g <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
    g <- simplify(g) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
    
    # create vertex stats
    V(g)$trial <- trial
    V(g)$genotype <- geno
    V(g)$day <- i 
    V(g)$sex = as.character(meta$sex[match(V(g)$name,meta$name)])
    V(g)$label <- V(g)$name # Add names to graph graphs (full 10 day plot)
    V(g)$family_group = as.character(meta$family_group[match(V(g)$name,meta$name)])
    V(g)$node_priority_access_score <- as.character(pas1$priority_access_score[match(V(g)$name,pas1$name)])
    V(g)$color = V(g)$sex #assign the "Sex" attribute as the vertex color
    if(geno == "C57") {
      V(g)$color = gsub("F","sienna1",V(g)$color) 
      V(g)$color = gsub("M","sienna",V(g)$color) 
    } else {
      V(g)$color = gsub("F","skyblue",V(g)$color) 
      V(g)$color = gsub("M","skyblue4",V(g)$color) 
    }
    V(g)$node_degree_centrality <- degree(g, mode="all") #node degree centrality, raw score. 
    V(g)$node_eigen_centrality <- eigen_centrality(g)$vector # measure of being connected to the well connected. not clear if this should be used or centr_eigen$vector
    V(g)$node_closeness_centrality <- closeness(g, mode = "all", normalized = TRUE, weights = NA) # measure of # steps required to access every other node. 
    V(g)$node_betweeness_centrality <- betweenness(g, directed = FALSE, weights = NA) # node BETWEENNESS = DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A VERTEX
    V(g)$node_edge_strength <- graph.strength(g) # node edge strength = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
    V(g)$node_page_rank <- page_rank(g)$vector
    V(g)$node_authority_score <- authority_score(g)$vector
    
    vertex_stats_list[[i]] <- igraph::as_data_frame(g, "vertices") # write daily vertex stats to list
    
    ## create net stats
    net_stats <- data.frame(matrix(ncol = 1,nrow = 1)) # create empty dataframe
    graph_centrality = centr_degree(g, mode = "all", loops = T, normalized = T) # https://igraph.org/r/doc/centr_degree.html
    net_stats$trial <- trial
    net_stats$genotype <- geno
    net_stats$day <- i
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
    net_stats$net_assortativity_family_group <- assortativity(g, as.factor(V(g)$family_group))
    
    net_stats[1] <- NULL #clean net_stats
    net_stats_list[[i]] <- net_stats ## write daily net stats to list
    
    ## Edge measures
    sort(edge_betweenness(g, directed = FALSE, weights = NA)) # EDGE BETWEENESS = DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A particular edge
    
    ## clustering dendrogram 
    ## https://www.rdocumentation.org/packages/igraph/versions/1.2.6/topics/plot_dendrogram
    ## height corresponds with strength of relationship 
    # plot_dendrogram(cluster_fast_greedy(g), use.modularity = F) 
    
    ## create basic graph
    svg(file=paste0("output/",trial,"_",geno,"_Day_", i, "_network.svg"), bg = "transparent")
    par(mar = c(0.4, 0.1, 2, 0.1))
    plot(g,
         layout = layout.fruchterman.reingold,
         vertex.size = scales::rescale(V(g)$node_edge_strength, to = c(5,25)), #rescale vertex strength to reasonable min/max size
         vertex.label= NA, # Remove vertex labels
         edge.width = scales::rescale(E(g)$weight, to = c(0.5,30)), # rescale edge weight to reasonable min/max size
         edge.color = "darkgray",
         edge.curved = 0.2,
    )
    title(paste0("Day ", i), cex.main=3) #
    dev.off()
    
    ## create advanced graph
    # svg(file=paste0("output/",trial,"_",geno,"_Day_", i, ".svg"), bg = "transparent")
    # par(mar = c(0.4, 0.1, 2, 0.1))
    # plot(g,
    #      layout = layout.fruchterman.reingold,
    #      # layout = layout_nicely,
    #      # vertex.size = V(g)$node_degree_centrality*2,
    #      vertex.size = scales::rescale(V(g)$node_edge_strength, to = c(5,25)), #rescale vertex strength to reasonable min/max size
    #      # vertex.label= V(g)$name, # include labels
    #      vertex.label= NA, # Remove vertex labels
    #      # vertex.label.font = 2,
    #      # vertex.label.color = "black",
    #      # vertex.label.cex = 1,
    #      # vertex.label.degree = 2,
    #      edge.width = scales::rescale(E(g)$weight, to = c(0.5,30)), # rescale edge weight to reasonable min/max size
    #      edge.color = "darkgray",
    #      edge.curved = 0.2,
    #      
    #      ### COMMUNITY CLUSTERING.
    #      # Decide if you actually want this. Maybe for full 10 day plot only.
    #      # mark.groups = cluster_fast_greedy(g),
    #      # mark.border = "darkgray"
    # )
    # title(paste0("Day ", i), cex.main=3) #
    # dev.off()
    # tkplot(g) # gui adjustments. cool!
    
    ## create circular graph
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
  
  # Get trial_net_stats_list
  trial_net_stats_list[[aa]] <- do.call("rbind", net_stats_list)
  
  ## Get vertex stats and COPY THE trial OUTPUT TO THE CLIPBOARD, paste directly into excel. 
  vertex_stats <- do.call("rbind", vertex_stats_list)
  
  trial_vertex_stats_list[[aa]] <- vertex_stats %>% 
    filter(!(name == "George")) %>% #T003: Anubis visually confirmed dead by seizure on day 5.  
    filter(!(name == "Anubis" & day >= 5)) %>% 
    filter(!(name == "Rae" & day >= 2)) %>%  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
    filter(!(name == "Hare" & day >= 2)) %>% #T004: Hare only appears day 1. Not recovered, presumed dead. 
    filter(!(name == "Isis" & day >= 3)) %>% #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
    filter(!(name == "Rose" & day >= 10)) %>% #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
    arrange(name, day)
}
alltrial_net_stats <- do.call("rbind", trial_net_stats_list)
alltrial_vertex_stats <- do.call("rbind", trial_vertex_stats_list)

write.csv(alltrial_net_stats, "data/ALLTRIAL_SNA_net_stats.csv")
write.csv(alltrial_vertex_stats, "data/ALLTRIAL_SNA_vertex_stats.csv")
