# 6_analyze_ALLTRIAL_MOVEBOUT_GBI
## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# load data and packages --------------------------------------------------
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

wd <- setwd("C:/Users/caleb/Box/0_CV_Shared_Analyses/7_LID_2020/RFID_analysis_v10")
output_fp <- paste("C:/Users/caleb/Desktop")
filenames <- list.files(wd, pattern = "*MOVEBOUT_GBI.csv")
myfiles = lapply(filenames, fread) ## READ IN ALL FILES
metadata <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
meta_short <- metadata %>% 
  select(trial, paddock, strain, sex, name, code, family_group)

# FIG4A-D_SXA-X_C57/Outbred social network structure over time -----------------------------------------
## TO DO!!!
## 1. you need to remove your outliers from the data set before you construct the matrix. This is going to alter a variety of your measures, so this is high priority
## 2. 

df <- as.data.frame(fread("T007_MOVEBOUT_GBI.csv", stringsAsFactors = TRUE)) ### CHANGE THIS
df <- subset(df, select = -c(V1))
net_stats_list <- list()
vertex_stats_list <- list()
node_stats_list <- list()
i=1
for(i in 1:10) { # loop through days
  gbi <- df %>%  # select mouse columns
    filter(day == i) %>% #comment out if you want the full network.
    # select(matches("*C57*"))
  select(matches("*NYOB*"))
  # select(matches("*-M-"))
  # select(matches("*-F-"))
  
  # colnames(gbi) <- gsub("C57-M-","",colnames(gbi)) #remove strain-sex info from colnames
  # colnames(gbi) <- gsub("C57-F-","",colnames(gbi))
  colnames(gbi)<-gsub("NYOB-M-","",colnames(gbi))
  colnames(gbi)<-gsub("NYOB-F-","",colnames(gbi))
  ids <- colnames(gbi)
  
  ## create network
  undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION, # TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
                              data_format = "GBI",
                              association_index = "SRI") #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
  g <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
  g <- simplify(g) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES. 
  
  ## Network  Measures
  net_stats <- data.frame(matrix(ncol = 1,nrow = 1)) # create empty dataframe
  net_stats$mean_distance <- mean_distance(g, directed = FALSE) # mean distance = avg. # of edges between any two nodes in network. 
  graph_centrality = centr_degree(g, mode = "all") # GET GRAPH CENTRALITY MEASURES AT NODE AND NETWORK LEVEL
  net_stats$graph_centrality <- graph_centrality$centralization
  net_stats$edge_density <- edge_density(g, loops = FALSE) # EDGE DENSITY = RATIO OF THE NUMBER OF EDGES AND THE NUMBER OF ALL POSSIBLE EDGES. Single Number
  net_stats$transitivity <- transitivity(g) # aka clustering coefficient, probability that adjacent nodes of a network are connected. 
    cluster <- cluster_infomap(g) #infomap method attempts to map the flow of information in a network, and the different clusters in which information may get remain for longer periods. Similar to walktrap, but not necessarily maximizing modularity, but rather the so-called “map equation”.
  net_stats$infomap_modularity <- modularity(cluster)
  # net_stats$reciprocity <- reciprocity(g) # for directed graphs only: propensity of each edge to be a mutual edge. 
  #clean nets_stats
  net_stats[1] <- NULL
  net_stats_list[[i]] <- net_stats
  
  ## Vertex Measures
  V(g)$day <- i 
  V(g)$degree <- degree(g, mode="all") #VERTEX DEGREE
  V(g)$strength <- graph.strength(g) # assign vertex strengths. STRENGTH = SUM OF ALL EDGE WEIGHTS FOR A SINGLE NODE (BETWEEN 0-1)
  V(g)$closeness <- closeness(g, mode = "all", normalized = TRUE, weights = NA) # CLOSENESS = MEASURES HOW MANY STEPS IT TAKES TO REACH EVERY OTHER VERTEX FROM A GIVEN VERTEX. # RUN ON GRAPHS THAT ARE WELL CONNECTED.
  V(g)$betweeness <- betweenness(g, directed = FALSE, weights = NA) # VERTEX BETWEENNESS = DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A VERTEX
  V(g)$eigen_centrality <- eigen_centrality(g)$vector
  V(g)$page_rank <- page_rank(g)$vector
  V(g)$authority_score <- authority_score(g)$vector
  V(g)$sex = as.character(meta_short$sex[match(V(g)$name,meta_short$name)])
  # V(g)$label <- V(g)$name # remove if you dont want names on graphs
  V(g)$color = V(g)$sex #assign the "Sex" attribute as the vertex color
  V(g)$color = gsub("F","red",V(g)$color) #Females will be red
  # V(g)$color = gsub("M","lightblue",V(g)$color) #Males will be blue
  V(g)$color = gsub("M","blue",V(g)$color) #Males will be blue
  vertex_stats <- igraph::as_data_frame(g, "vertices")
  vertex_stats_list[[i]] <- vertex_stats
  
  ## Edge measures
  # sort(edge_betweenness(g, directed = FALSE, weights = NA)) # EDGE BETWEENESS = DEFINED BY THE NUMBER OF GEODESICS (SHORTEST PATHS) GOING THROUGH A particular edge
  
  ## PLOT 1: standard network
  svg(file=paste0(output_fp,"/","P1_Day_", i, ".svg"), bg = "transparent") ## M+F plot
  par(mar = c(0.4, 0.1, 2, 0.1))
  plot(g, 
       layout = layout.fruchterman.reingold,
       # layout = layout_nicely,
       # vertex.size = V(g)$degree*2,
       vertex.size = scales::rescale(V(g)$graph.strength, to = c(5,30)), #rescale vertex strength to reasonable min/max size
       # vertex.label= V(g)$name, # include labels
       vertex.label= NA, # Remove vertex labels
       # vertex.label.font = 2,
       # vertex.label.color = "black",
       # vertex.label.cex = 1,
       # vertex.label.degree = 2,
       edge.width = E(g)$weight*100, #perhaps this should be scaled as well. 
       edge.color = "darkgray",
       edge.curved = 0.2,
       
       ### COMMUNITY CLUSTERING. Only works if all vertices are connected. 
       mark.groups = cluster_infomap(g),
       mark.border = "darkgray"
  )
  title(paste0("Day ", i), cex.main=3) #
  dev.off()
  # tkplot(g) # gui adjustments. cool!
  

  ## PLOT 2: Circular Network
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

# Get net Stats
net_stats <- do.call("rbind", net_stats_list)
write.table(net_stats, "clipboard", sep="\t", row.names=FALSE, col.names = FALSE) # COPY THE OUTPUT TO THE CLIPBOARD, paste directly into excel. 

#net_stats plot
df <- read_excel(paste(wd,"1_paper_documents", "Figs_Data_v2.xlsx", sep = "/"), sheet = "Fig4C")
df1 <- df %>% 
  group_by(strain, day) %>%
  summarise(mean = mean(centrality), # 
            sd = sd(centrality), 
            count = n(), 
            sem = (sd/(sqrt(count))))


(p <- ggplot(df1, aes(x=day, y=mean, color = strain)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(1,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
    # scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       # values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab("Day") +
    ylab("Network centrality") +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "top")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent") 

## Get vertex stats
vertex_stats <- do.call("rbind", vertex_stats_list)
vertex_stats <- merge(vertex_stats, meta_short, by.x = "name", by.y = "name")
vertex_stats <- vertex_stats %>% 
  dplyr::rename(sex = sex.x) %>% 
  select(trial, strain, sex, name, code, day, degree, strength, betweeness, closeness, eigen_centrality, page_rank, authority_score) %>% 
  arrange(name, day)
write.table(vertex_stats, "clipboard", sep="\t", row.names=FALSE, col.names = F) # COPY THE OUTPUT TO THE CLIPBOARD, paste directly into excel. 

## vertex_stats plot
df <- read_excel(paste(wd,"1_paper_documents", "Figs_Data_v2.xlsx", sep = "/"), sheet = "Fig4E")

#data cleaning (performed once! excel file updated)
clean <- df %>% 
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
  filter(!(name == "Rose" & day >= 10))  
write.table(clean, "clipboard-16384", sep="\t", row.names=FALSE, col.names = T) # COPY THE OUTPUT TO THE CLIPBOARD, paste directly into excel. 

df1 <- df %>% 
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(closeness), # degree, strength, betweeness, closeness, eigen_centrality, page_rank, authority_score
            sd = sd(closeness), 
            count = n(), 
            sem = (sd/(sqrt(count))))


(p <- ggplot(df1, aes(x=day, y=mean, color = strain_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(1,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab("Day") +
    ylab("Vertex closeness") +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          # legend.position = "right")
          legend.position = "")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent") 




# FIG4C_Cumulative sum of unique mice encountered ------------
## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
males <- meta_short %>% 
  filter(sex == "M") %>% 
  select(name) %>% 
  filter(!is.na(name))
male_list <- dplyr::pull(males, name)

females <- meta_short %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  select(name) %>% 
  filter(!is.na(name))
female_list <- dplyr::pull(females, name)

trial_stats <- list()
aa = 1
for(aa in 1:length(myfiles)){
  df <- myfiles[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  df <- subset(df, select = -c(V1))
  colnames(df)<-gsub("C57-M-","",colnames(df))
  colnames(df)<-gsub("C57-F-","",colnames(df))
  colnames(df)<-gsub("NYOB-M-","",colnames(df))
  colnames(df)<-gsub("NYOB-F-","",colnames(df))
  
  col_ids <- colnames(df[,10:ncol(df)]) ## get mouse column names starting at col 10
  bb = col_ids[1]
  all_mouse_list <- list()
  first_flag = 1
  for(bb in col_ids[1:length(col_ids)]) {
    df2 <- df %>% 
      filter((!!as.symbol(bb)) == 1) %>% # pull all rows where bb mouse is present. 
      mutate(name = bb) %>% 
      relocate(name)
    
    non_self_ids <- col_ids[!col_ids %in% bb] # remove current mouse from the next loop to compare to other animals
    non_self_ids[1]
    novel_mouse_rows <- list()
    second_flag = 1
    for(i in non_self_ids[1:length(non_self_ids)]) {
      df3 <- df2 %>% 
        filter((!!as.symbol(i)) == 1) %>% 
        mutate(novel_mouse_met = i)
      
      novel_mouse_rows[[second_flag]] <- df3[1,] #save first observed meeting of the focal and novel mouse to list
      second_flag = second_flag + 1
    }
    all_mouse_list[[first_flag]] <- do.call("rbind", novel_mouse_rows)
    first_flag <- first_flag + 1
  }
  df4 <- do.call("rbind", all_mouse_list)
  
  #remove na rows which are introduced when a mouse does not ever meet a particular other mouse. 
  df4[rowSums(is.na(df4)) > 0,]
  df5 <- df4[complete.cases(df4), ] 
  
  ## ADD RELEVANT METADATA INFORMATION. 
  df6 <- merge(df5, meta_short, by.x = "name", by.y = "name")
  df7 <- df6 %>%
    dplyr::rename(trial = trial.x) %>% 
    select(trial, paddock, strain, sex, name, code, novel_mouse_met, day, zone, field_time_start, field_time_stop, m_sum, f_sum, mf_sum, duration_s)
  
  trial_stats[[aa]] <- df7
}
df8 <- do.call("rbind", trial_stats)
df9 <- df8[with(df8, order(name, day)),]

# df9 has a dataframe ordered by first meeting time with each mouse in the paddock .

df10 <- df9 %>% 
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  group_by(strain_sex, name, day) %>% 
  tally() %>% # 
  complete(name, day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse, adds NAs. 
  replace(is.na(.), 0) %>% 
  mutate(csum = cumsum(n)) %>% #get cumulative # of novel mice met
  arrange(name, day) %>% 
  fill(csum) %>% ## fill cumulative sum data from last observed day
  dplyr::rename(novel_indivs_met = n, csum_novel_indivs_met = csum)

#data cleaning
df11 <- df10 %>% 
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
  filter(!(name == "Rose" & day >= 10))  

# output csv
write.csv(df11, paste(wd,"1_paper_documents", "Figure4C_data.csv", sep = '/'), row.names = FALSE)

#plot
df12 <- df11 %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(csum_novel_indivs_met), 
            sd = sd(csum_novel_indivs_met), 
            count = n(), 
            sem = (sd/(sqrt(count))))


(p <- ggplot(df12, aes(x=day, y=mean, color = strain_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(1,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab("Day") +
    ylab("Cumulative Unique Mice Met") +
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




