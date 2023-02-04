## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)

# http://estebanmoro.org/post/2015-12-21-temporal-networks-with-r-and-igraph-updated/
# https://kateto.net/network-visualization

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)
social_data = lapply(filenames, fread) ## READ IN ALL FILES


# TESTING!!! social_data: Hourly Social network evolution animation -----------------------------------------
# get trial social GBI data. 
df <- social_data[[4]] # goes in order of trials. 
print(unique(df$trial))
el_list <- list()
df$field_time_start <- as.POSIXct(df$field_time_start, format="%m/%d/%Y %H:%M")
df$field_time_stop <- as.POSIXct(df$field_time_stop, format="%m/%d/%Y %H:%M")
df <- df[order(field_time_start)]

# create hours column in df
df$time_hours <- floor(as.numeric(difftime(df$field_time_start, df$field_time_start[[1]], units = "hours"))) + 1

# move time columns. filter out rows with less than 2 mice. 
df <- df %>% 
  select(time_hours, everything()) 
# filter(mf_sum > 1)

# Create sampling period
df1 <- df %>% 
  mutate(samp_period = cut(time_hours, seq(0,length(df$time_hours), 2))) # set time interval

df2 <- transform(df1,time_step = as.numeric(factor(samp_period)))

df <- df2 %>% 
  select(time_hours, time_step, everything())

i=118
for(i in 1:(max(df$time_step))) { # loop through time intervals
  gbi <- df %>%  # select mouse columns
    filter(time_step == i) %>% #comment out if you want the full 10 day network and run from here below
    dplyr::select(matches(c("*C57*", "*NYOB*")))   # choose your strain.
  
  # remove strain-sex info from colnames of strain you are working with. 
  colnames(gbi) <- gsub("C57-M-","",colnames(gbi))
  colnames(gbi) <- gsub("C57-F-","",colnames(gbi))
  colnames(gbi)<-gsub("NYOB-M-","",colnames(gbi))
  colnames(gbi)<-gsub("NYOB-F-","",colnames(gbi))
  ids <- colnames(gbi)
  
  ## create network
  undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION, # TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
                              data_format = "GBI",
                              association_index = "SRI") #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
  g <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
  g <- simplify(g) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES.
  print(i)
  #turn adjacency matrix into edge list
  el <- as.data.frame(get.edgelist(g))
  el$time <- i
  colnames(el) <- c('id1', 'id2','time')
  el_list[[i]] <- el 
}
edges <- do.call("rbind", el_list)

#generate the full graph, plot the full graph
g <- graph.data.frame(edges,directed=F)
plot(g)

# set up vertex attributes
V(g)$sex = as.character(meta_short$sex[match(V(g)$name,meta_short$name)])
V(g)$label <- V(g)$name # Add names to graph graphs (full 10 day plot)
V(g)$family_group = as.character(meta_short$family_group[match(V(g)$name,meta_short$name)])
V(g)$color = V(g)$sex #assign the "Sex" attribute as the vertex color
V(g)$color = gsub("F","sienna1",V(g)$color) #Females will be red
V(g)$color = gsub("M","sienna",V(g)$color) #Males will be blue
layout.old <- norm_coords(layout.fruchterman.reingold(gt), xmin = -1, xmax = 1, ymin = -1, ymax = 1)

#total time of networks
total_time <- max(E(g)$time)

#This is the time interval for the animation. In this case is taken to be 1/10
#of the time (i.e. 10 snapshots) between adding two consecutive nodes
dt <- 0.1

png(file=paste0(output_fp, "/", "example%03d.png"), width=1920,height=1080)
time=56
for(time in seq(1,total_time,dt)){
  #remove edges which are not present
  # gt <- delete_edges(g, which(E(g)$time > time)) #deletes edges forward in time, but not backwards!
  gt <- delete_edges(g, which(E(g)$time > time | E(g)$time < floor(time)-12)) #deletes edges not in the relevant time period. Keep edges from the last 24 hours. 
  E(gt)$time #check that it worked
  
  #with the new graph, we update the layout a little bit
  layout.new <- layout_with_fr(gt, coords=layout.old, niter=10, start.temp=0.05, grid="nogrid")
  
  #plot the new graph
  plot(gt,
       layout=layout.new,
       # vertex.size=1+2*log(degree(gt)),
       vertex.size = scales::rescale(V(g)$node_edge_strength, to = c(5,25)), #rescale vertex strength to reasonable min/max size
       vertex.label= V(g)$name, # include labels
       vertex.frame.color=V(g)$color,
       # edge.width=1.5, 
       edge.width = scales::rescale(E(g)$weight, to = c(0.5,30)), # rescale edge weight to reasonable min/max size
       edge.color = "darkgray",
       edge.curved = 0.2,
       asp=9/16,margin=-0.15)
  layout.old <- layout.new
}
dev.off()

## terminal command
# ffmpeg -r 20 -i example%03d.png -b:v 20M output.mp4




