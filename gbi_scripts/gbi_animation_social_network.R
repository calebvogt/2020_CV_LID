## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

# load  & clean data ---------------------------------------------------------------
filenames <- list.files("data/", pattern = "*MOVEBOUT_GBI.csv", full.names = T)
social_data = lapply(filenames, fread) ## READ IN ALL FILES
meta <- read_excel("data/metadata.xlsx", skip=1)

# clean social data for triaged mice from social interaction bouts. 
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  
  df2 <- df %>% 
    select(!(V1)) %>% 
    filter(day %in% 1:10) ## filter range of days for analysis
  
  print(aa)
  social_data[[aa]] <- df2
}


# Create sliding window dynamic net animation -----------------------------------------
## http://estebanmoro.org/post/2015-12-21-temporal-networks-with-r-and-igraph-updated/
## https://kateto.net/network-visualization

df <- social_data[[7]] # goes in order of trials. 
print(unique(df$trial))
df$field_time_start <- as.POSIXct(df$field_time_start, format="%m/%d/%Y %H:%M")
df$field_time_stop <- as.POSIXct(df$field_time_stop, format="%m/%d/%Y %H:%M")
df <- df[order(field_time_start)]
df$time_hours <- floor(as.numeric(difftime(df$field_time_start, df$field_time_start[[1]], units = "hours"))) + 1 # create hours column in df
df <- df %>% # move time columns. filter out rows with less than 2 mice. 
  select(time_hours, everything()) %>% 
  filter(mf_sum > 1)

el_list <- list()
i=1
for(i in min(df$time_hours):(max(df$time_hours))) { # loop through time intervals. Weirdly, sometimes throws an error, but just keep running until it works, usually works after 2-3 times
  gbi <- df %>%  # select mouse columns
    filter(between(time_hours, i-24,i)) %>% ## this results in a true sliding window of 24 hours over the entire trial period
    dplyr::select(matches(c("*C57*", "NYOB*")))   # choose your strain from the column information. 
  
  colnames(gbi) <- gsub("C57-M-","",colnames(gbi))
  colnames(gbi) <- gsub("C57-F-","",colnames(gbi))
  colnames(gbi)<-gsub("NYOB-M-","",colnames(gbi))
  colnames(gbi)<-gsub("NYOB-F-","",colnames(gbi))
  ids <- colnames(gbi)
  
  undir_matrix <- get_network(association_data = gbi,    # ASNIPE FUNCTION, # TURN GBI DATA INTO UNDIRECTED WEIGHTED ADJACENCY MATRIX FOR USE IN IGRAPH. 
                              data_format = "GBI",
                              association_index = "SRI") #CHOOSE SIMPLE RATIO INDEX OR SAMPLING PERIODS ARRAY
  g <- graph.adjacency(undir_matrix, mode="undirected", weighted = TRUE, diag = FALSE) # CREATE IGRAPH OBJECT FROM ADJACENCY MATRIX. GRAPHS SHOULD BE UNDIRECTED. 
  g <- simplify(g) # SIMPLIFY IGRAPH OBJECT. REMOVES MULTIPLE EDGES AND LOOP EDGES.
  print(i)
  
  el <- as.data.frame(cbind(get.edgelist(g), E(g)$weight)) ## turn adjacency matrix into edge list with weights
  el$time <- i
  colnames(el) <- c('id1', 'id2','weight','time')
  el_list[[i]] <- el 
}
edges <- do.call("rbind", el_list)
g <- graph.data.frame(edges,directed=F) ## generate the full graph to adjust settings
# plot(g) ## plot the full graph, but note that this is not what will ultimately be used


## create vertex attributes
meta <- meta %>% mutate(node_color = paste(strain, sex,sep="-")) ## create node_color groups
meta <- meta %>% mutate(node_code = paste(sex, code,sep="-")) ## create node_color groups
V(g)$sex <- as.character(meta$sex[match(V(g)$name,meta$name)])
V(g)$code <- as.character(meta$code[match(V(g)$name,meta$name)])
V(g)$node_code <- as.character(meta$node_code[match(V(g)$name,meta$name)])
V(g)$node_color = as.character(meta$node_color[match(V(g)$name,meta$name)])
V(g)$color = V(g)$node_color ## Have to set the founder and invader attributes by color and not shape, due to issue with the log degree, only good way to make nodes disappear. 
V(g)$color = gsub("C57-F","sienna1",V(g)$color) 
V(g)$color = gsub("C57-M","sienna",V(g)$color) 
V(g)$color = gsub("NYOB-F","skyblue",V(g)$color) 
V(g)$color = gsub("NYOB-M","skyblue4",V(g)$color) 

## create 24h sliding window networks
g2 <- delete_edges(g,which(E(g)$time > 1)) ## get the layout for the first time step
layout.old <- norm_coords(layout.graphopt(g2), xmin = -1, xmax = 1, ymin = -1, ymax = 1) ## create normal layout
dt <- 1 ## set animation time interval. Check this and the delete_edges argument, related

png(file=paste0("output/animations/pngs/", "example%03d.png"), width=1920,height=1080) ## set up png
i=1
for(i in seq(min(E(g)$time),max(E(g)$time),dt)){ ##
  print(i)
  g2 <- delete_edges(g, which(E(g)$time !=i)) #deletes edges not in the relevant time period. note that time X represents that hour and the preceding 24h worth of data
  # g2 <- delete_edges(g, which(E(g)$time > i | E(g)$time < i)) #deletes edges not in the relevant time period. note that time X represents that hour and the preceding 24h worth of data
  E(g2)$time #check that it worked
  layout.new <- layout_with_fr(g2, ## with the new graph, we update the layout a little bit
                               coords=layout.old, 
                               niter=10, ## number of iterations of energy minimization
                               start.temp=0.05, ## maximum amount of movement along one axis 
                               grid="nogrid") 
  # valid_nodes <- which(degree(g2) != 0)
  # g3 <- induced.subgraph(g2, valid_nodes) ## create graph 
  
  plot(g2,
       layout=layout.new, ## provides coordinates 
       vertex.label= V(g2)$node_code, # include labels
       # vertex.size=10+3*log(degree(g2)), ## log scale means non-present nodes disappear due to infinite values not displayed. Doesn't appear to work with non-circular shapes.
       vertex.size = scales::rescale(degree(g2), to = c(5,30)), ## rescale vertex strength to reasonable min/max size
       edge.width = scales::rescale(as.numeric(E(g2)$weight), to = c(0.5,15)), # rescale edge weight to reasonable min/max size
       edge.color = "darkgray",
       edge.curved = 0.2,
       # asp=9/16, ## change the aspect ratio. defaults to circle
       # margin=-0.15, ## set the margin. 
  )
  # title(paste0("Night ", ceiling(i/24)), cex.main=3) ## get night 
  title(paste0("Night ", round(i/24+1,1)), cex.main=4) ## get night + decimal
  layout.old <- layout.new
}
dev.off()

## DELETE LAST EMPTY PNG IN THE FOLDER PRIOR TO MAKING VIDEO. 

## create mp4 using terminal command 
## open cmd window in png folder
## -r sets frames/pngs shown per second. 1 png per hour. 

ffmpeg -r 15 -i example%03d.png -b:v 20M output.mp4 
ffmpeg -r 15 -i example%03d.png -b:v 20M output.gif 

