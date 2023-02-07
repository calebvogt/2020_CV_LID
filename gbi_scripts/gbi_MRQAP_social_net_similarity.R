## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
social_data = lapply(filenames, fread) ## READ IN ALL FILES


# Figure 4f: MRQAP.DSP Daily Network Prediction of Day 10 Network -------------------
# get trial social GBI data. 
library(asnipe)
df <- social_data[[1]] # goes in order of trials. 
print(unique(df$trial))
net_stats_list <- list()
vertex_stats_list <- list()
node_stats_list <- list()

# BONUS: Do Day 1 networks predict Day 10 networks?
adj_list <- list()
aa = 1
for(aa in 1:10){
  gbi2 <- df %>%  # select mouse columns
    filter(day == aa) %>%
    dplyr::select(matches(c("*C57*","*NYOB*")))   # choose your strain.
  # convert gbi to ibg to adjacency matrix.
  ibg <- t(gbi2)
  adj_list[[aa]] <- ibg %*% t(ibg)
}

for(bb in 1:10) {
  print(bb)
  print(mrqap.dsp(adj_list[[10]]~adj_list[[bb]], directed = "undirected", test.statistic = "t-value", randomisations = 10))
}
# manually copy and paste outputs into excel. 

df <- read_excel("Figure_Data.xlsx", sheet = "Fig4f")


## Figure 4f
(p <- ggplot(df, aes(x=day, y=adjusted_r_squared, group = trial, color = strain, fill = strain, shape = significant)) + 
    geom_line(size = 0.7) + 
    geom_point(size = 3) +
    # geom_point(aes(color = strain), size = 4, shape = 21) +
    scale_color_manual(breaks = c("C57", "Outbred"), values=c("goldenrod1", "goldenrod4")) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    theme_classic() +
    xlab("Day") +
    # ylab("net_mean_dist") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.2,0.85))
          legend.position = "right")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=5, height=2.5, bg = "transparent") 




