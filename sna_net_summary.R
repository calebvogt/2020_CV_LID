## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)

df <- read.csv("data/ALLTRIAL_SNA_net_stats.csv")
df$trial <- as.factor(df$trial)
df$genotype <- as.factor(df$genotype)
df$day <- as.numeric(df$day)

# sna_net_num_graph_components.png -------------------------------------------
df1 <- df %>% 
  group_by(genotype, day) %>%
  summarise(mean = mean(net_components_num), sd = sd(net_components_num), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color = genotype)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"), values=c("goldenrod1", "goldenrod4")) +
    theme_classic() +
    xlab("Day") +
    ylab("# network components") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.2,0.85))
          legend.position = "")
)
# ggsave(p, filename = "output/sna_net_num_graph_components.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_net_num_graph_components.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## STATS 
m1 = lmer(net_components_num ~ genotype*log(day) + (1|trial), data = df) 
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# sna_net_edge_density.png ------------------------------------------------
df1 <- df %>% 
  group_by(genotype, day) %>%
  summarise(mean = mean(net_edge_density), sd = sd(net_edge_density), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color = genotype)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"), values=c("goldenrod1", "goldenrod4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Network edge density") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.2,0.85))
          legend.position = "")
)
# ggsave(p, filename = "output/sna_net_edge_density.png", device = "png", bg = "white")
ggsave(p, filename = "output/sna_net_edge_density.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## STATS 
m1 = lmer(net_edge_density ~ genotype*log(day) + (1|trial), data = df) 
summary(m1)
AIC(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# sna_net_eigen_centrality.png ------------------------------------------------
df1 <- df %>% 
  group_by(genotype, day) %>%
  summarise(mean = mean(net_eigen_centrality), sd = sd(net_eigen_centrality), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color = genotype)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"), values=c("goldenrod1", "goldenrod4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Network eigenvector centrality") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.2,0.85))
          legend.position = "")
)
# ggsave(p, filename = "output/sna_net_eigen_centrality.png", device = "png", bg = "white")
ggsave(p, filename = "output/sna_net_eigen_centrality.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## STATS 
m1 = lmer(net_eigen_centrality ~ genotype*log(day) + (1|trial), data = df) 
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# # daily models: change for days 1 through 10. 
# i = 4
# for(i in 1:10){
#   d = lm(net_components_num ~ genotype, data = subset(df, day == i)) 
#   #get daily contrast estimates
#   contrasts <- emmeans(d, pairwise ~ genotype)
#   print(i)
#   print(contrasts$contrasts)
#   # anova(d)
#   # summary(d)
# }

### REPORTED SIG RESULTS
# net_edge_density #sig
# net_components_num # sig
# net_eigen_centrality # sig  #https://igraph.org/r/doc/centr_eigen.html

### UNREPORTED SIG RESULTS
# net_centrality # sig # https://igraph.org/r/doc/centr_degree.html
# net_mean_dist # sig
# net_modularity_infomap_group_n #sig
# net_modularity_fast_greedy_group_n #sig
# net_transitivity # NS, except for day
# net_modularity_infomap # NS
# net_modularity_fast_greedy # NS


