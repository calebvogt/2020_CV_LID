## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)


df <- read.csv("data/ALLTRIAL_SNA_vertex_stats.csv")
df$trial <- as.factor(df$trial)
df$genotype <- as.factor(df$genotype)
df$day <- as.numeric(df$day)
df$sex <- as.factor(df$sex)

# sna_vertex_degree_centrality.png ----------------------------------------
df1 <- df %>% 
  mutate(genotype_sex = paste0(genotype, "-", sex)) %>% 
  group_by(genotype_sex, day) %>%
  summarise(mean = mean(node_degree_centrality), sd = sd(node_degree_centrality), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color = genotype_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "Outbred-F", "Outbred-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Node degree centrality") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.10,0.9))
          # legend.position = "left")
          legend.position = "")
)
# ggsave(p, filename = "output/sna_vertex_degree_centrality.png", device = "png", bg = "white")
ggsave(p, filename = "output/sna_vertex_degree_centrality.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# STATS
m1 = lmer(node_degree_centrality ~ genotype*sex*log(day) + (1|trial), data = df) 
summary(m1)

m2 = lmer(node_degree_centrality ~ genotype*sex*log(day) + (1|trial) + (log(day)|name), data = df) 
summary(m2)

m3 = lmer(node_degree_centrality ~ genotype*sex*day + (1|trial) + (day|name), data = df) 
summary(m3)

AIC(m1,m2,m3)
write.table(summary(m2)$coef, "clipboard-16384", sep="\t", row.names=TRUE, col.names = TRUE)

# # daily models: change for days 1 through 10. 
# # i = 4
# for(i in 1:10){
#   d = lmer(node_degree_centrality ~ genotype*sex + (1|trial), data = subset(df, day == i)) 
#   #get daily contrast estimates
#   contrasts <- emmeans(d, pairwise ~ genotype*sex)
#   print(i)
#   print(contrasts$contrasts)
#   # anova(d)
#   # summary(d)
# }


