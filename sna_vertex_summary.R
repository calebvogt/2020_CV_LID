## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)


## import gbi_summary
gbi <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_summary.csv")

# Figure 4SX: Node measures and stats ---------------------------------------
# note, closeness not always a great measure on disconnected graphs, but that is in fact what we are trying to describe. 
# see here. net closeness limited to the largest sub-network. likely that node closeness is also limited. dividing by infinity, which is why we get such low values. 
# https://toreopsahl.com/2010/03/20/closeness-centrality-in-networks-with-disconnected-components/
# not one we are currently implementing. 
df <- read_excel("Data_Stats_v6.xlsx", sheet = "Fig4_node_data")

### REPORTED NODE MEASURES
# node_edge_strength

### UNREPORTED NODE MEASURES
# node_eigen_centrality
# node_closeness_centrality
# node_betweeness_centrality
# node_authority_score

## Choose your node measure and make your graphs.
df1 <- df %>% 
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(node_edge_strength), 
            sd = sd(node_edge_strength), 
            count = n(), 
            sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color = strain_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab("Day") +
    ylab("node_edge_strength") +
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
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") 

# STATS
df$trial <- as.factor(df$trial)
df$sex <- as.factor(df$sex)
df$strain <- as.factor(df$strain)
df$day <- as.numeric(df$day)
mod1 = lmer(node_edge_strength ~ strain*sex*log(day) + (1|trial) + (1+log(day)|name), data = df) 
AIC(mod1)
anova(mod1)
summary(mod1)
write.table(summary(mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
emmeans(mod1, pairwise ~ strain*sex*log(day))
ggpredict(mod1, terms = c("strain", "sex", "day"), type = "fe") %>% 
  plot()

# daily models: change for days 1 through 10. 
# i = 4
for(i in 1:10){
  d = lmer(closeness ~ strain*sex + (1|trial), data = subset(df, day == i)) 
  #get daily contrast estimates
  contrasts <- emmeans(d, pairwise ~ strain*sex)
  print(i)
  print(contrasts$contrasts)
  # anova(d)
  # summary(d)
}



