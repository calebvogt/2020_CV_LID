## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(lme4)
library(lmerTest)

metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

## Load displace data
displace <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT_GBI_displace.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# displace_duration_home_wins ---------------------------------------
df2 <- displace %>% 
  group_by(strain, winner, winner_loc) %>% 
  filter(winner_loc == "home")

## Boxplot 
(p <- ggplot(df2, aes(x=as.factor(day), y=duration_s, fill = as.factor(strain))) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    xlab("Day") +
    ylab("Duration (s) of home win displacment events") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "right") 
)
# ggsave(p, filename = "output/displace_duration_home_wins.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/displace_duration_home_wins.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# displace_duration_home_losses ---------------------------------------
df2 <- displace %>% 
  group_by(strain, loser, loser_loc) %>% 
  filter(loser_loc == "home")

## Boxplot 
(p <- ggplot(df2, aes(x=as.factor(day), y=duration_s, fill = as.factor(strain))) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    xlab("Day") +
    ylab("Duration (s) of home losses displaces") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "right") 
)
# ggsave(p, filename = "output/displace_duration_home_losses.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/displace_duration_home_losses.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


