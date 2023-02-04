## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(data.table)
library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(ggpubr)

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# Number of zone visits and time in zone correlation -----------------------
# number of zone visits
df <- move_data %>% 
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(trial, strain, sex, name, zone) %>% 
  tally()

# duration of time spent in zone
df2 <- move_data %>% 
  mutate(duration_min = duration_s / 60)  %>% 
  group_by(trial, strain, sex, name, zone) %>% 
  tally(sum(duration_min))

df3 <- as.data.frame(cbind(df, df2$n))
df3$strain_sex <- paste0(df3$strain,"-",df3$sex)
colnames(df3) <- c("trial", "strain", "sex", "name", "zone", "num_visits", "total_time_min", "strain_sex")


(p <- ggplot(data = df3, aes(x = num_visits, y = total_time_min, color = strain_sex)) +
    geom_point() +
    scale_color_brewer(palette = "PuOr") + 
    geom_smooth(method="lm") + 
    xlab("Number of zone visits") +
    ylab("Time spent in zone (min)") +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    stat_cor(aes(color = strain_sex), label.x = 4, method = "pearson", p.accuracy = 0.001) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
)
# ggsave(p, filename = "output/move_zone_visit_time_correlation.png", device = "png", bg = "white")
# ggsave(p, filename = "output/move_zone_visit_time_correlation.svg", device = "svg", width=4.5, height=3, bg = "transparent")



