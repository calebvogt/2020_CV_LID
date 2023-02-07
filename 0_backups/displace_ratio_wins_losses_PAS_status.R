## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(lme4)
library(lmerTest)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

## Load displace data
displace <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT_GBI_displace.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# displace_ratio_win_loss_PAS_status ---------------------------------------------
wins <- displace %>% 
  group_by(trial,winner, winner_PAS_status, day) %>% 
  tally() %>% 
  rename(win_count = n, name = winner)

losses <- displace %>% 
  group_by(trial, loser, loser_PAS_status, day) %>% 
  tally() %>% 
  rename(loss_count = n, name = loser)

df2 <- merge(wins, losses[,c("name", "day", "loss_count")], by = c("name", "day"))

df3 <- df2 %>% 
  mutate(WL_ratio = win_count / loss_count) %>% 
  rename(PAS_status = winner_PAS_status)

## Boxplot 
(p <- ggplot(df3, aes(x=as.factor(day), y=WL_ratio, fill = as.factor(PAS_status))) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    xlab("Day") +
    ylab("Win-Loss ratio") +
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
# ggsave(p, filename = "output/displace_ratio_win_loss_PAS_status_boxplot.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/displace_ratio_win_loss_PAS_status_boxplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


## lineplot
df4 <- df3 %>% 
  group_by(PAS_status, day) %>% 
  summarise(mean = mean(WL_ratio), sd = sd(WL_ratio),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df4, aes(x = day, y = mean, color = PAS_status)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
    # values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Win-Loss ratio") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = "right")
)
# ggsave(p, filename = "output/displace_ratio_win_loss_PAS_status_lineplot.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/displace_ratio_win_loss_PAS_status_plot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m1 = lmer(WL_ratio ~ PAS_status + (1|name), data = df3) 
summary(m1)

# write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

