## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

move_data <- merge(move_data, metadata[,c("name", "zone_drop")], by = "name", all = FALSE) # bring in metadata

df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, strain_sex, strain,sex,name,zone_drop, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) %>% 
  group_by(name) %>% 
  slice_max(percent_time, n = 2) %>% 
  mutate(rank_order = rank(desc(percent_time))) %>% 
  complete(rank_order = 1:2, fill = list(n = 0, percent_time = 0)) %>% 
  mutate(rank_order = as.factor(rank_order))

df1 <- df %>% 
  dplyr::select(trial, name, strain, sex, rank_order, zone, n, percent_time) %>% 
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  filter(rank_order == 1) %>% 
  dplyr::rename(duration_min = n)

# write.table(df1, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
# df <- read_excel("Figure_Data.xlsx", sheet = "Fig2d")

table(round(df1$percent_time, 6), df1$strain_sex)
mode(df1$percent_time)

(p <- ggplot(df1, aes(x=strain_sex, y=percent_time, fill = strain_sex)) + 
    geom_violin(width=1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    
    xlab("") +
    ylab("Prop. time in top occupied zone") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.position = "none")
)
ggsave(p, filename = "output/move_data_prop_time_top_zone.png", device = "png", bg = "white")
ggsave(p, filename = "output/move_data_prop_time_top_zone.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# STATS
df1[sapply(df1, is.character)] <- lapply(df1[sapply(df1, is.character)],as.factor)
m1 = lmer(asin(sqrt(percent_time)) ~ strain*sex + (1|trial), data = df1) # random effects
AIC(m1)
anova(m1)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
emmeans(m1, pairwise ~ strain*sex, adjust = 'none')

summary(m1)
#C57-F untransformed estimate
sin(0.94079)^2
#C57-M untransformed estimate
sin(0.26718+0.94079)^2
#OB-F untransformed estimate
sin(0.25090 +0.94079)^2
#OB-M untransformed estimate
sin(0.34030+0.94079)^2

## notes
## there is a theoretical reason for going with X as it is proportion data. logit transformations cannot be used if there are 0 or 1
## in the case where there are 0s or 1s in proportion data, you need to do an arcsin transformatoin of the percentage data. 
## if you do the asin transformation, you need to UNDO it when you want to report the estimate values in the text
## to UNDO this, you need to do the sin(estimate)^2 

## you can use log transformations if there are 0s or 1s, you just add a constant to the all the data, for count data its common to do: log(x + 1)
## If its a time-based variable its standard to use a constant the is half the min value above 0
## i.e. sin(0.94)^2

