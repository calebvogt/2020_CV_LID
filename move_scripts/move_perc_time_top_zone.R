## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, strain_sex, strain,sex,name, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(prop_time = n / sum(n)) %>% 
  arrange(desc(prop_time)) %>% 
  group_by(name) %>% 
  slice_max(prop_time, n = 2) %>% 
  mutate(rank_order = rank(desc(prop_time))) %>% 
  complete(rank_order = 1:2, fill = list(n = 0, prop_time = 0)) %>% 
  mutate(rank_order = as.factor(rank_order)) %>% 
  dplyr::select(trial, name, strain, sex, rank_order, zone, n, prop_time) %>% 
  filter(rank_order == 1) %>% 
  dplyr::rename(duration_min = n)

df1 <- df %>% 
  select(name, rank_order, zone) %>% 
  rename(top_zone = zone) ## top zone across all days

## get prop time in most occupied zone across all days. 
df2 <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, strain_sex, strain,sex,name, noon_day, zone) %>% 
  tally(sum(duration_min)) %>% 
  rename(duration_min = n) %>% 
  ungroup() %>% 
  group_by(trial, strain_sex, strain,sex,name, noon_day) %>%
  complete(zone = 1:8, fill = list(duration_min = 0)) %>%  # fill in days where mouse doesnt appear with 0s.
  mutate(prop_time = duration_min / sum(duration_min)) %>% #grouped by noon_day
  merge(., df1, by="name") %>% 
  filter(zone == top_zone)

# move_perc_time_top_zone_M.png --------------------------------------------------------------
df3 <- df2 %>% 
  filter(sex == "M") %>% 
  group_by(strain_sex, noon_day) %>%
  summarise(mean = mean(prop_time*100), sd = sd(prop_time*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df3, aes(x = noon_day, y = mean, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% time in most occupied zone") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "" )
)
# ggsave(p, filename = "output/move_perc_time_top_zone_M.png", device = "png", bg = "white")
# ggsave(p, filename = "output/move_perc_time_top_zone_M.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats on the proportion data
df2$trial <- as.factor(df2$trial)
df2$strain <- as.factor(df2$strain)
df2$noon_day <- as.numeric(df2$noon_day)

m1=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(1|name), data = subset(df2, sex=="M")) 
m2=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(1|trial), data = subset(df2, sex=="M")) 
m3=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(1|name)+(1|trial), data = subset(df2, sex=="M")) 
m4=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(log(noon_day)|name), data = subset(df2, sex=="M")) 
m5=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(log(noon_day)|name)+(1|trial), data = subset(df2, sex=="M")) 
AIC(m1,m2,m3,m4,m5)

summary(m4)
qqnorm(resid(m4))
qqline(resid(m4))
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# move_perc_time_top_zone_F.png --------------------------------------------------------------
df3 <- df2 %>% 
  filter(sex == "F") %>% 
  group_by(strain_sex, noon_day) %>%
  summarise(mean = mean(prop_time*100), sd = sd(prop_time*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df3, aes(x = noon_day, y = mean, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% time in most occupied zone") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "" )
)
# ggsave(p, filename = "output/move_perc_time_top_zone_F.png", device = "png", bg = "white")
# ggsave(p, filename = "output/move_perc_time_top_zone_F.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats on the proportion data
df2$trial <- as.factor(df2$trial)
df2$strain <- as.factor(df2$strain)
df2$noon_day <- as.numeric(df2$noon_day)

m1=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(1|name), data = subset(df2, sex=="F")) 
m2=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(1|trial), data = subset(df2, sex=="F")) 
m3=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(1|name)+(1|trial), data = subset(df2, sex=="F")) 
m4=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(log(noon_day)|name), data = subset(df2, sex=="F")) 
m5=lmer(asin(sqrt(prop_time))~strain*log(noon_day)+(log(noon_day)|name)+(1|trial), data = subset(df2, sex=="F")) 
summary(m5)
qqnorm(resid(m5))
qqline(resid(m5))

AIC(m1,m2,m3,m4,m5)
write.table(summary(m5)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# all strain and sex combinations -----------------------------------------
df3 <- df2 %>% 
  group_by(strain_sex, noon_day) %>%
  summarise(mean_n = mean(prop_time), sd_n = sd(prop_time),count = n(), se_n = (sd_n/(sqrt(count))))

(p <- ggplot(df3, aes(x = noon_day, y = mean_n, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    # values=c("red1", "red4", "steelblue", "steelblue4")) +
    # values=c("goldenrod1", "goldenrod4", "slateblue", "slateblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% time in most occupied zone") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "" )
)
# ggsave(p, filename = "output/move_data_perc_time_most_occupied_zone_daily.png", device = "png", bg = "white")
# ggsave(p, filename = "output/move_data_perc_time_most_occupied_zone_daily.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS
df2$trial <- as.factor(df2$trial)
df2$sex <- as.factor(df2$sex)
df2$strain <- as.factor(df2$strain)
df2$noon_day <- as.numeric(df2$noon_day)
# df2$c57F <- ifelse(df2$strain=="C57"&df2$sex=="F",1,0)

m1 = lmer(asin(sqrt(prop_time)) ~ strain*sex*log(noon_day) + (1|trial) + (1+log(noon_day)|name), data = df2) 
AIC(m1)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)

write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# # daily models: change for days 1 through 10. 
# d = lmer(asin(sqrt(prop_time)) ~ strain*sex + (1|trial), data = subset(df2, noon_day == 1)) 
# em.d <- emmeans(d, pairwise ~ strain*sex) # tukey adjusted
# write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


############### percentage across all days
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial, strain_sex, strain,sex,name, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(prop_time = n / sum(n)) %>% 
  arrange(desc(prop_time)) %>% 
  group_by(name) %>% 
  slice_max(prop_time, n = 2) %>% 
  mutate(rank_order = rank(desc(prop_time))) %>% 
  complete(rank_order = 1:2, fill = list(n = 0, prop_time = 0)) %>% 
  mutate(rank_order = as.factor(rank_order)) %>% 
  dplyr::select(trial, name, strain_sex, strain, sex, rank_order, zone, n, prop_time) %>% 
  filter(rank_order == 1) %>% 
  dplyr::rename(duration_min = n)

table(round(df$prop_time, 6), df$strain_sex)
mode(df$prop_time)

(p <- ggplot(df, aes(x=strain_sex, y=prop_time, fill = strain_sex)) + 
    geom_violin(width=1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    scale_y_continuous(limits = c(0,1)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    
    xlab("") +
    ylab("% time in most occupied zone") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.position = "none")
)
ggsave(p, filename = "output/move_data_perc_time_most_occupied_zone.png", device = "png", bg = "white")
ggsave(p, filename = "output/move_data_perc_time_most_occupied_zone.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# STATS
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
m1 = lmer(asin(sqrt(prop_time)) ~ strain*sex + (1|trial), data = df) # random effects
AIC(m1)
anova(m1)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

em.d <- emmeans(m1, pairwise ~ strain*sex, adjust = 'none')
write.table(em.d$c, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

summary(m1)
#C57-F untransformed estimate
sin(0.94079)^2
#C57-M untransformed estimate
sin(0.26718+0.94079)^2
#OB-F untransformed estimate
sin(0.25090 +0.94079)^2
#OB-M untransformed estimate
sin(0.34030+0.94079)^2


