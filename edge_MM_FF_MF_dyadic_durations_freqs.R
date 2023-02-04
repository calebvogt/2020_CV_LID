## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(magrittr)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

data <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_edgelist.csv")
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

df0 <- data %>% 
  merge(., metadata[,c("name","code", "sex", "strain")], by.x="ID1", by.y="name") %>% 
  rename(ID1_code=code, ID1_sex=sex, field_time_start=start, field_time_stop=stop) %>% 
  merge(., metadata[,c("name","code", "sex")], by.x="ID2", by.y="name") %>% 
  rename(ID2_code=code, ID2_sex=sex) %>% 
  select(trial, strain, day, zone, field_time_start, field_time_stop, duration_s, 
         ID1, ID1_code, ID1_sex,ID2, ID2_code, ID2_sex) %>% 
  mutate(grouptype=paste0(ID1_sex,ID2_sex)) %>% 
  mutate(strain=ifelse(strain=="NYOB", "Outbred","C57"))
  
# test <- df8 %>% 
#   group_by(trial, ID1, ID2,day) %>% 
#   summarize(sum_duration_s = sum(duration_s))

# edge_MM_time_half_elapsed  ----------------------------------------------------------------------
df2 <- df0 %>% 
  group_by(trial) %>% 
  mutate(trial_time_min = as.numeric(difftime(field_time_start,min(field_time_start),units="min")+1)) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  mutate(trial_time_adj = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(grouptype == "MM") %>% 
  mutate(duration_s = replace(duration_s, duration_s == 0, 1)) ## replace 0s (due to loss of ms in field time values) with 1s. 

## get trial 50% time estimate
df3 <- df2 %>% 
  group_by(trial) %>% 
  arrange(field_time_start) %>% 
  mutate(trial_duration_min = sum((duration_s)/60),
         csum_min = cumsum((duration_s)/60),
         perc_time_elapsed = (csum_min/trial_duration_min)*100) %>% 
  slice(which.min(abs(perc_time_elapsed - 50))) ## get 50 percent

df3 %>% 
  group_by(strain) %>% 
  summarize(mean = mean(trial_time_min))


##graphs
ggplot(df3, aes(x=strain, y= trial_time_min)) +
  geom_boxplot() +
  geom_point()
boxplot(trial_time_min ~ strain, data = subset(df3, df3$trial!="T003"))

## stats
hist(df3$trial_time_min)
hist(log(df3$trial_time_min)) ## better to log

t.test(log(trial_time_min) ~ as.factor(strain),df3)
m0 <- lm(log(trial_time_min)~as.factor(strain), df3) 
summary(m0)
m0 %>% broom::tidy() %>% 
  clipr::write_clip()

write.table(summary(m0)$coefficients, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## Get summary statistics
df4 <- df3 %>% 
  group_by(strain) %>% 
  summarise(mean = mean(trial_time_min), sd = sd(trial_time_min),count = n(), se = (sd/(sqrt(count))), mode = mode(trial_time_min))


# edge_MM_durations.png  ----------------------------------------------------------------------
df <- df0 %>% 
  group_by(trial) %>% 
  arrange(field_time_start) %>% 
  mutate(trial_time = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## note time is a prop out of 10. time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(grouptype == "MM") %>%
  filter(!(duration_s == 0)) # doesnt change anything

df1 <- df %>% 
  group_by(strain, day) %>%
  summarise(mean = mean((duration_s/60)), sd = sd((duration_s/60)),count = n(), se = (sd/(sqrt(count))), mode = mode((duration_s/60)))

## boxplot
(p <- ggplot(df, aes(x=as.factor(day), y=(duration_s/60), fill = as.factor(strain))) +
    geom_boxplot() +
    ylim(0,40)
)

## lineplot
(p <- ggplot(df1, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Male-Male social bout durations (min)") + 
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
ggsave(p, filename = "output/edge_MM_durations_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/edge_MM_durations_lineplot.png", device = "svg", width=5, height=2.15, bg = "transparent")

## scatterplot
(p <- ggplot(data=df, aes(y=(duration_s/60), x=trial_time, color=strain))+
    geom_point(size = 0.9) +
    facet_wrap(~strain) +
    scale_color_manual(breaks = c("C57", "Outbred"), values=c("sienna", "skyblue4")) +
    scale_x_continuous(limits = c(1,11), breaks = seq(1, 10, by = 1)) +
    ylim(0,40) +
    ylab("Male-Male social bout durations (min)") + 
    xlab("Day") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.9,0.9))
  # legend.position = "")
)
ggsave(p, filename = "output/edge_MM_durations_scatterplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/edge_MM_durations_scatterplot.svg", device = "svg", width=5, height=2.15, bg = "transparent")

## Stats
m1=lmer((duration_s/60)~strain*trial_time+(1|trial), data=df)
summary(m1)
AIC(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# edge_MF_durations.png  ----------------------------------------------------------------------
df <- df0 %>% 
  group_by(trial) %>% 
  arrange(field_time_start) %>% 
  mutate(trial_time = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## note time is a prop out of 10. time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(grouptype == "MF" | grouptype == "FM") %>%
  filter(!(duration_s == 0)) # doesnt change anything

df1 <- df %>% 
  group_by(strain, day) %>%
  summarise(mean = mean((duration_s/60)), sd = sd((duration_s/60)),count = n(), se = (sd/(sqrt(count))), mode = mode((duration_s/60)))


## boxplot
(p <- ggplot(df, aes(x=as.factor(day), y=(duration_s/60), fill = as.factor(strain))) +
    geom_boxplot() +
    ylim(0,40)
)

## lineplot
(p <- ggplot(df1, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("edgelist Male-female  social bout durations (min)") + 
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
# ggsave(p, filename = "output/edge_MF__durations_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/edge_MF__durations_lineplot.png", device = "svg", width=5, height=2.15, bg = "transparent")

## scatterplot
(p <- ggplot(data=df, aes(y=(duration_s/60), x=trial_time, color=strain))+
    geom_point(size = 0.9) +
    facet_wrap(~strain) +
    scale_color_manual(breaks = c("C57", "Outbred"), values=c("goldenrod1", "goldenrod4")) +
    scale_x_continuous(limits = c(1,11), breaks = seq(1, 10, by = 1)) +
    ylim(0,40) +
    ylab("Male-female  social bout durations (min)") +  
    xlab("Day") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.9,0.9))
  # legend.position = "")
)
ggsave(p, filename = "output/edge_MF__durations_scatterplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_concat_ff__durations.svg", device = "svg", width=5, height=2.15, bg = "transparent")

## Stats
m1=lmer((duration_s/60)~strain*trial_time+(1|trial), data=df)
summary(m1)
AIC(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# edge_FF_durations.png  ----------------------------------------------------------------------
df <- df0 %>% 
  group_by(trial) %>% 
  arrange(field_time_start) %>% 
  mutate(trial_time = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## note time is a prop out of 10. time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(grouptype == "FF") %>% 
  filter(!(duration_s == 0)) # doesnt change anything

df1 <- df %>% 
  group_by(strain, day) %>%
  summarise(mean = mean((duration_s/60)), sd = sd((duration_s/60)),count = n(), se = (sd/(sqrt(count))), mode = mode((duration_s/60)))

## scatterplot
(p <- ggplot(data=df, aes(y=(duration_s/60), x=trial_time, color=strain))+
    geom_point(size = 0.5) +
    ylab("Female-female spatiotemporal association bout durations (min)") + 
    xlab("Day") +
    facet_wrap(~strain) +
    scale_color_manual(breaks = c("C57", "Outbred"), 
                       values=c("sienna1", "skyblue"), 
                       labels=c("C57 females", "Outbred females")) +
    scale_x_continuous(limits = c(1,11), breaks = seq(1, 10, by = 1)) +
    ylim(0,40) +
    labs(y=paste0("Female-female spatiotemporal","\n", "association bout durations (min)")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.key.size=unit(1,"line"), ## set legend spacing
          legend.position = c(0.9,0.9)) +
          # legend.position = "")
    guides(color=guide_legend(override.aes = list(size=1.5)))
)
# ggsave(p, filename = "output/edge_FF_durations_scatterplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/edge_FF_durations_scatterplot.svg", device = "svg", width=5, height=2.15, bg = "transparent")

## Stats
m1=lmer((duration_s/60)~strain*trial_time+(1|trial), data=df)
summary(m1)
AIC(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

options(scipen = 999)


## lineplot
(p <- ggplot(df1, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Female-female social bout durations (min)") + 
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
ggsave(p, filename = "output/edge_FF_durations_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/edge_FF_durations_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## boxplot
(p <- ggplot(df, aes(x=as.factor(day), y=(duration_s/60), fill = as.factor(strain))) +
    geom_boxplot() +
    ylim(0,40)
)


