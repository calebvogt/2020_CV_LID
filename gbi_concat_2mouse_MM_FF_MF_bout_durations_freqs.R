## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## here, dyadic means that only two mice were estimated to be in the zone. 

library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom)
library(clipr)

data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT_GBI_concat.csv", stringsAsFactors = FALSE))

# gbi_concat_MM_dyadic_time_half_elapsed  ----------------------------------------------------------------------
df2 <- data %>% 
  group_by(trial) %>% 
  mutate(trial_time_min = as.numeric(difftime(field_time_start,min(field_time_start),units="min")+1)) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  mutate(trial_time_adj = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(m_sum == 2) %>% ## only looks at times when there are two males, indepe
  mutate(grouptype = "mm") %>% 
  mutate(duration_s = replace(duration_s, duration_s == 0, 1)) ## replace 0s (due to loss of ms in field time values) with 1s. 

df2$strain <- ifelse(df2$trial == "T001", "C57", 
                     ifelse(df2$trial == "T002", "C57", 
                            ifelse(df2$trial == "T003", "C57",
                                   ifelse(df2$trial == "T004", "Outbred",
                                          ifelse(df2$trial == "T005", "Outbred",
                                                 ifelse(df2$trial == "T006", "C57",
                                                        ifelse(df2$trial == "T007", "Outbred", NA)))))))
## get trial 50% time estimate
df3 <- df2 %>% 
  group_by(trial) %>% 
  arrange(field_time_start) %>% 
  mutate(trial_duration_min = sum((duration_s)/60),
         csum_min = cumsum((duration_s)/60),
         perc_dyadic_time_elapsed = (csum_min/trial_duration_min)*100) %>% 
  slice(which.min(abs(perc_dyadic_time_elapsed - 55)))

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

# gbi_concat_MM_dyadic_frequencies  ----------------------------------------------------------------------
df2 <- data %>% 
  group_by(trial) %>% 
  mutate(trial_time_min = as.numeric(difftime(field_time_start,min(field_time_start),units="min")+1)) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  mutate(trial_time_adj = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(m_sum == 2) %>% 
  mutate(grouptype = "mm") %>% 
  mutate(duration_s = replace(duration_s, duration_s == 0, 1))  ## replace 0s (due to loss of ms in field time values) with 1s. 
  
df2$strain <- ifelse(df2$trial == "T001", "C57", 
                     ifelse(df2$trial == "T002", "C57", 
                            ifelse(df2$trial == "T003", "C57",
                                   ifelse(df2$trial == "T004", "Outbred",
                                          ifelse(df2$trial == "T005", "Outbred",
                                                 ifelse(df2$trial == "T006", "C57",
                                                        ifelse(df2$trial == "T007", "Outbred", NA)))))))

df3 <- df2 %>% 
  group_by(trial, strain, day) %>% 
  tally() %>% 
  mutate(mm_bouts = n)

## boxplot
ggplot(df3, aes(x=as.factor(day), y=mm_bouts, fill = as.factor(strain))) +
  geom_boxplot()

## lineplot
df4 <- df3 %>% 
  group_by(strain, day) %>%
  summarise(mean = mean(mm_bouts), sd = sd(mm_bouts),count = n(), se = (sd/(sqrt(count))), mode = mode(mm_bouts))

(p <- ggplot(df4, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    xlab("Day") +
    ylab("Male-male dyadic social bout trial frequency") + 
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
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
# ggsave(p, filename = "output/gbi_concat_mm_dyadic_frequency_line.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_concat_mm_dyadic_frequency_line.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m0 <- lmer(log(mm_bouts) ~ strain*log(day) + (1|trial), data=df3) 
summary(m0)
qqnorm(resid(m0))
qqline(resid(m0))
write.table(summary(m0)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# gbi_concat_MM_dyadic_durations ----------------------------------------------------------------------
df2 <- data %>% 
  group_by(trial) %>% 
  mutate(trial_time_min = as.numeric(difftime(field_time_start,min(field_time_start),units="min")+1)) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  mutate(trial_time_adj = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(m_sum == 2) %>% 
  mutate(grouptype = "mm") %>% 
  mutate(duration_s = replace(duration_s, duration_s == 0, 1)) ## replace 0s (due to loss of ms in field time values) with 1s. 

df2$strain <- ifelse(df2$trial == "T001", "C57", 
                    ifelse(df2$trial == "T002", "C57", 
                           ifelse(df2$trial == "T003", "C57",
                                  ifelse(df2$trial == "T004", "Outbred",
                                         ifelse(df2$trial == "T005", "Outbred",
                                                ifelse(df2$trial == "T006", "C57",
                                                       ifelse(df2$trial == "T007", "Outbred", NA)))))))


## boxplot
(p <- ggplot(df2, aes(x=as.factor(day), y=(duration_s/60), fill = as.factor(strain))) +
    geom_boxplot() +
    ylim(0,40)
)


## lineplot
df3 <- df2 %>% 
  group_by(strain, day) %>%
  summarise(mean = mean((duration_s/60)), sd = sd((duration_s/60)),count = n(), se = (sd/(sqrt(count))), mode = mode((duration_s/60)))

(p <- ggplot(df3, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    xlab("Day") +
    ylab("Male-male dyadic social bout durations (min)") + 
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "right" )
)
# ggsave(p, filename = "output/gbi_concat_mm_dyadic_durations_line.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_concat_mm_dyadic_durations_line.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


## scatterplot
(p <- ggplot(data=df2, aes(x=trial_time_adj, y=(duration_s/60), color=strain))+
    geom_point(size = 0.9) +
    xlab("Day") +
    labs(y=paste0("Male-male dyadic spatiotemporal","\n", "association bout durations (min)")) +
    facet_wrap(~strain) +
    scale_color_manual(breaks = c("C57", "Outbred"), 
                       values=c("sienna", "skyblue4"),
                       labels=c("C57 males", "Outbred males")) +
    scale_x_continuous(limits = c(1,11), breaks = seq(1, 10, by = 1)) +
    scale_y_continuous(limits = c(0,150), breaks = seq(0,150, by = 25)) +
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
          legend.position = c(0.5,0.9)) +
    guides(color=guide_legend(override.aes = list(size=1.5)))
)
# ggsave(p, filename = "output/gbi_concat_MM_dyadic_durations_scatterplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_concat_MM_dyadic_durations_scatterplot.svg", device = "svg", width=5, height=2.15, bg = "transparent")

## Stats
hist(df2$duration_s)
hist(log(df2$duration_s)) ## normalize with log

m0 <- lmer(log(duration_s/60) ~ strain*log(trial_time_adj) + (1|trial), data=df2)
m1 <- lmer(log(duration_s/60) ~ strain+log(day) + (1|trial),  data=df2) 
m2 <- lmer(log(duration_s/60) ~ strain*log(day) + (1|trial),  data=df2) 
qqnorm(resid(m2))
qqline(resid(m2))
AIC(m2)
summary(m2)
write.table(summary(m2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# daily models: change for days 1 through 10 and report
t.test((duration_s/60) ~ strain, subset(df2,day==2)) 
m1 <- lm((duration_s/60) ~ strain, subset(df2,day==2)) %>% summary()
coef(m1)


# gbi_concat_FM_dyadic_durations.png  ----------------------------------------------------------------------
df2 <- data %>% 
  group_by(trial) %>% 
  mutate(trial_time = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(m_sum == 1 & f_sum == 1) %>% 
  mutate(grouptype = "mf") %>%
  filter(!(duration_s == 0)) # doesnt change anything

df2$strain <- ifelse(df2$trial == "T001", "C57", 
                     ifelse(df2$trial == "T002", "C57", 
                            ifelse(df2$trial == "T003", "C57",
                                   ifelse(df2$trial == "T004", "Outbred",
                                          ifelse(df2$trial == "T005", "Outbred",
                                                 ifelse(df2$trial == "T006", "C57",
                                                        ifelse(df2$trial == "T007", "Outbred", NA)))))))


## scatterplot
(p <- ggplot(data=df2, aes(y=(duration_s/60), x=trial_time, color=strain))+
    geom_point(size = 0.9) +
    xlab("Day") +
    labs(y=paste0("Female-male dyadic spatiotemporal","\n", "association bout durations (min)")) +
    facet_wrap(~strain) +
    scale_color_manual(breaks = c("C57", "Outbred"), 
                       values=c("goldenrod1", "goldenrod4")) +
    scale_x_continuous(limits = c(1,11), breaks = seq(1, 10, by = 1)) +
    ylim(0,40) +
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
    guides(color=guide_legend(override.aes = list(size=1.5)))
)
# ggsave(p, filename = "output/gbi_concat_FM_dyadic_durations_scatterplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_concat_FM_dyadic_durations_scatterplot.svg", device = "svg", width=5, height=2.15, bg = "transparent")

## Stats
hist(df2$duration_s)
hist(log(df2$duration_s)) ## normalize with log

m1 <- lmer((duration_s/60)~strain*trial_time+(1|trial), data=df2)
m2 <- lmer(log(duration_s/60) ~ strain*log(day) + (1|trial),  data=df2) 
m3 <- lmer(log(duration_s/60) ~ strain*trial_time + (1|trial),  data=df2) 
m4 <- lmer(log(duration_s/60) ~ strain*log(trial_time) + (1|trial),  data=df2) 
AIC(m1, m2,m3,m4)
summary(m3)
write.table(summary(m3)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# gbi_concat_FF_dyadic_durations.png  ----------------------------------------------------------------------
df2 <- data %>% 
  group_by(trial) %>% 
  mutate(trial_time = (((as.numeric(difftime(field_time_start,min(field_time_start),units="secs")+1))/864000)*10)+1) %>% ## time from first observed movebout GBI flocking event, including single animal movements
  ungroup() %>% 
  filter(f_sum == 2) %>%  ## original, duration of time with only two females present in the box.
  mutate(grouptype = "ff") %>%
  filter(!(duration_s == 0)) # doesnt change anything

df2$strain <- ifelse(df2$trial == "T001", "C57", 
                     ifelse(df2$trial == "T002", "C57", 
                            ifelse(df2$trial == "T003", "C57",
                                   ifelse(df2$trial == "T004", "Outbred",
                                          ifelse(df2$trial == "T005", "Outbred",
                                                 ifelse(df2$trial == "T006", "C57",
                                                        ifelse(df2$trial == "T007", "Outbred", NA)))))))

df3 <- df2 %>% 
  group_by(strain, day) %>%
  summarise(mean = mean((duration_s/60)), sd = sd((duration_s/60)),count = n(), se = (sd/(sqrt(count))), mode = mode((duration_s/60)))

## scatterplot
(p <- ggplot(data=df2, aes(y=(duration_s/60), x=trial_time, color=strain))+
    geom_point(size = 0.5) +
    xlab("Day") +
    labs(y=paste0("Female-female dyadic spatiotemporal","\n", "association bout durations (min)")) +
    facet_wrap(~strain) +
    scale_color_manual(breaks = c("C57", "Outbred"), 
                       values=c("sienna1", "skyblue"), 
                       labels=c("C57 females", "Outbred females")) +
    scale_x_continuous(limits = c(1,11), breaks = seq(1, 10, by = 1)) +
    ylim(0,40) +
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
# ggsave(p, filename = "output/gbi_concat_FF_dyadic_durations_scatterplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_concat_FF_dyadic_durations_scatterplot.svg", device = "svg", width=5, height=2.15, bg = "transparent")

## Stats
hist(df2$duration_s)
hist(log(df2$duration_s)) ## normalize with log

m1 <- lmer((duration_s/60)~strain*trial_time+(1|trial), data=df2)
m2 <- lmer(log(duration_s/60) ~ strain*log(day) + (1|trial),  data=df2) 
m3 <- lmer(log(duration_s/60) ~ strain*trial_time + (1|trial),  data=df2) 
m4 <- lmer(log(duration_s/60) ~ strain*log(trial_time) + (1|trial),  data=df2) 
AIC(m1, m2,m3,m4)
summary(m4)
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## lineplot
(p <- ggplot(df3, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
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


## boxplot
(p <- ggplot(df2, aes(x=as.factor(day), y=(duration_s/60), fill = as.factor(strain))) +
    geom_boxplot() +
    ylim(0,40)
)


# deprecated analysis ------------------------------------------------------------
## based off an old version of the MOVEBOUT_GBI_summary, with 153 grouping cutoff. 
## this would explain the different values....

# setwd("C:/Users/caleb/Box/7_LID_2020/1_MS_files/Fig3.jan18.2022")
# ff <- as.data.frame(fread("data/liddell2020.ff_copy.csv", stringsAsFactors = FALSE))
mm <- as.data.frame(fread("data/liddell2020.mm_copy.csv", stringsAsFactors = FALSE))
# ff=read_csv("liddell2020.ff.csv")
# fm=read_csv("liddell2020.fm.csv")
# mm=read_csv("liddell2020.mm.csv")
mm$TrialTime <- mm$TrialTime+1
# ff$TrialTime <- ff$TrialTime+1
# fm$TrialTime <- fm$TrialTime+1


## graph
(p <- ggplot(data=mm, aes(y=(duration_s/60), x=TrialTime, color=Strain))+
    geom_point(size = 0.9) +
    ylab("mikeMale-male social bout durations (min)") + 
    xlab("Day") +
    facet_wrap(~Strain) +
    scale_color_manual(breaks = c("C57", "NY"), values=c("sienna", "skyblue4")) +
    scale_x_continuous(limits = c(0.8,11), breaks = seq(0, 10, by = 1)) +
    # ylim(0,40) +
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
# ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=7, height=2.5, bg = "transparent") 

mm.mod1=lmer((duration_s/60)~Strain*TrialTime+(1|Trial), data=mm)
anova(mm.mod1)
summary(mm.mod1)
# write.table(summary(mm.mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)




# Fig. S3d: FF
# write.table(ff, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
# ff <- read_excel("Figure_Data.xlsx", sheet = "FigS3d")

(p <- ggplot(data=ff, aes(y=(duration_s/60), x=TrialTime, color=Strain))+
    geom_point(size = 0.9) +
    ylab("Female-female social bout durations (min)") + 
    xlab("Day") +
    facet_wrap(~Strain) +
    scale_color_manual(breaks = c("C57", "NY"), values=c("sienna1", "skyblue")) +
    scale_x_continuous(limits = c(0.8,11), breaks = seq(0, 10, by = 1)) +
    ylim(0,40) +
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
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=7, height=2.5, bg = "transparent") 

## Stats
ff.mod1=lmer((duration_s/60)~Strain*TrialTime+(1|Trial), data=ff)
anova(ff.mod1)
summary(ff.mod1)
write.table(summary(ff.mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)




# Fig. S3c: FM, mixed sex interactions
fm <- read_excel("Figure_Data.xlsx", sheet = "FigS3c")

(p <- ggplot(data=fm, aes(y=(duration_s/60), x=TrialTime, color=Strain))+
    geom_point(size = 0.9) +
    ylab("Female-male social bout durations (min)") + 
    xlab("Day") +
    facet_wrap(~Strain) +
    scale_color_manual(breaks = c("C57", "NY"), values=c("goldenrod1", "goldenrod4")) +
    scale_x_continuous(limits = c(0.8,11), breaks = seq(0, 10, by = 1)) +
    ylim(0,40) +
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
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=7, height=2.5, bg = "transparent") 

##stats
write.table(fm, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
fm.mod1=lmer((duration_s/60)~Strain*TrialTime+(1|Trial), data= fm)
anova(fm.mod1)
summary(fm.mod1)
write.table(summary(fm.mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

