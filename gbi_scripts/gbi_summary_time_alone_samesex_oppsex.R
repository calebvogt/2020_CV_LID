## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)

# load_data ---------------------------------------------------------------
gbi <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_summary.csv")
gbi$name <- gsub("C57-M-","",gbi$name)
gbi$name <- gsub("C57-F-","",gbi$name)
gbi$name <- gsub("NYOB-M-","",gbi$name)
gbi$name <- gsub("NYOB-F-","",gbi$name)

# gbi_summary_perc_total_time_alone ------------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]), ## only one male or one female
            social_s=sum(duration_s[mf_sum>1]), ## more than one male or one female
            mf_s=sum(duration_s[m_sum>0&f_sum>0]), ## mixed sex any number of males and females greater than 0 
            mm_s=sum(duration_s[m_sum>1&f_sum==0]), ## males with no females
            ff_s=sum(duration_s[m_sum==0&f_sum>1]))  %>% ## females with no males
  mutate(samesex_s=if(sex=="M") mm_s else ff_s, ## time with only members of the same sex
         oppsex_s=mf_s, .after=social_s) %>% ## time in groups with 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.

df1 <- df %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(alone_s/sum_s*100), sd = sd(alone_s/sum_s*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day,y=mean,color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, by = 25)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% time alone") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "")
)
# ggsave(p, filename = "output/gbi_summary_perc_time_alone.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_perc_time_alone.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m1 = lmer(asin(sqrt(alone_s/sum_s)) ~ strain*sex + (1|trial), data = df)
summary(m1)
emmeans(m1, pairwise ~ strain*sex)

m2 = lmer(asin(sqrt(alone_s/sum_s)) ~ strain+sex + (1|trial), data = df)
summary(m2)

m3 = lmer(asin(sqrt(alone_s/sum_s)) ~ strain*sex*day + (1|trial), data = df)
summary(m3)

m4 = lmer(asin(sqrt(alone_s/sum_s)) ~ strain*sex*log(day) + (1|trial), data = df)
summary(m4)

AIC(m1,m2,m3,m4)
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# gbi_summary_perc_total_time_alone_violin --------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% ## drop daily measures
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s)

## get numbers of mice that spent majority of their total time alone. 
df1 <- df %>% 
  group_by(strain_sex,name) %>%
  mutate(total_perc_alone = alone_s/sum_s*100, 
         majority_alone=ifelse(total_perc_alone>50,1,0)) %>% 
  group_by(strain_sex,majority_alone) %>% 
  tally()

## plot
(p <- df %>% 
    filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
           !(name == "Anubis"), #T003: Anubis visually confirmed dead by seizure on day 5.  
           !(name == "Rae"), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
           !(name == "Hare"), #T004: Hare only appears day 1. Not recovered, presumed dead. 
           !(name == "Isis"), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
           !(name == "Rose")) %>% #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
    ggplot(., aes(strain_sex, y=alone_s/sum_s*100, fill=strain_sex))+ 
    geom_violin(alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    # geom_hline(yintercept=asin(sqrt(0.5))) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    ylab("% total time alone") +
    xlab("") +
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
          legend.position = "")
)
# ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") 

## STATS
m1 = lmer(asqrtAlone ~ strain*sex + (1|trial), data = time.Social.filtered)
AIC(m1)
anova(m1)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
emmeans(m1, pairwise ~ strain*sex)



# gbi_summary_perc_total_time_samesex -------------------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]), ## only one male or one female
            social_s=sum(duration_s[mf_sum>1]), ## more than one male or one female
            mf_s=sum(duration_s[m_sum>0&f_sum>0]), ## mixed sex any number of males and females greater than 0 
            mm_s=sum(duration_s[m_sum>1&f_sum==0]), ## males with no females
            ff_s=sum(duration_s[m_sum==0&f_sum>1]))  %>% ## females with no males
  mutate(samesex_s=if(sex=="M") mm_s else ff_s, ## time with only members of the same sex
         oppsex_s=mf_s, .after=social_s) %>% ## time in groups with 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.

df1 <- df %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(samesex_s/sum_s*100), sd = sd(samesex_s/sum_s*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day,y=mean,color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, by = 25)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% total time in same-sex associations") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "")
)
# ggsave(p, filename = "output/gbi_summary_perc_total_time_samesex.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_perc_total_time_samesex.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m1 = lmer(asin(sqrt(samesex_s/sum_s)) ~ strain*sex + (1|trial), data = df)
summary(m1)

m2 = lmer(asin(sqrt(samesex_s/sum_s)) ~ strain+sex + (1|trial), data = df)
summary(m2)

m3 = lmer(asin(sqrt(samesex_s/sum_s)) ~ strain*sex*day + (1|trial), data = df)
summary(m3)

m4 = lmer(asin(sqrt(samesex_s/sum_s)) ~ strain*sex*log(day) + (1|trial), data = df)
summary(m4)

AIC(m1,m2,m3,m4)
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# gbi_summary_perc_total_time_samesex_violin --------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% ## drop daily measures
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s)

## plot
(p <- df %>% 
    filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
           !(name == "Anubis"), #T003: Anubis visually confirmed dead by seizure on day 5.  
           !(name == "Rae"), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
           !(name == "Hare"), #T004: Hare only appears day 1. Not recovered, presumed dead. 
           !(name == "Isis"), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
           !(name == "Rose")) %>% #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
    ggplot(., aes(strain_sex, y=samesex_s/sum_s*100, fill=strain_sex))+ 
    geom_violin(alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    # geom_hline(yintercept=asin(sqrt(0.5))) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    ylab("% total time in same-sex associations") +
    xlab("") +
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
          legend.position = "")
)

### mikes code, deprecated
### Figure S3b. Proportion of time spent with asqrt Same Sex
(p <- ggplot(data=time.Social.filtered, aes(x=strain_sex, y=asqrtSS, fill=strain_sex)) +
    geom_violin(alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    # geom_hline(yintercept=asin(sqrt(0.5))) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    ylab("Prop. time in same sex groups (arcsine)") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = "")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") 


## STATS
m1 = lmer(asqrtSS ~ strain*sex + (1|trial), data = time.Social.filtered)
AIC(m1)
anova(m1)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
emmeans(m1, pairwise ~ strain*sex)



# gbi_summary_perc_social_time_samesex -------------------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]), ## only one male or one female
            social_s=sum(duration_s[mf_sum>1]), ## more than one male or one female
            mf_s=sum(duration_s[m_sum>0&f_sum>0]), ## mixed sex any number of males and females greater than 0 
            mm_s=sum(duration_s[m_sum>1&f_sum==0]), ## males with no females
            ff_s=sum(duration_s[m_sum==0&f_sum>1]))  %>% ## females with no males
  mutate(samesex_s=if(sex=="M") mm_s else ff_s, ## time with only members of the same sex
         oppsex_s=mf_s, .after=social_s) %>% ## time in groups with 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.

df1 <- df %>% 
  filter(!(social_s==0)) %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(samesex_s/social_s*100), sd = sd(samesex_s/social_s*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day,y=mean,color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, by = 25)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% social time in same-sex associations") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "")
)
# ggsave(p, filename = "output/gbi_summary_perc_social_time_samesex.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_perc_social_time_samesex.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# gbi_summary_perc_total_time_oppsex --------------------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]), ## only one male or one female
            social_s=sum(duration_s[mf_sum>1]), ## more than one male or one female
            mf_s=sum(duration_s[m_sum>0&f_sum>0]), ## mixed sex any number of males and females greater than 0 
            mm_s=sum(duration_s[m_sum>1&f_sum==0]), ## males with no females
            ff_s=sum(duration_s[m_sum==0&f_sum>1]))  %>% ## females with no males
  mutate(samesex_s=if(sex=="M") mm_s else ff_s, ## time with only members of the same sex
         oppsex_s=mf_s, .after=social_s) %>% ## time in groups with 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.


df1 <- df %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(oppsex_s/sum_s*100), sd = sd(oppsex_s/sum_s*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day,y=mean,color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, by = 25)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% time in opposite-sex associations") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "")
)
# ggsave(p, filename = "output/gbi_summary_perc_total_time_oppsex.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_summary_perc_total_time_oppsex.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m1 = lmer(asin(sqrt(oppsex_s/sum_s)) ~ strain*sex + (1|trial), data = df)
summary(m1)

m2 = lmer(asin(sqrt(oppsex_s/sum_s)) ~ strain+sex + (1|trial), data = df)
summary(m2)

m3 = lmer(asin(sqrt(oppsex_s/sum_s)) ~ strain*sex*day + (1|trial), data = df)
summary(m3)

m4 = lmer(asin(sqrt(oppsex_s/sum_s)) ~ strain*sex*log(day) + (1|trial), data = df)
summary(m4)

m5 = lmer(asin(sqrt(oppsex_s/sum_s)) ~ strain*sex*log(day) + (1|trial) + (1|name), data = df)
summary(m5)

AIC(m1,m2,m3,m4,m5)
write.table(summary(m5)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# gbi_summary_perc_total_time_oppsex_violin --------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% ## drop daily measures
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis"), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae"), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare"), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis"), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose")) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

## plot
(p <- ggplot(df, aes(strain_sex, y=oppsex_s/sum_s*100, fill=strain_sex))+ 
    geom_violin(alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    # geom_hline(yintercept=asin(sqrt(0.5))) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    ylab("% total time in opposite-sex associations") +
    xlab("") +
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
          legend.position = "")
)

### mikes code, deprecated
### Fig. S3a: Proportion of time spent with mixed sex
(p <- ggplot(data=time.Social.filtered, aes(x=strain_sex, y=asqrtFM, fill=strain_sex)) +
   geom_violin(alpha = 1) +
   geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
   # geom_hline(yintercept=asin(sqrt(0.5))) +
   scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                     values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
   ylab("Prop. time in mixed sex groups (arcsine)") +
   xlab("") +
   theme_classic() +
   theme(axis.text.x = element_text(color = "black", size = 8),
         axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
         axis.text.y = element_text(color = "black", size = 8),
         axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
         plot.background = element_rect(fill = "transparent", color = NA),
         panel.background = element_rect(fill = "transparent"), 
         legend.title = element_blank(),
         legend.text = element_text(size=6),
         legend.background = element_rect(fill='transparent'),
         legend.position = "")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") 

## STATS
m1 = lmer(asqrtFM ~ strain*sex + (1|trial), data = time.Social.filtered)
AIC(m1)
anova(m1)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
emmeans(m1, pairwise ~ strain*sex)


# gbi_summary_perc_social_time_oppsex -------------------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]), ## only one male or one female
            social_s=sum(duration_s[mf_sum>1]), ## more than one male or one female
            mf_s=sum(duration_s[m_sum>0&f_sum>0]), ## mixed sex any number of males and females greater than 0 
            mm_s=sum(duration_s[m_sum>1&f_sum==0]), ## males with no females
            ff_s=sum(duration_s[m_sum==0&f_sum>1]))  %>% ## females with no males
  mutate(samesex_s=if(sex=="M") mm_s else ff_s, ## time with only members of the same sex
         oppsex_s=mf_s, .after=social_s) %>% ## time in groups with 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.

df1 <- df %>% 
  filter(!(social_s==0)) %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(oppsex_s/social_s*100), sd = sd(oppsex_s/social_s*100),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day,y=mean,color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(0,100), breaks = seq(0, 100, by = 25)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% social time in opposite-sex associations") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "")
)
# ggsave(p, filename = "output/gbi_summary_perc_social_time_oppsex.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_summary_perc_social_time_oppsex.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")
## Stats
m1 = lmer(asin(sqrt(oppsex_s/social_s)) ~ strain*sex + (1|trial), data = df)
summary(m1)

m2 = lmer(asin(sqrt(oppsex_s/social_s)) ~ strain+sex + (1|trial), data = df)
summary(m2)

m3 = lmer(asin(sqrt(oppsex_s/social_s)) ~ strain*sex*day + (1|trial), data = df)
summary(m3)

m4 = lmer(asin(sqrt(oppsex_s/social_s)) ~ strain*sex*log(day) + (1|trial), data = df)
summary(m4)

m5 = lmer(asin(sqrt(oppsex_s/social_s)) ~ strain*sex*log(day) + (1|trial) + (1|name), data = df)
summary(m5)

AIC(m1,m2,m3,m4,m5)
write.table(summary(m5)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)






# gbi_summary_perc_social_time_oppsex_violin --------------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% ## drop daily measures
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis"), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae"), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare"), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis"), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose")) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

## get numbers of mice that spent majority of their social time in oppsex associations
df1 <- df %>% 
  group_by(strain_sex,name) %>%
  mutate(social_perc_oppsex = oppsex_s/social_s*100, 
         majority_social_oppsex=ifelse(social_perc_oppsex>50,1,0)) %>% 
  group_by(majority_social_oppsex) %>% 
  tally()

## plot
(p <- ggplot(df, aes(strain_sex, y=oppsex_s/social_s*100, fill=strain_sex))+ 
    geom_violin(alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    # geom_hline(yintercept=asin(sqrt(0.5))) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    ylab("% social time in opposite-sex associations") +
    xlab("") +
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
          legend.position = "")
)


### mikes code, deprecated
### Fig. S3a: Proportion of time spent with mixed sex
(p <- ggplot(data=time.Social.filtered, aes(x=strain_sex, y=asqrtFM, fill=strain_sex)) +
    geom_violin(alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    # geom_hline(yintercept=asin(sqrt(0.5))) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    ylab("Prop. time in mixed sex groups (arcsine)") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = "")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") 

## STATS
m1 = lmer(asqrtFM ~ strain*sex + (1|trial), data = time.Social.filtered)
AIC(m1)
anova(m1)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
emmeans(m1, pairwise ~ strain*sex)

