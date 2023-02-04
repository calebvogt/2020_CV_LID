## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggpubr)

# load data ---------------------------------------------------------------
## import gbi_summary
gbi <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_summary.csv")
gbi$name <- gsub("C57-M-","",gbi$name)
gbi$name <- gsub("C57-F-","",gbi$name)
gbi$name <- gsub("NYOB-M-","",gbi$name)
gbi$name <- gsub("NYOB-F-","",gbi$name)

## import pas values
pas0 <- read.csv("data/priority_access_scores_males.csv")
pas <- pas0 %>% 
  group_by(name) %>% 
  filter(noon_day == max(noon_day)) %>% 
  ungroup() %>% 
  mutate(pas = csum_daily_capture_penalty) %>% 
  group_by(trial) %>% 
  mutate(rank_order = rank(-pas)) ## create rank order
  

# gbi_summary_PAS_status_hilo_oppsex_csum_total_time_lineplot ----------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code,day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s)

df1 <- df %>% 
  filter(sex=="M") %>% 
  merge(., pas[,c("name", "pas", "rank_order")], by = "name", all.x =T) %>% 
  mutate(status = ifelse(pas>0,"High", "Low"),
         strain_status=paste0(strain,"-",status)) %>% 
  group_by(trial,strain,name) %>% 
  mutate(csum_oppsex_min = cumsum(oppsex_s/60)) %>%  
  ungroup() %>% 
  complete(name, day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, day) %>%
  group_by(name) %>% 
  fill(trial,strain, sex, oppsex_s,csum_oppsex_min) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.

df2 <- df1 %>% 
  group_by(strain, status,strain_status, day) %>% 
  summarise(mean = mean(csum_oppsex_min), sd = sd(csum_oppsex_min), count = n(), se = (sd/(sqrt(count))))

# PLOT
(p <- ggplot(df2, aes(x=day,y=mean,group = strain_status)) + 
    geom_line(aes(linetype = status, color = strain),size = 0.5) + 
    geom_point(aes(color = strain),size = 0.5) +
    geom_errorbar(aes(color = strain, ymin = mean - se, ymax = mean + se), width = 0.2) +
    xlab("Day") +
    ylab("Cumulative total time with females (min)") +
    scale_x_continuous(breaks = seq(1,11,by=1), limits = c(1,10.2)) +
    scale_color_manual(breaks = c("C57", "NYOB"),
                       values=c("sienna", "skyblue4")) +
    # scale_linetype_manual(values=c("solid","twodash"))+
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          # legend.position = c(0.15,0.8))
          legend.position = "")
)
# ggsave(p, filename = "output/gbi_summary_PAS_status_hilo_oppsex_csum_total_time_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_PAS_status_hilo_oppsex_csum_total_time_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m1 <- lmer(csum_oppsex_min ~ strain*status*log(day) + (1|trial), data=df1) 
m2 <- lmer(csum_oppsex_min ~ strain_status*log(day) + (1|trial), data=df1)
m3 <- lmer(csum_oppsex_min ~ strain_status*log(day) + (1|trial), data=subset(df1, strain == "C57")) 
AIC(m1,m2)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## range of total time low status OB males overlap with females 
df3 <- df1 %>% 
  filter(strain=="NYOB" & status=="Low") %>% 
  group_by(name) %>% 
  filter(csum_oppsex_min == max(csum_oppsex_min)) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis"), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" ), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare"), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis"), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose")) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.
  


# gbi_summary_PAS_rank_oppsex_prop_total_time_lineplot --------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s)

df1 <- df %>% 
  filter(sex=="M") %>% 
  merge(., pas[,c("name", "pas", "rank_order")], by = "name", all.x =T) 

##Plot
(p <- ggplot(df1, aes(x=rank_order, y=oppsex_s/sum_s*100, color=strain)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_color_manual(breaks = c("C57", "NYOB"),
                       values=c("sienna", "skyblue4")) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    xlab("Ranked male Priority Access Score") +
    ylab("% total time with females") +
    stat_cor(aes(color = strain, label=paste(..rr.label.., ..p.label..,sep = "~`,`~")), 
             label.x = 3, method = "pearson", p.accuracy = 0.0001) +
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
          # legend.position = c(0.9,0.9))
          legend.position="")
)
# ggsave(p, filename = "output/gbi_summary_PAS_rank_oppsex_perc_total_time_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_PAS_rank_oppsex_perc_total_time_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## stats
m1 <- lmer(asin(sqrt(oppsex_s/sum_s)) ~ rank_order*strain + (1|trial), data=df1) ## arcsinsqrt transform the proportion data
m2 <- lm(asin(sqrt(oppsex_s/sum_s)) ~ rank_order*strain, data=df1)
m3 <- lmer(oppsex_s/sum_s ~ rank_order*strain + (1|trial), data=df1) 
m4 <- lmer(asin(sqrt(oppsex_s/sum_s)) ~ rank_order + (1|trial), data=subset(df1, strain=="NYOB")) ## arcsinsqrt transform the proportion data
m4 <- lm(oppsex_s/sum_s*100 ~ rank_order, data=subset(df1, strain=="NYOB")) ## arcsinsqrt transform the proportion data
summary(m4)
m5 <- lmer(asin(sqrt(oppsex_s/sum_s)) ~ rank_order + (1|trial), data=subset(df1, strain=="C57")) ## arcsinsqrt transform the proportion data
summary(m5)
AIC(m1,m2,m3)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))

write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# gbi_summary_PAS_rank_oppsex_total_time_lineplot --------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s)

df1 <- df %>% 
  filter(sex=="M") %>% 
  merge(., pas[,c("name", "pas", "rank_order")], by = "name", all.x =T) 


##Plot
(p <- ggplot(df1, aes(x=rank_order, y=oppsex_s/3600, color=strain)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_color_manual(breaks = c("C57", "NYOB"),
                       values=c("sienna", "skyblue4")) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    xlab("Ranked male Priority Access Score") +
    ylab("Prop. social time with females (arcsine)") +
    stat_cor(aes(color = strain), label.x = 4, method = "pearson", p.accuracy = 0.001) +
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
          # legend.position = c(0.9,0.9))
          legend.position="")
)
ggsave(p, filename = "output/gbi_summary_MF_time_PAS_rank_lineplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_summary_MF_time_PAS_rank_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## stats
m1 <- lmer(asqrt_prop_f_social_time ~ rank_order*strain + (1|trial), data=df)
m2 <- lm(f_time_min ~ rank_order*strain, data=df)
AIC(m1,m2)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


## boxplot
(p <- ggplot(df, aes(x=as.factor(strain), y=asqrt_prop_f_social_time, fill = as.factor(status))) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    ggtitle("Day 1-10") + # adjust days if desired
    stat_compare_means(aes(group = status))
)
# ggsave(p, filename = "output/gbi_summary_M_PAS_status_F_time_boxplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_M_PAS_status_F_time_boxplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")



# gbi_summary_PAS_status_himedlo_oppsex_prop_total_time_lineplot --------------------------------------
df<-gbi %>% 
  mutate(strain_sex=paste0(strain,"-",sex), .after=trial) %>% 
  group_by(trial, strain_sex,strain, sex, name, code) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=if(sex=="M") mm_s else ff_s,
         oppsex_s=mf_s, .after=social_s)

df1 <- df %>% 
  filter(sex=="M") %>% 
  merge(., pas[,c("name", "pas", "rank_order")], by = "name", all.x =T) %>% 
  mutate(status = ifelse(pas>5,1, ifelse(pas < -5,3, 2))) # 1 = high, 2 = medium, 3 = low
  
##Plot
(p <- ggplot(df1, aes(x=status, y=oppsex_s/sum_s*100, color=strain)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_color_manual(breaks = c("C57", "NYOB"),
                       values=c("sienna", "skyblue4")) +
    # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    xlab("Priority Access Score status") +
    ylab("% total time with females") +
    stat_cor(aes(color = strain, label=paste(..rr.label.., ..p.label..,sep = "~`,`~")), 
             label.x = 3, method = "pearson", p.accuracy = 0.0001) +
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
          # legend.position = c(0.9,0.9))
          legend.position="")
)
ggsave(p, filename = "output/gbi_summary_PAS_rank_oppsex_perc_total_time_lineplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_summary_PAS_rank_oppsex_perc_total_time_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## stats
m1 <- lmer(asin(sqrt(oppsex_s/sum_s)) ~ rank_order*strain + (1|trial), data=df1) ## arcsinsqrt transform the proportion data
m2 <- lm(asin(sqrt(oppsex_s/sum_s)) ~ rank_order*strain, data=df1)
m3 <- lmer(oppsex_s/sum_s ~ rank_order*strain + (1|trial), data=df1) 
m4 <- lmer(asin(sqrt(oppsex_s/sum_s)) ~ rank_order + (1|trial), data=subset(df1, strain=="NYOB")) ## arcsinsqrt transform the proportion data
m4 <- lm(oppsex_s/sum_s*100 ~ rank_order, data=subset(df1, strain=="NYOB")) ## arcsinsqrt transform the proportion data
summary(m4)
m5 <- lmer(asin(sqrt(oppsex_s/sum_s)) ~ rank_order + (1|trial), data=subset(df1, strain=="C57")) ## arcsinsqrt transform the proportion data
summary(m5)
AIC(m1,m2,m3)
summary(m1)
qqnorm(resid(m1))
qqline(resid(m1))

write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# gbi_summary_male_PAS_rank_vs_prop_female_time ---------------------------
df <- gbi %>% 
  # filter(day %in% (6:10)) %>% ## adjust days
  group_by(name) %>% 
  mutate(total_all_time_min = sum(duration_s/60)) %>% 
  ungroup() %>% 
  filter(sex == "M", 
         mf_sum >1) %>% 
  group_by(name) %>% 
  mutate(total_social_time_min = sum(duration_s/60)) %>% 
  filter(sex == "M",
         f_sum > 0) %>% 
  group_by(trial, strain, name) %>% 
  mutate(f_time_min = sum(duration_s/60)) %>% 
  merge(., pas[,c("name", "pas", "rank_order")], by = "name", all.x =T) %>% 
  distinct(name, .keep_all = T) %>% 
  filter(!(name == "Anubis")) %>% 
  mutate(prop_f_all_time = f_time_min / total_all_time_min, 
         asqrt_prop_f_all_time = asin(sqrt(prop_f_all_time)), 
         prop_f_social_time = f_time_min / total_social_time_min, 
         asqrt_prop_f_social_time = asin(sqrt(prop_f_social_time)), 
         status = ifelse(pas > 0, "High", "Low"))
  
##Plot
(p <- ggplot(df, aes(x=rank_order, y=asqrt_prop_f_social_time, color=strain)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_color_manual(breaks = c("C57", "NYOB"),
                       values=c("sienna", "skyblue4")) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    xlab("Ranked male Priority Access Score") +
    ylab("Prop. social time with females (arcsine)") +
    stat_cor(aes(color = strain), label.x = 4, method = "pearson", p.accuracy = 0.001) +
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
          # legend.position = c(0.9,0.9))
          legend.position="")
)
ggsave(p, filename = "output/gbi_summary_MF_time_PAS_rank_lineplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_summary_MF_time_PAS_rank_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## stats
m1 <- lmer(asqrt_prop_f_social_time ~ rank_order*strain + (1|trial), data=df)
m2 <- lm(f_time_min ~ rank_order*strain, data=df)
AIC(m1,m2)
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


## boxplot
(p <- ggplot(df, aes(x=as.factor(strain), y=asqrt_prop_f_social_time, fill = as.factor(status))) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    ggtitle("Day 1-10") + # adjust days if desired
    stat_compare_means(aes(group = status))
)
# ggsave(p, filename = "output/gbi_summary_M_PAS_status_F_time_boxplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/gbi_summary_M_PAS_status_F_time_boxplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")





# Mike's code, outdated. male Priority Access Score Ranking and time spent in mixed sex groups --------
# setwd("C:/Users/caleb/Box/7_LID_2020/1_MS_files/Fig3.jan18.2022")
# setwd("C:/Users/caleb vogt/Box/7_LID_2020/1_MS_files/Fig3.jan18.2022")

## excel file mike gave me doesnt have any of the formulas saved. ffs. 

terrSocial <- read.csv("data/terr.TimeSocial.csv") 
terrSocial <- terrSocial %>%
  mutate(strain_sex = paste0(strain, "-", sex))

# output to clipboard
# write.table(terrSocial, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
# terrSocial <- read_excel("Figure_Data.xlsx", sheet = "Fig3d")

##Plot
(p <- ggplot(terrSocial, aes(x=rankTerr, y=asqrtMixedSex, color=strain_sex)) +
    # geom_point(aes(color=strain_sex)) +
    geom_point() +
    geom_smooth(method="lm") +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    xlab("Male Priority Access Score, Ranked") +
    ylab("Prop. time with females (arcsine)") +
    # stat_cor(aes(color = strain_sex), label.x = 4, method = "pearson", p.accuracy = 0.001) +
    # stat_cor(label.x = 4, method = "pearson", p.accuracy = 0.001) +
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
  # legend.position = "left")
  # legend.position = "")
)
# ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") 

terrSocial.mod1 = lmer(asqrtMixedSex ~ rankTerr*strain + (1|trial), data=terrSocial)
AIC(terrSocial.mod1)
anova(terrSocial.mod1) ## mike reported these outputs
summary(terrSocial.mod1)
write.table(summary(terrSocial.mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


