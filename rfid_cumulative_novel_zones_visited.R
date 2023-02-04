## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

df <- rfid_data %>% 
  select(trial, strain_sex,strain, sex, name,noon_day,zone) %>%
  group_by(trial, strain_sex,  name, zone) %>% 
  distinct(zone, .keep_all = T) %>% # get distinct zones ever visited, keep associated day it was visited
  group_by(trial, strain_sex,strain, sex, name, noon_day) %>% #regroup by day to count novel zones visited per day. 
  tally() %>% #tally the number of unique zones visited within a single day
  mutate(csum_novel_zones = cumsum(n)) %>%  ## get the csum of 
  complete(noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, noon_day) %>% 
  fill(csum_novel_zones) %>% ## fill cumulative sum data from last observed day
  select(trial, strain_sex,strain, sex,  name, noon_day,csum_novel_zones) %>% 
  na.omit(csum_novel_zones) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

# summary -----------------------------------------------------------------
df %>% 
  filter(noon_day==10) %>% 
  group_by(strain_sex,strain,sex) %>% 
  summarise(mean = mean(csum_novel_zones), sd = sd(csum_novel_zones),count = n(), se = (sd/(sqrt(count))))

# cumulative sum novel zones visited over time ----------------------------
df1 <- df %>% 
  group_by(strain_sex,strain,sex, noon_day) %>% 
  summarise(mean = mean(csum_novel_zones), sd = sd(csum_novel_zones),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=noon_day, y=mean, color = strain_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean+ se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1,11,by=1)) +
    scale_y_continuous(limits =c(1,8), breaks = seq(1,8,by=1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day") +
    ylab("Cumulative novel zones visited") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.position = "none") +
    guides(color=guide_legend(override.aes=list(fill=NA)))
)
# ggsave(p, filename = "output/rfid_cumulative_novel_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_cumulative_novel_zones_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
df$csum_novel_zones <- as.numeric(df$csum_novel_zones)

m2 = lmer(csum_novel_zones ~ strain*sex + (1|trial), data = subset(df, noon_day == 10))
summary(m2)
write.table(summary(m2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

em.m2 = emmeans(m2, pairwise ~ sex*strain)
write.table(em.m2$emmeans, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
write.table(em.m2$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


#logistic regression model of differences in all zones visited by strain*sex interaction
df2 <- df
df2$every_zone <- ifelse(df2$csum_novel_zones == 8, 1, 0)

m3 = glmer(every_zone ~ strain_sex + (1|trial), data = subset(df2, noon_day == 10), family ="binomial") #family binomial because response variable is 0 or 1. 
summary(m3)
m4 = glmer(every_zone ~ (1|trial), data = subset(df2, noon_day == 10), family ="binomial") #family binomial because response variable is 0 or 1. 
summary(m4)

aov(m3, m4)
write.table(summary(m3)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## get numbers per group
df2 %>% 
  filter(noon_day==10) %>% 
  group_by(strain_sex, every_zone) %>% 
  count()


