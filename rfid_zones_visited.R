## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(plotrix)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df0 <- rfid_data %>% 
  filter(!(read_tag == "982.126057708741" & scan.date >= "07/24/2020")) %>%  #T005: Jeeves tag falls out of body on top of antenna)
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

df <- df0 %>% # 
  filter(noon_day %in% 1:10) %>%
  mutate(strain_sex = paste(strain, sex, sep = "-")) %>% 
  select(trial, strain_sex, sex, strain, name, noon_day, zone) %>% 
  group_by(trial, strain_sex, sex,strain, name, noon_day, zone) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(trial, strain_sex, sex, strain,name, noon_day, .drop = FALSE) %>%
  tally() %>% # get number of unique zone visits
  complete(noon_day = 1:max(noon_day), fill = list(n = 0)) %>% # fill in days where mouse doesnt appear with 0s
  dplyr::rename(unique_zones_visited = n) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10))

# Summary --------------------------------------------------------------
## average zones visited, all animals
df1 <- df %>% 
  group_by(name) %>%
  summarise(mean_zone = mean(unique_zones_visited), sd_zone = sd(unique_zones_visited),count = n(), se_zone = (sd_zone/(sqrt(count))))
mean(df1$mean_zone)
plotrix::std.error(df1$mean_zone)
median(df1$mean_zone)
range(df1$mean_zone)

## average zones visited, males
df1 <- df %>% 
  filter(sex=="M") %>% 
  group_by(name) %>%
  summarise(mean_zone = mean(unique_zones_visited), sd_zone = sd(unique_zones_visited),count = n(), se_zone = (sd_zone/(sqrt(count))))

median(df1$mean_zone)
range(df1$mean_zone)


# All sex and strain groups -----------------------------------------------
df1 <- df %>% 
  group_by(strain_sex, noon_day) %>%
  summarise(mean_zone = mean(unique_zones_visited), sd_zone = sd(unique_zones_visited),count = n(), se_zone = (sd_zone/(sqrt(count))))

(p <- ggplot(df1, aes(x = noon_day, y = mean_zone, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_zone - se_zone, ymax = mean_zone + se_zone), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Zones visited") +
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
# ggsave(p, filename = "output/rfid_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS
df$trial <- as.factor(df$trial)
df$sex <- as.factor(df$sex)
df$strain <- as.factor(df$strain)
df$noon_day <- as.numeric(df$noon_day)
df$c57F <- ifelse(df$strain=="C57"&df$sex=="F",1,0)

m0 = lmer(unique_zones_visited ~ strain*sex*noon_day + (1|trial) + (noon_day|name), data = df)
summary(m0)

m1 = lmer(unique_zones_visited ~ strain*sex*log(noon_day) + (1|trial) + (log(noon_day)|name), data = df) 
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
emmeans(m1, pairwise ~ strain*sex*log(noon_day)) #tukey
ggpredict(m1, terms = c("strain", "sex", "noon_day"), type = "fe") %>% 
  plot()

m2 = lmer(unique_zones_visited ~ strain*sex+log(noon_day) + (1|trial) + (log(noon_day)|name), data = df)
summary(m2)

AIC(m0,m1,m2)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# Removing C57s demonstrates they are solely responsible for effect. 
m3 = lmer(unique_zones_visited ~ strain_sex*noon_day + (1|name), data = subset(df, df$strain_sex != "C57-F"))
summary(m3)
m4 = lmer(unique_zones_visited ~ c57F*log(noon_day) + (1|trial) + (1|name), data = df) # keeps C57 females in model, but evaluates them separately.
summary(m4)
m5 = lmer(unique_zones_visited ~ c57F*log(noon_day) + (1|trial) + (log(noon_day)|name), data = df) 
summary(m5)

AIC(m3,m4,m5)
write.table(summary(m5)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# daily models: change for days 1 through 10. 
d = lmer(unique_zones_visited ~ strain*sex + (1|trial), data = subset(df, noon_day == 2)) 
em.d <- emmeans(d, pairwise ~ strain*sex) # tukey adjusted
em.d$contrasts
write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

