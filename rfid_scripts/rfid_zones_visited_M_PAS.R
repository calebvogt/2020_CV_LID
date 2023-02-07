## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)
library(glmmTMB) ## alternative package for linear models. 

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
pas_data <- as.data.frame(fread("data/priority_access_scores_males.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
pas_data <- pas_data %>% 
  filter(noon_day == 10)

rfid_data <- rfid_data %>% 
  filter(!(read_tag == "982.126057708741" & scan.date >= "07/24/2020")) %>%  #T005: Jeeves tag falls out of body on top of antenna)
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

df <- rfid_data %>% # 
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

## removes females
df3 <- merge(df, pas_data[,c("name", "csum_daily_capture_penalty")], by = "name")
df4 <- df3 %>% 
  mutate(terr_status = ifelse(csum_daily_capture_penalty>0, "high", "low"))

## rfid_unique_zones_visited_male_PAS_status.png
df5 <- df4 %>% 
  group_by(terr_status, noon_day) %>% 
  summarise(mean = mean(unique_zones_visited), sd = sd(unique_zones_visited),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df5, aes(x = noon_day, y = mean, color = terr_status)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_color_manual(breaks = c("high", "low"), values=c("red", "blue")) +
    theme_classic() +
    xlab("Day") +
    ylab("Unique zones visited") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.6,0.9))
          legend.position = "")
)
# ggsave(p, filename = "output/rfid_unique_zones_visited_male_PAS_status.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_unique_zones_visited_male_PAS_status.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
m0 = lmer(unique_zones_visited~noon_day+(1|name), data=df4)
m1 = lmer(unique_zones_visited~noon_day+(noon_day|name), data=df4) ## possible that different individuals have different slopes, thus (random slope | animal)
m2 = lmer(unique_zones_visited~terr_status*noon_day+(noon_day|name), data=df4)
m3 = glmmTMB::glmmTMB(unique_zones_visited~terr_status*noon_day + (noon_day|name)+(1|trial), data = df4) 
AIC(m0,m1,m2,m3) ## random slopes better. 
qqnorm(resid(m2))
qqline(resid(m2))
summary(m2)

write.table(summary(m2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


## exploratory plots
plot(unique_zones_visited~csum_daily_capture_penalty, subset(df3, noon_day.x == 10))
mean(subset(df3,noon_day.x==10&csum_daily_capture_penalty<0)$unique_zones_visited)
length(subset(df3, noon_day.x==10&csum_daily_capture_penalty<0)$unique_zones_visited)
sd(subset(df3, noon_day.x==10&csum_daily_capture_penalty<0)$unique_zones_visited)/sqrt(20)
mean(subset(df3, noon_day.x==10&csum_daily_capture_penalty>0)$unique_zones_visited)
length(subset(df3, noon_day.x==10&csum_daily_capture_penalty>0)$unique_zones_visited)
sd(subset(df3, noon_day.x==10&csum_daily_capture_penalty>0)$unique_zones_visited)/sqrt(46)
##
