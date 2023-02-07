## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(lubridate)
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

# options(scipen = 3)
data <- df0 %>% 
  filter(sex=="M") %>% 
  group_by(trial, zone, noon_day) %>% 
  mutate(total_zone_reads=n()) %>% 
  group_by(trial, strain, name, zone, noon_day, total_zone_reads) %>% 
  tally() %>% 
  rename(mus_zone_reads = n) %>% 
  mutate(mus_prop_totalzone_reads = (mus_zone_reads/total_zone_reads)) %>% 
  group_by(trial, zone, noon_day) %>% 
  mutate(rank_order = rank(-mus_prop_totalzone_reads)) %>% ## create rank order
  mutate(trial_zone_day = paste(trial,zone, noon_day, sep = "_")) %>% 
  arrange(trial_zone_day, rank_order)


# % male-sourced reads within each zone belonging to top single male on Day 3 --------
## pools C57 and outbred males. 
df <- data %>% 
  filter(rank_order==1, noon_day == 3) %>% ## 
  group_by(noon_day) %>% 
  summarize(mean(mus_prop_totalzone_reads),
            sd = sd(mus_prop_totalzone_reads),
            count = n(),
            se = (sd/sqrt(count)))
df


# % male-sourced reads within each zone belonging to top single male by strain --------
df <- data %>% 
  filter(rank_order==1) %>% ## 
  mutate(trial_zone = paste(trial, zone, sep="_"))

## boxplot
(p <- ggplot(df, aes(x=as.factor(noon_day), y=mus_prop_totalzone_reads, fill = as.factor(strain))) +
    geom_boxplot()
)

## lineplot of percentages
df2 <- df %>% 
  group_by(strain, noon_day) %>%
  summarise(mean = mean(mus_prop_totalzone_reads*100), 
            sd = sd(mus_prop_totalzone_reads*100),
            count = n(), 
            se = (sd/(sqrt(count))), 
            mode = mode(mus_prop_totalzone_reads*100))

(p <- ggplot(df2, aes(x = noon_day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    xlab("Day") +
    ylab("% male-sourced zone RFID reads from rank 1 male") + 
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_y_continuous(limits = c(70,100), breaks = seq(70, 100, by = 10)) +
    scale_color_manual(breaks = c("C57", "NYOB"),
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
# ggsave(p, filename = "output/rfid_perc_malesourced_zone_reads_rank_1_male.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_perc_malesourced_zone_reads_rank_1_male.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
df$trial <- as.factor(df$trial)
df$strain <- as.factor(df$strain)
df$noon_day <- as.numeric(df$noon_day)

m1 = lmer(asin(sqrt(mus_prop_totalzone_reads)) ~ strain*log(noon_day) + (1|name), data = df) 
AIC(m1)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


