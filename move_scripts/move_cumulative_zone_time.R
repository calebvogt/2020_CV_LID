## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# sapply(df, function(x) sum(is.na(x)))
df <- move_data %>% 
  dplyr::select(trial, strain, sex, name, noon_day,duration_s) %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(trial,strain, sex, name, noon_day) %>% 
  tally(sum(duration_min)) %>% 
  mutate(csum_zone_min = cumsum(n)) %>%  #get cumulative # of novel mice met
  ungroup() %>% 
  complete(name, noon_day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, noon_day) %>%
  fill(trial,strain, sex, n,csum_zone_min) %>% ## fill cumulative sum data from last observed day
  mutate(group = paste0(strain, "-", sex)) %>% 
  dplyr::select(trial,group, strain, sex, name, noon_day, csum_zone_min)

# data cleaning. 
df <- df %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 

# output stats df1
# write.table(df, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)

## Graph
df1 <- df %>% 
  group_by(group, noon_day) %>% 
  summarise(mean_n = mean(csum_zone_min), sd_n = sd(csum_zone_min),count = n(), se_n = (sd_n/(sqrt(count))))

# PLOT
(p <- ggplot(df1, aes(x=noon_day, y=mean_n, color = group)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(breaks = seq(1,11,by=1), limits = c(1,10.2)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day") +
    ylab("Cumulative zone time (min)") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          # legend.position = c(0.15,0.8))
          legend.position = "") +
    guides(color=guide_legend(override.aes=list(fill=NA)))
)
# ggsave(p, filename = "output/move_cumulative_zone_time.png", device = "png", bg = "white")
# ggsave(p, filename = "output/move_cumulative_zone_time.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# STATS
summary(df)
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)

mod1 = lmer(csum_zone_min ~ sex*strain + (1|trial), data = subset(df, noon_day == 10))
anova(mod1)
summary(mod1)
write.table(summary(mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
# 
# em.mod1 = emmeans(mod1, pairwise ~ sex*strain)
# write.table(em.mod1$emmeans, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
# write.table(em.mod1$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)



