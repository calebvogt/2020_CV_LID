## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
rfid_data <- rfid_data %>% 
  filter(!(read_tag == "982.126057708741" & scan.date >= "07/24/2020")) %>%  #T005: Jeeves tag falls out of body on top of antenna)
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
gc()

df <- rfid_data
ids <- unique(df$name)
data_list <- list()
aa = ids[60]
for(aa in ids[1:length(ids)]){
  df1 <- df %>% 
    filter(name == aa) 
  
  # delete consecutive repeat antenna hits, only keep rows where zones change. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$zone)] 
  df2$dist <- NA
  n <- nrow(df2)
  if(n==1){
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  } else{
    df2$dist[2:n] <- sqrt((df2$zone_x[2:n] - df2$zone_x[1:n-1])^2 + (df2$zone_y[2:n] - df2$zone_y[1:n-1])^2)
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  }
}
df3 <- do.call("rbind", data_list)

df4 <- df3 %>% 
  group_by(trial, strain_sex,strain, sex, name, noon_day) %>% 
  tally(sum(dist)) %>% 
  complete(noon_day = 1:10, fill = list(n = 0)) %>%  # fill in days where mouse doesnt appear with 0s.
  dplyr::rename(dist = n) %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
 
df5 <- df4 %>% 
  group_by(strain_sex, noon_day) %>%
  summarise(mean_n = mean(dist), sd_n = sd(dist),count = n(), se_n = (sd_n/(sqrt(count))))
  
(p <- ggplot(df5, aes(x = noon_day, y = mean_n, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    # values=c("red1", "red4", "steelblue", "steelblue4")) +
    # values=c("goldenrod1", "goldenrod4", "slateblue", "slateblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Minimum distance travelled (m)") +
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
ggsave(p, filename = "output/rfid_distance_travelled.png", device = "png", bg = "white")
ggsave(p, filename = "output/rfid_distance_travelled.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS  
df7 <- df4 %>% 
  group_by(trial, strain, sex, name, noon_day) %>% 
  tally(sum(dist)) %>% 
  group_by(trial, strain, sex, name, noon_day) %>% 
  dplyr::rename(sum_dist = n)

# write.table(df7, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
# df7 <- read_excel("Figure_Data.xlsx", sheet = "Fig2a")

df7$trial <- as.factor(df7$trial)
df7$strain <- as.factor(df7$strain)
df7$sex <- as.factor(df7$sex)
df7$name <- as.factor(df7$name)
df7$c57_F<-ifelse(df7$strain!="C57"|df7$sex!="F","other","pC57_F")

mod1 = lmer(sum_dist ~ strain*sex*log(noon_day) + (1|trial) + (log(noon_day)|name), data = df7) 
mod2 = lmer(sum_dist ~ strain*sex*log(noon_day) + (log(noon_day)|name), data = df7) 
mod3 = lmer(sum_dist ~ strain*sex*log(noon_day) + (1|trial) + (1|name), data = df7) 
mod4 = lmer(sum_dist ~ strain+sex+log(noon_day) + (1|trial) + (log(noon_day)|name), data = df7) 
AIC(mod1,mod2, mod3)
qqnorm(resid(mod1))
anova(mod1)
summary(mod1)
emmeans(mod1, pairwise ~ strain*sex*log(noon_day))
write.table(summary(mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# daily models: change for days 1 through 10 and report
d = lmer(sum_dist ~ strain*sex + (1|trial), data = subset(df7, noon_day == 10)) 
em.d <- emmeans(d, pairwise ~ strain*sex) #tukey adjusted
write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
  
## Total distance travelled
df6a <- df4 %>% 
  group_by(trial, strain_sex, name) %>% 
  tally(sum(dist)) %>% 
  dplyr::rename(total_dist = n)

# write.table(df6a, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)

(p <- ggplot(df6a, aes(x=strain_sex, y=total_dist, fill = strain_sex)) + 
    geom_violin(width=1, alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    scale_y_continuous(limits = c(0,2000)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    # scale_fill_discrete(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                        # values=c("sienna1", "sienna", "skyblue", "skyblue4"),
                        # labels=c('C57 Female', 'C57 Male', 'Outbred Female', 'Outbred Male')) +
    xlab("") +
    ylab("Total estimated dist. travelled (m)") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.position = "none")
)
ggsave(p, filename = "output/rfid_distance_travelled_total.png", device = "png", bg = "white")
ggsave(p, filename = "output/rfid_distance_travelled_total.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

m1 <- lmer(total_dist ~ strain_sex + (1|trial), data = df6a)
summary(m1)
anova(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

em.d <- emmeans(m1, pairwise ~ strain_sex, adjust = 'none') 
write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
