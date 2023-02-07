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

# edge_perc_FF_time_cagemates -------------------------------------------------
## group dyad spatiotemporal overlap by day
df <- data %>%
  group_by(trial, ID1, ID2,day) %>%
  summarize(sum_duration_s = sum(duration_s))

df0 <- df %>% 
  merge(., metadata, by.x = "ID1", by.y = "name") %>%
  mutate(trial = trial.x, ID1_strain = strain, ID1_sex = sex, ID1_family_group = family_group, ID1_cage = cage, ID1_zone_drop = zone_drop) %>% 
  select(trial, day, sum_duration_s, ID1, ID2, ID1_strain, ID1_sex, ID1_family_group, ID1_cage, ID1_zone_drop) %>% 
  merge(., metadata, by.x = "ID2", by.y = "name") %>%
  mutate(trial = trial.x, ID2_strain = strain, ID2_sex = sex, ID2_family_group = family_group, ID2_cage = cage, ID2_zone_drop = zone_drop) %>% 
  select(trial, day, sum_duration_s, ID1, ID2, 
         ID1_strain, ID1_sex, ID1_family_group, ID1_cage, ID1_zone_drop,
         ID2_strain, ID2_sex, ID2_family_group, ID2_cage, ID2_zone_drop) %>% 
  filter(ID1_sex == "F" & ID2_sex == "F") ## choose female subjects

x <- unique(df0$ID1)
y <- unique(df0$ID2)
mice <- c(x,y)
mice <- unique(mice)
mouse_list <- list()
aa = mice[1]
for(aa in mice[1:length(mice)]) { ## for each mouse, rank their top social partners per day
  df <- df0 %>% 
    filter(ID1 == aa | ID2 == aa)
  days <- unique(df$day)
  
  day_list <- list()
  bb =1
  for(bb in days[1:length(days)]) {
    df2 <- df %>% 
      filter(ID1 == aa | ID2 == aa) %>% 
      filter(day == bb) %>%
      mutate(focal = aa, ID1 = na_if(ID1, aa), ID2 = na_if(ID2, aa)) %>% 
      mutate(partner = coalesce(ID1,ID2)) %>% 
      arrange(desc(sum_duration_s)) %>% 
      mutate(rank_order = c(1:nrow(.))) %>% 
      select(trial, day, focal, partner, rank_order, sum_duration_s)
    
    day_list[[bb]] <- df2
    print(paste(aa, "day", bb, "finshed"))
  }
  mouse_list[[aa]] <- do.call("rbind", day_list)
}
df <- do.call("rbind", mouse_list)
cols <- colnames(df)

## determine percent time with cage mates vs non cagemates
df2 <- df %>% 
  merge(., metadata, by.x = "focal", by.y = "name") %>%
  mutate(trial = trial.x, focal_strain = strain, focal_sex = sex, focal_family_group = family_group, focal_cage = cage, focal_zone_drop = zone_drop) %>% 
  select(trial, day, rank_order,sum_duration_s, focal, partner, focal_strain, focal_sex, focal_family_group, focal_cage, focal_zone_drop) %>%  
  merge(., metadata, by.x = "partner", by.y = "name") %>%
  mutate(trial = trial.x, 
         partner_strain = strain, 
         partner_sex = sex, 
         partner_family_group = family_group, 
         partner_cage = cage, 
         partner_zone_drop = zone_drop) %>% 
  select(trial, focal, partner, day,rank_order, sum_duration_s,focal_strain, focal_sex, focal_family_group, focal_cage, focal_zone_drop,
         partner_strain, partner_sex, partner_family_group, partner_cage, partner_zone_drop) %>%  
  arrange(focal, day) %>% 
  group_by(focal,day) %>% 
  mutate(same_cage = ifelse(focal_cage == partner_cage,"Y","N"), 
         focal_strain_sex = paste(focal_strain,focal_sex,sep="-")) %>% 
  mutate(total_daily_time = sum(sum_duration_s)) %>% 
  relocate(same_cage, total_daily_time, .after = sum_duration_s) %>% 
  select(focal, partner, focal_strain_sex, day, same_cage, sum_duration_s, total_daily_time) %>% 
  group_by(focal,day, same_cage) %>% 
  mutate(total_cage_noncage_time = sum(sum_duration_s),
        perc_cage_noncage_time = (total_cage_noncage_time / total_daily_time)*100) %>% 
  mutate(perc_cagemate_time = ifelse(same_cage == "Y", perc_cage_noncage_time, abs(perc_cage_noncage_time-100))) %>% 
  ungroup()

## isolate desired variables. note each female has their own percent time spent with cagemates
df3 <- df2 %>% 
  merge(., metadata[,c("name","trial", "cage", "num_cagemates")], by.x = "focal", by.y = "name") %>%
  select(trial,focal_strain_sex,focal,cage,num_cagemates,day,perc_cagemate_time) %>%
  mutate(perc_cagemate_time = signif(perc_cagemate_time, digits = 4)) %>% 
  unique() %>% 
  mutate(name_day = paste(focal,day,sep="_"))

## grouped values, line plot. Can also change delta directly to perc_cagemate_time
df4 <- df3 %>% 
  group_by(focal_strain_sex, day) %>%
  summarise(mean = mean(perc_cagemate_time), sd = sd(perc_cagemate_time),count = n(), se = (sd/(sqrt(count))))

## grouped lineplot, perc_FF_cagemate_time
(p <- ggplot(df4, aes(x = day, y = mean, color = focal_strain_sex)) +
    geom_line(size = 0.75) + 
    # geom_hline(yintercept = 0, linetype="dashed") + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    labs(y=paste0("% female-female social time","\n", "spent with cagemates")) +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.41,0.8))
          legend.position = "")
)
ggsave(p, filename = "output/edge_perc_FF_cagemate_time_lineplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/edge_perc_FF_cagemate_time_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats perc_cagemate_time
m1<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*day+(1|focal), data=df3)
m2<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*log(day)+(1|focal), data=df3)
m3<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*log(day)+(log(day)|focal), data=df3)
m4<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*log(day)+(log(day)|focal)+(1|trial), data=df3) 

qqnorm(resid(m4))
qqline(resid(m4))
AIC(m1,m2,m3,m4)
summary(m4)
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

