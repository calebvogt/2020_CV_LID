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

# Expected FF time with cagemates  ---------------------
## females
df <- metadata %>% 
  filter(sex=="F") %>% 
  group_by(trial, strain, sex, cage) %>% 
  tally()

## C57
C57_mean_cagemates <- mean(c(4,4,4,4,4, ## T1
              4,4,4,4,4, 
              2,2,2, ## T2
              2,2,2,
              3,3,3,3,
              3,3,3,3, ## T3
              2,2,2,
              2,2,2,
              2,2,2, ## T6
              3,3,3,3,
              2,2,2))
C57_expected <- C57_mean_cagemates/9*100

## Outbred
OB_mean_cagemates <- mean(c(1,1, ## T4
                            3,3,3,3,
                            3,3,3,3,
                            1,1, ## T5
                            3,3,3,3,
                            1,1,
                            1,1,
                            1,1, ## T7
                            3,3,3,3,
                            3,3,3,3))
OB_expected <- OB_mean_cagemates/9*100


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
df2a <- df2 %>% 
  merge(., metadata[,c("name","trial", "cage", "num_cagemates")], by.x = "focal", by.y = "name") %>%
  select(trial,focal_strain_sex,focal,cage,num_cagemates,day,perc_cagemate_time) %>%
  mutate(perc_cagemate_time = signif(perc_cagemate_time, digits = 4)) %>% 
  unique() %>% 
  mutate(name_day = paste(focal,day,sep="_"))

## create data frame of daily expected # of females in the enclosure and the # of available cagemates. 
## triaging check: rae absent as expected, Isis only on day 1 & 2 as expected, rose absent on day 10 as expected
## rae = T003, cage 22, day 2-10, take away one. 
## rose = T003, cage 18, day 10, take away one
## isis = T004, cage 24, day 3-10, take away one
## add dashed lines with the null expectations each day. 
## strain average # of daily_num_cagemates / number of females in the enclosure on that day-1 * 100

df2b <- metadata %>% 
  filter(sex=="F") %>% 
  select(trial, strain, sex, name, code,cage,  num_cagemates) %>% 
  group_by(name) %>% 
  summarize(trial, strain, sex, name, code,cage, num_cagemates, day=rep(1:10, each=1)) %>% ## create idealized paddock population levels
  filter(!(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. 
         !(name == "Rose" & day >= 10)) %>% #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  ## adjust number of cagemates 
  mutate(adj_num_cagemates = ifelse(cage == 22 & day %in% 2:10, num_cagemates-1, num_cagemates), ## solve cage num for loss of rae
         adj_num_cagemates1 = ifelse(cage==18 & day==10, adj_num_cagemates-1,adj_num_cagemates), ## solve cage num for loss of rose
         adj_num_cagemates_final=ifelse(cage==24 & day %in% 3:10, adj_num_cagemates1-1,adj_num_cagemates1)) %>% ## solve cage num for loss of isis. 
  select(-adj_num_cagemates, - adj_num_cagemates1) %>% 
  ## adjust number of total females
  mutate(num_females = 10, ## num of females released into the trial
         adj_num_females=ifelse(trial=="T004" & day %in% 3:10, num_females-1,num_females), ## solve female num for loss of isis
         adj_num_females1=ifelse(trial=="T003" & day %in% 2:10, adj_num_females-1,adj_num_females), ## solve female num for loss of rae
         adj_num_females_final=ifelse(trial=="T003" & day == 10, adj_num_females1-1,adj_num_females1)) %>% ## solve female num for loss of rose
  select(-adj_num_females, - adj_num_females1) %>% 
  group_by(strain, day) %>% 
  mutate(daily_avg_num_cagemates = mean(adj_num_cagemates_final)) %>%  ## strain level mean
  ungroup() %>% 
  mutate(expected_perc_cagemate_time = daily_avg_num_cagemates/(adj_num_females_final-1)*100) %>%
  select(trial, strain, name, day, expected_perc_cagemate_time) %>% 
  mutate(name_day = paste(name,day,sep="_"))
  

df3 <- df2a %>% 
  merge(., df2b[,c("name_day", "expected_perc_cagemate_time"), by = "name_day"]) %>% 
  select(-name_day) %>% 
  mutate(delta = perc_cagemate_time - expected_perc_cagemate_time)
  
## grouped values, line plot. Can also change delta directly to perc_cagemate_time
df4 <- df3 %>% 
  group_by(focal_strain_sex, day) %>%
  summarise(mean = mean(delta), sd = sd(delta),count = n(), se = (sd/(sqrt(count))))

# ## Boxplot of raw values
# (p <- ggplot(df4, aes(x=as.factor(day), y=perc_cagemate_time, fill = as.factor(focal_strain_sex))) + 
#     geom_boxplot(outlier.shape=NA) +
#     geom_point(position=position_jitterdodge(jitter.width = 0.2),
#                pch=21) +
#     xlab("Day") +
#     ylab("daily_perc_FF_cagemate_time") +
#     theme_classic() +
#     theme(axis.text.x = element_text(color = "black", size = 8),
#           axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
#           axis.text.y = element_text(color = "black", size = 8),
#           axis.title.y = element_text(color = "black", size = 8, face = "bold"),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           panel.background = element_rect(fill = "transparent"), 
#           legend.title = element_blank(),
#           legend.text = element_text(size=8), 
#           legend.position = "right") 
# )
# # ggsave(p, filename = "output/edge_perc_FF_time_cagemates_boxplot.png", device = "png", bg = "white") ## Change to F
# # ggsave(p, filename = "output/edge_perc_FF_time_cagemates_boxplot.png", device = "svg", width=2.5, height=2.15, bg = "transparent")


# ## Individual values, lineplot
# (p <- ggplot(df4, aes(x=as.factor(day), y=perc_cagemate_time,group = focal, color = focal_strain_sex)) + 
#     geom_line(size =1, alpha = 0.5) +
#     xlab("Day") +
#     ylab("daily_perc_FF_cagemate_time") +
#     theme_classic() +
#     theme(axis.text.x = element_text(color = "black", size = 8),
#           axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
#           axis.text.y = element_text(color = "black", size = 8),
#           axis.title.y = element_text(color = "black", size = 8, face = "bold"),
#           plot.background = element_rect(fill = "transparent", color = NA),
#           panel.background = element_rect(fill = "transparent"), 
#           legend.title = element_blank(),
#           legend.text = element_text(size=8), 
#           legend.position = "right") 
# )
# # ggsave(p, filename = "output/edge_perc_FF_time_cagemates_lineplot_individual.png", device = "png", bg = "white") ## Change to F
# # ggsave(p, filename = "output/edge_perc_FF_time_cagemates_lineplot_individual.png", device = "svg", width=2.5, height=2.15, bg = "transparent")

## grouped lineplot, % female cagemate preference
(p <- ggplot(df4, aes(x = day, y = mean, color = focal_strain_sex)) +
    geom_line(size = 0.75) + 
    geom_hline(yintercept = 0, linetype="dashed") + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab(paste("% female cagemate preference \n (Observed-Expected)")) +
    # labs(y=paste0("% female-female social time","\n", "spent with cagemates")) +
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
ggsave(p, filename = "output/edge_perc_FF_cagemate_preference_lineplot.png", device = "png", bg = "white")
ggsave(p, filename = "output/edge_perc_FF_cagemate_preference_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats % female preference obs-expected

## C57 comp to OB and 0
m8<-lmer(delta~focal_strain_sex+(1|focal)+(1|trial), data=filter(df3, day %in% 2:10))
AIC(m8)
summary(m8)
write.table(summary(m8)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## NYOB comp to 0
df3$nyob_first<-factor(df3$focal_strain_sex, levels=c('NYOB-F', 'C57-F')) 
m9<-lmer(delta~nyob_first+(1|focal)+(1|trial), data=filter(df3, day %in% 2:10))
summary(m9)
write.table(summary(m9)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## Stats perc_cagemate_time_deprecated
m1<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*day+(1|focal), data=df3)
m2<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*log(day)+(1|focal), data=df3)
m3<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*log(day)+(log(day)|focal), data=df3)
m4<-lmer(asin(sqrt(perc_cagemate_time/100))~focal_strain_sex*log(day)+(log(day)|focal)+(1|trial), data=df3) 

qqnorm(resid(m4))
qqline(resid(m4))
AIC(m1,m2,m3,m4)
summary(m4)
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)







