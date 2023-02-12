## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(magrittr)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

el <- read.csv("data/ALLTRIAL_MOVEBOUT_edgelist.csv")
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)


# test <- df8 %>% 
#   group_by(trial, ID1, ID2,day) %>% 
#   summarize(sum_duration_s = sum(duration_s))
# top OppSex repeated top ranked partner -------------------------------------------------
df0 <- el %>% 
  merge(., metadata, by.x = "ID1", by.y = "name") %>%
  mutate(trial = trial.x, ID1_strain = strain, ID1_sex = sex, ID1_family_group = family_group, ID1_cage = cage, ID1_zone_drop = zone_drop) %>% 
  select(trial, day, sum_duration_s, ID1, ID2, ID1_strain, ID1_sex, ID1_family_group, ID1_cage, ID1_zone_drop) %>% 
  merge(., metadata, by.x = "ID2", by.y = "name") %>%
  mutate(trial = trial.x, ID2_strain = strain, ID2_sex = sex, ID2_family_group = family_group, ID2_cage = cage, ID2_zone_drop = zone_drop) %>% 
  select(trial, day, sum_duration_s, ID1, ID2, 
         ID1_strain, ID1_sex, ID1_family_group, ID1_cage, ID1_zone_drop,
         ID2_strain, ID2_sex, ID2_family_group, ID2_cage, ID2_zone_drop)


x <- unique(df0$ID1)
y <- unique(df0$ID2)
mice <- c(x,y)
mice <- unique(mice)

mouse_list <- list()
aa = mice[1]
aa = "Audrey"
for(aa in mice[1:length(mice)]) {
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
      select(trial, day, focal, partner, sum_duration_s) %>% 
      merge(., metadata, by.x = "focal", by.y = "name") %>%
      mutate(trial = trial.x, focal_sex = sex) %>% 
      select(trial, day, sum_duration_s, focal, focal_sex, partner) %>%  
      merge(., metadata, by.x = "partner", by.y = "name") %>%
      mutate(trial = trial.x, partner_sex = sex) %>% 
      select(trial, day, sum_duration_s, focal, focal_sex, partner, partner_sex) %>% 
      filter(!(focal_sex == partner_sex))
    
    if(nrow(df2)>0) {
      df3 <- df2  %>% 
        arrange(desc(sum_duration_s)) %>% 
        mutate(rank_order = c(1:nrow(.))) %>% 
        select(trial, day, focal, partner, rank_order, sum_duration_s)
      
      day_list[[bb]] <- df3
      print(paste(aa, "day", bb, "finshed"))
    }
    ## note that nas are produced for animals 
  }
  mouse_list[[aa]] <- do.call("rbind", day_list)
}
df <- do.call("rbind", mouse_list)

cols <- colnames(df)
df2 <- df %>% 
  filter(rank_order==1) %>% 
  group_by(focal) %>% 
  mutate(same_partner = ifelse(partner == lag(partner),1,0)) %>%
  merge(., metadata, by.x = "focal", by.y = "name") %>%
  mutate(strain_sex = paste(strain,sex,sep="-"), trial = trial.x) %>% 
  select(trial, day, focal, strain_sex, strain, sex, partner, rank_order, sum_duration_s, same_partner) %>% 
  filter(!(day==1)) %>% 
  filter(!is.na(same_partner))

df3 <- df2 %>% 
  group_by(strain_sex, day) %>% 
  summarise(count_tot=n(),
            p = mean(as.numeric(as.character(same_partner))),
            count_same=sum(same_partner)) %>% 
  mutate(perc_same = (count_same / count_tot) * 100,
         se = sqrt(p*(1-p)/count_tot)*100)


(p <- ggplot(df3, aes(x = day, y = perc_same, color = strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = perc_same - se, ymax = perc_same + se), width = 0.2) +
    # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    # values=c("red1", "red4", "steelblue", "steelblue4")) +
    # values=c("goldenrod1", "goldenrod4", "slateblue", "slateblue4")) +
    theme_classic() +
    ggtitle("OppSex") +
    xlab("Day") +
    ylab("% mice with repeated previously observed rank 1 OppSex partner") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          # legend.position = "",
          legend.title = element_blank(),
          legend.text = element_text(size=8))
)
ggsave(p, filename = "output/edge_data_perc_same_top_social_partner_OppSex.png", device = "png", bg = "white")
ggsave(p, filename = "output/edge_data_perc_same_top_social_partner_OppSex.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


## complete list of possible social partners, total of 19, and shuffle the list 10x times. for the cage mates, 
## get the same_cage percentage on shuffled data, do the same thing 1000x, and then plots histogram. for 
## for each animal 
