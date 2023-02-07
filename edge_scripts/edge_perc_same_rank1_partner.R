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

# top repeated ranked partner list for all partners -------------------------------------------------
x <- unique(el$ID1)
y <- unique(el$ID2)
mice <- c(x,y)
mice <- unique(mice)

mouse_list <- list()
aa = mice[1]
aa = "Roach"
for(aa in mice[1:length(mice)]) {
  df <- el %>% 
    filter(ID1 == aa | ID2 == aa)
  days <- unique(df$day)
  
  day_list <- list()
  bb =6
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
df2 <- df %>% 
  merge(., metadata, by.x = "focal", by.y = "name") %>%
  mutate(trial = trial.x, focal_strain = strain, focal_sex = sex, focal_family_group = family_group, focal_cage = cage, focal_zone_drop = zone_drop) %>% 
  select(trial, day, rank_order,sum_duration_s, focal, partner, focal_strain, focal_sex, focal_family_group, focal_cage, focal_zone_drop) %>%  
  merge(., metadata, by.x = "partner", by.y = "name") %>%
  mutate(trial = trial.x, partner_strain = strain, partner_sex = sex, partner_family_group = family_group, partner_cage = cage, partner_zone_drop = zone_drop) %>% 
  select(trial, focal, partner, day,rank_order, sum_duration_s,focal_strain, focal_sex, focal_family_group, focal_cage, focal_zone_drop,
         partner_strain, partner_sex, partner_family_group, partner_cage, partner_zone_drop) %>%  
  filter(rank_order==1) %>% 
  group_by(focal) %>% 
  arrange(focal, day) %>% 
  mutate(same_partner = ifelse(partner == lag(partner),1,0), ## is this rank 1 partner the same as the previous day rank 1 partner?
         same_cage = ifelse(focal_cage == partner_cage,1,0)) %>% ## is this rank 1 partner a cage mate?
  mutate(focal_strain_sex = paste(focal_strain,focal_sex,sep="-")) %>% 
  relocate(focal_strain_sex,same_partner, same_cage, .after = sum_duration_s) %>% 
  filter(!(day==1))

## for each strain sex, calculate the percent rank 1 partners who are the same as the last seen rank 1 partner, per day
df3 <- df2 %>% 
  group_by(focal_strain_sex, day) %>% 
  summarise(count_tot=n(),
            p = mean(as.numeric(as.character(same_partner))),
            count_same=sum(same_partner)) %>%
  mutate(perc_same = (count_same / count_tot) * 100,
         se = sqrt(p*(1-p)/count_tot)*100)

(p <- ggplot(df3, aes(x = day, y = perc_same, color = focal_strain_sex)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = perc_same - se, ymax = perc_same + se), width = 0.3) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% mice with repeated last top ranked partner") + ## Change
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
ggsave(p, filename = "output/edge_data_perc_same_top_social_partner.png", device = "png", bg = "white")
ggsave(p, filename = "output/edge_data_perc_same_top_social_partner.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# same_partner stats
# df2$trial <- as.factor(df2$trial)
# df2$sex <- as.factor(df2$sex)
# df2$strain <- as.factor(df2$strain)
# df2$noon_day <- as.numeric(df2$noon_day)


## These models need to be evaluated... 
m1 = glmer(same_partner ~ focal_strain_sex*log(day)  +(1|trial)+ (1+log(day)|focal), family="binomial", data=df2) 
m2 = glmer(same_partner ~ focal_strain_sex*log(day)  + (1|trial) + (1|focal), family="binomial", data=df2) 
# when we experienced convergence issues we simplified the random effect structure until the convergence disappeared. 

m3 = glmer(same_partner ~ focal_strain+log(day)  +(1|trial)+ (1+log(day)|focal), family="binomial", data=df2) # strain only because males and females are not independent
## MM,FF, OppSex models with just focal_strain


AIC(m1,m2)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m3)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
