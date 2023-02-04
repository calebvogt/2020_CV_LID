## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(lme4)
library(lmerTest)

metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

## Load displace data
displace <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT_GBI_displace.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# displace_daily_home_win_rate -----------------------------------------------
home_wins <- displace %>% 
  filter(!is.na(winner_zones_owned)) %>% ## only get home wins from mice who own zones
  group_by(winner, winner_loc, day) %>% ## get daily home wins
  tally() %>% 
  filter(winner_loc == "home") %>% 
  rename(home_wins = n, name = winner) %>% 
  na.omit()

home_losses <- displace %>% 
  filter(!is.na(loser_zones_owned)) %>% ## only get home losses from mice who own zones
  group_by(loser, loser_loc, day) %>% ## get daily home losses
  tally() %>% 
  filter(loser_loc == "home") %>% 
  rename(home_losses = n,  name = loser) %>% 
  na.omit()

df2 <- merge(home_wins, home_losses, by = c("name", "day"), all = T) ## per day

df3 <- df2 %>% 
  replace(is.na(.), 0) %>% ## change NAs to 0s, or mice that one at home but never lost at home. 
  mutate(home_win_rate = home_wins / (home_wins+home_losses)) %>% ## win rate = won events / all events
  merge(., metadata[, c("name", "trial", "strain")]) %>% 
  select(trial, strain, name, day, home_wins, home_losses, home_win_rate)

## Boxplot 
(p <- ggplot(df3, aes(x=as.factor(day), y=home_win_rate, fill = as.factor(strain))) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    xlab("Day") +
    ylab("Win rate at home") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "right") 
)
# ggsave(p, filename = "output/displace_home_win_rate_boxplot.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/displace_home_win_rate_boxplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## lineplot
df4 <- df3 %>% 
  group_by(strain, day) %>% 
  summarise(mean = mean(home_win_rate), sd = sd(home_win_rate),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df4, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
    # values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Home win-rate") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = "right")
)
# ggsave(p, filename = "output/displace_home_win_rate_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/displace_home_win_rate_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# displace_daily_away_win_rate -----------------------------------------------
away_wins <- displace %>% 
  filter(!is.na(winner_zones_owned)) %>% ## only get away wins from mice who own zones
  group_by(winner, winner_loc, day) %>% 
  tally() %>% 
  filter(winner_loc == "away") %>% 
  rename(away_wins = n, name = winner) %>% 
  na.omit()

away_losses <- displace %>% 
  filter(!is.na(loser_zones_owned)) %>% ## only get away losses from mice who own zones
  group_by(loser, loser_loc, day) %>% 
  tally() %>% 
  filter(loser_loc == "away") %>% 
  rename(away_losses = n,  name = loser) %>% 
  na.omit()

df2 <- merge(away_wins, away_losses, by = c("name", "day"), all = T)

df3 <- df2 %>% 
  replace(is.na(.), 0) %>% 
  mutate(away_win_rate = away_wins / (away_wins + away_losses)) %>% 
  merge(., metadata[, c("name", "trial", "strain")]) %>% 
  select(trial, strain, name, day, away_wins, away_losses, away_win_rate)

## Boxplot 
(p <- ggplot(df3, aes(x=as.factor(day), y=away_win_rate, fill = as.factor(strain))) + 
    geom_boxplot(outlier.shape=NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.2),
               pch=21) +
    xlab("Day") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "right") 
)
# ggsave(p, filename = "output/displace_away_win_rate_boxplot.png", device = "png", bg = "white") ## Change to F
# ggsave(p, filename = "output/displace_away_win_rate_boxplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## lineplot
df4 <- df3 %>% 
  group_by(strain, day) %>% 
  summarise(mean = mean(away_win_rate), sd = sd(away_win_rate),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df4, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
    # values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("away_win_rate") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = "right")
)
# ggsave(p, filename = "output/displace_away_win_rate_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/displace_away_win_rate_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# displace_overall_winrate_home_away ------------------------------------
## get home win rate
home_wins <- displace %>% 
  filter(!is.na(winner_zones_owned)) %>% ## only get home wins from mice who own zones
  group_by(winner, winner_loc) %>% ## get overall home wins
  tally() %>% 
  filter(winner_loc == "home") %>% 
  rename(home_wins = n, name = winner) %>% 
  na.omit()

home_losses <- displace %>% 
  filter(!is.na(loser_zones_owned)) %>% ## only get home losses from mice who own zones
  group_by(loser, loser_loc) %>% ## get all home losses
  tally() %>% 
  filter(loser_loc == "home") %>% 
  rename(home_losses = n,  name = loser) %>% 
  na.omit()

df2 <- merge(home_wins, home_losses, by = c("name"), all = T) ## all days

df3 <- df2 %>% 
  replace(is.na(.), 0) %>% ## change NAs to 0s, or mice that one at home but never lost at home. 
  mutate(home_win_rate = home_wins / (home_wins+home_losses)) %>% ## win rate = won events / all events
  merge(., metadata[, c("name", "trial", "strain")]) %>% 
  select(trial, strain, name, home_wins, home_losses, home_win_rate)

home_win_rate <- df3 %>% 
  mutate(type = "home") %>% 
  rename(win_rate = home_win_rate) %>% 
  select(name, type, win_rate)

## get away win rate
away_wins <- displace %>% 
  filter(!is.na(winner_zones_owned)) %>% ## only get away wins from mice who own zones
  group_by(winner, winner_loc) %>% 
  tally() %>% 
  filter(winner_loc == "away") %>% 
  rename(away_wins = n, name = winner) %>% 
  na.omit()

away_losses <- displace %>% 
  filter(!is.na(loser_zones_owned)) %>% ## only get away losses from mice who own zones
  group_by(loser, loser_loc) %>% 
  tally() %>% 
  filter(loser_loc == "away") %>% 
  rename(away_losses = n,  name = loser) %>% 
  na.omit()

df2 <- merge(away_wins, away_losses, by = c("name"), all = T)

df3 <- df2 %>% 
  replace(is.na(.), 0) %>% 
  mutate(away_win_rate = away_wins / (away_wins + away_losses)) %>% 
  merge(., metadata[, c("name", "trial", "strain")]) %>% 
  select(trial, strain, name, away_wins, away_losses, away_win_rate)

away_win_rate <- df3 %>% 
  mutate(type = "away") %>% 
  rename(win_rate = away_win_rate) %>% 
  select(name, type, win_rate)

## Territorial male win rates at home vs away 
df <- merge(home_win_rate, away_win_rate[, c("name", "win_rate")], by = c("name"), keep.all = T)
df4 <- rbind(home_win_rate, away_win_rate)
df5 <- merge(df4, metadata[,c("trial","strain", "name")], by = "name", keep=T)

## connected individual plot
df6 <- df5 %>% 
  group_by(name) %>%
  dplyr::filter(!n() == 1) %>% ## remove animals that do not have observed fights both at home and away. 
  ungroup()

(p <- ggplot(df6) +
    geom_point(aes(x=type,y=win_rate, col=strain), size = 1.5) +
    geom_line(aes(x=type,y=win_rate,group=name,col=strain), size = 0.75) + 
    scale_x_discrete(limits = c("home", "away")) +
    scale_color_manual(breaks = c("C57",  "NYOB"),
                        values=c("sienna", "skyblue4")) +
    ylab("Overall win rate for territory holding males") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = "")
)
# ggsave(p, filename = "output/displace_overall_win_rate_home_away.png", device = "png", bg = "white")
# ggsave(p, filename = "output/displace_overall_win_rate_home_away.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

## Stats
## LMM of overall win rates of territorial males at home and away. 
m2 <- lmer(win_rate~type+(1|name),data=df6)
summary(m2)
write.table(summary(m2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

## get mean and SE of win rates at home and away. 
df6 %>% 
  group_by(type) %>% 
  summarise(mean = mean(win_rate), sd = sd(win_rate),count = n(), se = (sd/(sqrt(count))))


