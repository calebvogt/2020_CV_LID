## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(readxl)
library(gganimate)
library(gifski)
library(av)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")
rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
rfid_data <- rfid_data %>% 
  filter(noon_day %in% 1:12) %>% # traps dropped on 5/19 or day 13. so data up until day 12 is good. 
  mutate(sex_phase = paste(sex, phase, sep = "_"))


# GIF: 12 HOUR ALL TRIAL MALE PADDOCK ACTIVITY -----------------------
trial <- unique(clean$Trial)

## FILTER BY SEX, 20 GRIDS IS TOO MUCH. 
df <- clean %>% 
  filter(Sex == "M")

## LOOP ACROSS TRIALS
i=trial[4]
for(i in trial[1:length(trial)]) {
  current_trial <- print(i)
  
  move_df <- subset(df, Trial == current_trial)
  
  sex <- unique(df$Sex)
  
  # PLOT 4: ALL INDIVIDUAL MOVMENT ANIMATION ANIMATION
  p4 <- ggplot(move_df, aes(zone_x, zone_y, color = Name)) +
    ggtitle("Group Movement") +
    geom_point(show.legend = TRUE, alpha = 0.7, size = 2) +
    xlim("A","B") +
    ylim("A","B","C","D") +
    geom_jitter(width = 0.1, height = 0.1) +
    scale_color_viridis_d() +
    labs(x = "none", y = "none") +
    facet_wrap(~Name) +
    theme(plot.background = element_blank(),
          panel.grid.major = element_line(colour="black", size = 0.5),
          panel.grid.minor = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 20))
  
  # RENDER THE ANIMATION
  p4a <- p4 +
    transition_time(Field_Time) +
    labs(subtitle = "Time: {frame_time}")
  
  # SAVE THE ANIMATION
  animate(p4a, duration = 60, fps = 24, width = 500, height = 500, renderer = gifski_renderer())
  anim_save(paste(current_trial,sex, "Group_Movement_FULL.gif", sep ="_"), path = output_fp)
}



