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


df <- rfid_data
ids <- unique(df$name)
ids <- ids[! ids %in% c('Carter')]
i=ids[1]
for(i in ids[14:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  
  # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  move_df <- filter(df, name == current_mouse)
  
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(move_df$trial), unique(move_df$drop), unique(move_df$sex), sep = "_")
  
  (p2 <- ggplot(move_df, aes(field_time, zone, color = factor(name))) + 
    geom_line(na.rm=TRUE, color="red", size=1) +
    geom_point() +
    ggtitle(paste(current_mouse_info, current_mouse, "zone_activity", sep = "_")) +
    scale_color_viridis_d() +
    labs(x = "Date", y = "Zone") +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
  )
  
  # RENDER THE ANIMATION
  p2a <- p2 +
    geom_point(aes(group=seq_along(field_time))) +
    transition_reveal(field_time)
  
  # SAVE THE ANIMATION
  animate(p2a, duration = 10, fps = 10, renderer = gifski_renderer())
  anim_save(paste("rfid_data",current_mouse_info, current_mouse, "zone_activity.gif", sep ="_"), path = output_fp)
}
