## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(readxl)
library(data.table)
library(gganimate)
library(ggrepel)
library(gifski)
library(av)
library(RColorBrewer)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")
df <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df1 <- df %>% 
  filter(trial == "T001" & sex == "M")

(p <- ggplot(df1, aes(x = field_time, y = zone, group = name, color = name)) + 
    ggtitle("T001") +
    geom_line() +
    geom_point() +
   scale_color_viridis_d(begin = 0.1, end = 0.8, option = "inferno") +
    labs(x = "Date", y = "Zone") +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
)
## To do: figure out how to add labels but just to final observed read. geom_text_repel tries to label all points, which crashes
# ggsave(p, filename = "output/rfid_zone_reads_grouped_lineplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_zone_reads_grouped_lineplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# RENDER THE ANIMATION
p2a <- p +
  geom_text_repel(aes(label = name), 
                  force =1, 
                  force_pull = 2) +
  geom_point(aes(group=seq_along(field_time))) +
  transition_reveal(field_time) 

# SAVE THE ANIMATION
animate(p2a, duration = 20, fps = 10, width = 800, height = 800, renderer = gifski_renderer()) #12 seconds per day at 15fps
anim_save("rfid_T001_male_zone_activity.gif", path = output_fp)

