## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(readxl)
library(scales)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")
rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
rfid_data <- rfid_data %>% 
  filter(noon_day %in% 1:10)


# T001
rfid_data %>% filter(trial == "T001", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T001_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T001", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T001_F_zone_use.png", device = "png", bg = "white")


# T002
rfid_data %>% filter(trial == "T002", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T002_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T002", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T002_F_zone_use.png", device = "png", bg = "white")


# T003
rfid_data %>% filter(trial == "T003", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T003_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T003", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T003_F_zone_use.png", device = "png", bg = "white")

# T004
rfid_data %>% filter(trial == "T004", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T004_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T004", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T004_F_zone_use.png", device = "png", bg = "white")


# T005
rfid_data %>% filter(trial == "T005", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T005_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T005", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T005_F_zone_use.png", device = "png", bg = "white")

# T006
rfid_data %>% filter(trial == "T006", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T006_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T006", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T006_F_zone_use.png", device = "png", bg = "white")


# T007
rfid_data %>% filter(trial == "T007", sex == "M") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name)
ggsave(filename = "output/rfid_data_T007_M_zone_use.png", device = "png", bg = "white")

rfid_data %>% filter(trial == "T007", sex == "F") %>% 
  ggplot(., aes(x = field_time, y = zone)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
  facet_wrap(~name, scales = "free")
ggsave(filename = "output/rfid_data_T007_F_zone_use.png", device = "png", bg = "white")




# individual plots
ids <- unique(df$name)
i=ids[20]
for(i in ids[1:length(ids)]) {
  current_mouse <- print(i)
  move_df <- subset(df, name == current_mouse)
  current_mouse_info <- paste(unique(move_df$trial), unique(move_df$drop), unique(move_df$sex), sep = " ")
  
  p <- ggplot(move_df) +
    aes(x = field_time, y = zone, color = reader.id) +
    geom_point(na.rm=TRUE, size=1) +
    ggtitle(paste(current_mouse_info, current_mouse, sep = " ")) +
    xlab("Date") + 
    ylab("Zone") +
    scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
  
  p
  
  ggsave(filename=paste("rfid_zonereads",current_mouse_info, current_mouse, sep = "_"), 
         plot=p, 
         width = 5, 
         height = 4, 
         dpi = 300, 
         units = "in", 
         device='png', 
         path = output_fp)
} 

