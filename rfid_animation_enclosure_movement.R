## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(gganimate)
library(gifski)
library(av)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

## Change zone xy coordinates. change this mess....replace in dplyr or something
rfid_data$zone_x <- ifelse(grepl("1", rfid_data$zone), 1, 
                    ifelse(grepl("2", rfid_data$zone), 3,
                           ifelse(grepl("3", rfid_data$zone), 1,
                                  ifelse(grepl("4", rfid_data$zone), 3,
                                         ifelse(grepl("5", rfid_data$zone), 1,
                                                ifelse(grepl("6", rfid_data$zone), 3,
                                                       ifelse(grepl("7", rfid_data$zone), 1,
                                                              ifelse(grepl("8", rfid_data$zone), 3,
                                                                     "none"))))))))

rfid_data$zone_y <- ifelse(grepl("1", rfid_data$zone), 1, 
                    ifelse(grepl("2", rfid_data$zone), 1,
                           ifelse(grepl("3", rfid_data$zone), 3,
                                  ifelse(grepl("4", rfid_data$zone), 3,
                                         ifelse(grepl("5", rfid_data$zone), 5,
                                                ifelse(grepl("6", rfid_data$zone), 5,
                                                       ifelse(grepl("7", rfid_data$zone), 7,
                                                              ifelse(grepl("8", rfid_data$zone), 7,
                                                                     "none"))))))))

rfid_data$zone_x <- as.numeric(rfid_data$zone_x)
rfid_data$zone_y <- as.numeric(rfid_data$zone_y)

# choose who is included in the animation
df <- rfid_data %>% 
  # filter(name=="Abrams") %>% 
  # filter(trial == "T001") %>% 
  # filter(sex == "M") %>% 
  # filter(noon_day %in% 1:5) %>% 
  mutate(zone_trans = zone) %>% 
  select(trial, strain,name, code, sex, time_sec, zone, zone_x, zone_y) %>% 
  arrange(name, time_sec) %>% 
  mutate(time = 1:nrow(.)) %>% 
  mutate(zone_xj = jitter(zone_x, factor = 0.5), 
         zone_yj = jitter(zone_y, factor = 0.5)) %>% 
  # mutate(next_zone_xj = lead(zone_xj),
  # next_zone_yj = lead(zone_yj)) %>% ## created for segment code attempt. drop. 
  mutate(time_day = time_sec/86400+1) %>% 
  mutate(time_sec_round = round(time_sec))

# get data of xy boxes points and define boxes
# box0 <- unique(rfid_data[c("zone","zone_x","zone_y")])
# box <- box0 %>% 
#   mutate(x1=zone_x+1, x2=zone_x-1, y1=zone_y+1, y2=zone_y-1) 

## pull transitions from participants
# ids <- unique(df$name)
# data_list <- list()
# aa = ids[1]
# for(aa in ids[1:length(ids)]){
#   df1 <- df %>% 
#     filter(name == aa) 
#   data_list[[aa]] <- as.data.table(df1)[, .SD[1], by = rleid(df1$zone_trans)] # delete consecutive repeat antenna hits, only keep rows where zones change. 
# }
# df2 <- do.call("rbind", data_list)

# df3 <- df %>%  ## changed temp from df2 to df
#   select(name, sex, code, time_sec, zone, zone_x, zone_y) %>% 
#   arrange(name, time_sec) %>% 
#   mutate(time = 1:nrow(.)) %>% 
#   mutate(zone_xj = jitter(zone_x, factor = 0.5), 
#          zone_yj = jitter(zone_y, factor = 0.5)) %>% 
#   # mutate(next_zone_xj = lead(zone_xj),
#          # next_zone_yj = lead(zone_yj)) %>% ## created for segment code attempt. drop. 
#   mutate(time_day = time_sec/86400+1) %>% 
#   mutate(time_sec_round = round(time_sec))
#   
# testing -----------------------------------------------------------------
## loop for individuals
ids <- unique(df$code)
ids
i=ids[1]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  
  # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  df1 <- subset(df, code == current_mouse)
 
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(df1$trial), unique(df1$strain), unique(df1$sex), unique(df1$code), sep = "_")
  current_mouse_info <- current_mouse_info[current_mouse_info != "NA_NA_NA_NA"]
  print(current_mouse_info)
  
  p <- ggplot(data=df1, aes(x=zone_xj, y=zone_yj, color=code, label=code, shape =sex)) +
    # geom_rect(aes(xmin = zone_x+1, xmax = zone_x-1, ymin = zone_y+1, ymax = zone_y-1)) + ## draw rectangles. 
    geom_point(size=4) + 
    geom_path(size=1,alpha =0.1) +
    geom_text(nudge_x = 0.2)+
    # geom_path(size=1.5) +
    # xlim(0,15.24) +
    # ylim(0,38.1) +
    geom_vline(aes(xintercept=2)) + 
    geom_hline(aes(yintercept=2)) +
    geom_hline(aes(yintercept=4)) +
    geom_hline(aes(yintercept=6)) +
    scale_x_continuous(breaks = seq(1,4,1), limits=c(0,4)) +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(0,8)) +
    transition_reveal(along=time_sec_round) + ## has to be an integer value.
    # transition_manual(along=time_sec_round, time_sec_round_fade) +
    labs(title = 'Day: {round(frame_along/86400+1,digits=2)}') +
    labs(x = "none", y = "none") +
    theme(plot.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
  animate(p,duration=20, fps = 20, height = 400, width=200, renderer=gifski_renderer())
  anim_save(paste("output/animations/",current_mouse_info, "_enclosure_movement.gif", sep =""))
}


# create animation --------------------------------------------------------

## create background only
(p <- ggplot() +
    geom_rect(data=box, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0.5) +
    # geom_text(data=box, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=antenna), size=4, vjust = -2.2) + 
    xlim(0,15.24) +
    ylim(0,38.1) +
    scale_color_viridis_d() +
    labs(x = "none", y = "none") +
    theme(plot.background = element_blank(),
          # panel.grid.major = element_line(colour="black", size = 0.5),
          # panel.grid.minor = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 3),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
)

# add tracks
p1 <- p +
  # geom_segment(data = df3, aes(x = zone_xj, y = zone_yj, xend = next_zone_xj, yend = next_zone_yj, colour = name), size = 1) +
  # geom_point(data=df3, aes(x=zone_xj, y=zone_yj, colour = name, shape = name), size = 2) +
  # geom_path(data=df3, aes(x=zone_xj, y=zone_yj, colour = name, group = name)) +
  # labs(title = 'Time: {frame_along}') +
  transition_reveal(along = time_day) +
  labs(title = 'Time: {frame_along}') +
  # transition_states(field_time)
  # transition_time(field_time) +
  # transition_components(TimeNum, exit_length = 40) +
  # ease_aes(x = 'sine-out', y = 'sine-out') +
  # shadow_trail(distance = 0.01, size = 0.3)
  # shadow_mark(past = T, future = F, alpha = 0.3)
  # shadow_wake(0.05, size = 2, alpha = TRUE, wrap = FALSE, #exclude_layer = c(2, 3),
  # falloff = 'sine-in', exclude_phase = 'enter')
  # exit_fade()
  # theme(legend.title = element_blank(),
  #       legend.position = "bottom",
#       plot.title = element_text(size = 20))

animate(p1, fps = 10)





# individual plots --------------------------------------------------------
## loop for individuals
ids <- unique(df$name)
ids
i=ids[8]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)

  
  # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  df1 <- subset(df, name == current_mouse)
  ## take only the transitions between novel zones. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$ant_rleid)]
  df3 <- df2 %>% 
    complete(noon_to_noon_day = 1:20)
  
  move_df <- df3 %>% 
    mutate(time = 1:nrow(df3)) %>% 
    mutate(zone_x_j = jitter(df3$zone_x), 
           zone_y_j = jitter(df3$zone_y))
  
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(move_df$treatment), unique(move_df$strain), unique(move_df$sex), unique(move_df$name), sep = "_")
  current_mouse_info <- current_mouse_info[current_mouse_info != "NA_NA_NA_NA"]
  print(current_mouse_info)
  
  ## create backdrop
  (p <- ggplot() +
      geom_rect(data=box, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=antenna_type), color="black", alpha=0.5) +
      # geom_text(data=box, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=antenna), size=4, vjust = -2.2) + 
      xlim(0,15.24) +
      ylim(0,38.1) +
      scale_color_viridis_d() +
      labs(x = "none", y = "none") +
      theme(plot.background = element_blank(),
            # panel.grid.major = element_line(colour="black", size = 0.5),
            # panel.grid.minor = element_line(colour = "black"),
            panel.border = element_rect(colour = "black", fill = NA, size = 3),
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.title = element_blank(),
            legend.text = element_blank(),
            legend.position = "none",
            # legend.text = element_text(size = 15),
            plot.title = element_text(size = 20))
  )
  
  #add animation and titles. ## improvements. pause on days. 
  p2a <- p +
    geom_point(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=seq_along(as.factor(time))), size = 4) + 
    geom_line(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=1), color="red", size=1, alpha = 0.3) +
    labs(title = paste(current_mouse_info, "Day: {move_df$noon_to_noon_day[which.min(abs(move_df$time-frame_along))]}"), 
         subtitle = paste("Total Zone Transitions: {frame_along}")) +
    transition_reveal(time)
  # 
  # ## testing. Make transitions between days even in time.... not quite sure yet how to do this. 
  # (p2a <- p +
  #   geom_point(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=seq_along(as.factor(time))), size = 4) + 
  #   geom_line(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=1), color="red", size=1, alpha = 0.3) +
  #   labs(title = paste(current_mouse_info, "Day: {move_df$noon_to_noon_day[which.min(abs(move_df$time-frame_along))]}"), 
  #        subtitle = paste("Total Zone Transitions: {frame_along}")) +
  #   transition_states(noon_to_noon_day, state_length = 1, transition_length = 1)
  # )
  
  # SAVE GIF
  animate(p2a, width = 600, height = 600, renderer = gifski_renderer())
  anim_save(paste(current_mouse_info, "activity.gif", sep ="_"), path = output_fp)
}
