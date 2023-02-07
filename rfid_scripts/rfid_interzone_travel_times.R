## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

df <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

df$zone_copy <- df$zone # create copy of zone, as rleid deletes the focal column
ids <- unique(df$name)
day_list <- list()
data_list <- list()
aa = ids[4]
for(aa in ids[1:length(ids)]){
  df1 <- df %>% 
    filter(name == aa) 
  
  bb = 1
  for(bb in 1:max(df1$noon_day)) { #get daily inter-travel times. 
    df2 <- df1 %>% 
      filter(noon_day == bb)
    
    df3 <- as.data.table(df2)[, .SD[1], by = rleid(df2$zone_copy)] # delete consecutive repeat antenna hits, only keep rows where zone change. 
    day_list[[bb]] <- df3 %>% 
      select(trial, sex, name, zone, field_time) %>% 
      mutate(diff = field_time - lag(field_time), 
             diff_secs = as.numeric(diff, units = 'secs'), 
             diff_mins = as.numeric(diff, units = 'mins'), 
             diff_hours = as.numeric(diff, units = 'hours'))
  }
  data_list[[aa]] <- do.call("rbind", day_list)
}
df4 <- do.call("rbind", data_list)
df5 <- na.omit(df4)
min(df5$diff_secs)
mean(df5$diff_mins)
max(df5$diff_hours)

# Graph
# sdat <- summary(df5$diff_mins)
# summStr <- paste(names(sdat), format(sdat, digits = 0), collapse = "; ")
# op <- par(mar = c(7,4,4,2) + 0.1)
hist(df5$diff_mins, 
     xlim = c(0, 1000),
     breaks = 10000,
     # main = stuff,
     main = "",
     xlab = "Inter-zone travel time (min)"
)
# title(sub = summStr, line = 5.5)
# par(op)
# export as svg. 

