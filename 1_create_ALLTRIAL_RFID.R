## create_ALLTRIAL_RFID_DATA
## Caleb C. Vogt, Cornell University

library(tidyverse)
library(readxl)
library(data.table)

setwd(paste(getwd(), "data", sep = "/"))
wd <- getwd()
metadata <- read_excel("LID_2020_metadata.xlsx", sheet = 1, skip = 1)
metadata$tag_1 <- format(metadata$tag_1, digits = 15)
metadata$tag_2 <- format(metadata$tag_2, digits = 15)
weather <- read_excel("LID_2020_metadata.xlsx", sheet = 2, skip = 0)
weather$weather_time <- as.POSIXct(weather$Weather_Time, format = "%m/%d/%Y %H:%M")
folders <- list.dirs(wd, recursive=FALSE, full.names = TRUE)
folders_short <- list.dirs(wd, recursive=FALSE, full.names = FALSE)

i=folders[1]
flag=1
all_trial_list <- list()
for (i in folders[1:length(folders)]) {
  setwd(i)
  trial_var <- folders_short[flag]
  print(paste("Processing", trial_var))
  txt <- list.files(pattern = "*.txt")
  df <- do.call(bind_rows, lapply(txt, function(x) read_table(file = x, col_names = TRUE, col_types = NULL,
                                                                 locale = default_locale(), na = "NA", skip = 4, n_max = Inf,
                                                                 progress = show_progress(), comment = "")))
  df <- df[,1:8] 
  colnames(df) <- c("scan.date", "scan.time", "download.date","download.time", "reader.id", "antenna.id","hex.tag.id" , "dec.tag.id")
  df <- df[!is.na(df$dec.tag.id), ]
  df <- unique(df) #remove repeated read downloads
  df$read_tag <- df$dec.tag.id
  df1 <- merge(df, metadata, by.x ="dec.tag.id", by.y = "tag_1")
  df2 <- merge(df, metadata, by.x = "dec.tag.id", by.y = "tag_2")
  df3 <- bind_rows(df1, df2)
  df3 <- df3[ , c("trial", "paddock", "strain", "sex", "ID", "name", "code", "cage", "zone_drop", "field_age",
                  "family_group", "pre_mass", "post_mass", "body_mm", "tail", "testes",
                  "reader.id", "antenna.id", "scan.date", "scan.time", "read_tag")] 
  df3$field_time <- as.POSIXct(paste(df3$scan.date, df3$scan.time), format="%m/%d/%Y %H:%M:%OS")
  df4 <- df3[order(df3$field_time), ]
  df4$antenna.id <- as.numeric(df4$antenna.id)
  df4$ID <- as.numeric(df4$ID)
  df4$reader.id <- as.numeric(df4$reader.id)
  df4$zone <- df4$antenna.id
  df4$full_ids <- paste0(df4$strain,"-", df4$sex,"-", df4$ID)
  
  # note that if more than 10 antennas are present, grepl will take the first digit and change it
  # CREATE X AND Y COORDINATES FOR THE ZONES FOR PLOTTING ON A GRID
  df4$zone_x <- ifelse(grepl("1", df4$antenna.id), "3.75", 
                   ifelse(grepl("2", df4$antenna.id), "11.25",
                          ifelse(grepl("3", df4$antenna.id), "3.75",
                                 ifelse(grepl("4", df4$antenna.id), "11.25",
                                        ifelse(grepl("5", df4$antenna.id), "3.75",
                                               ifelse(grepl("6", df4$antenna.id), "11.25",
                                                      ifelse(grepl("7", df4$antenna.id), "3.75",
                                                             ifelse(grepl("8", df4$antenna.id), "11.25",
                                                                    "none"))))))))
  
  
  df4$zone_y <- ifelse(grepl("1", df4$antenna.id), "7.6", 
                   ifelse(grepl("2", df4$antenna.id), "7.6",
                          ifelse(grepl("3", df4$antenna.id), "15.2",
                                 ifelse(grepl("4", df4$antenna.id), "15.2",
                                        ifelse(grepl("5", df4$antenna.id), "22.8",
                                               ifelse(grepl("6", df4$antenna.id), "22.8",
                                                      ifelse(grepl("7", df4$antenna.id), "30.4",
                                                             ifelse(grepl("8", df4$antenna.id), "30.4",
                                                                    "none"))))))))
  
  
  # REMOVE OBSERVATIONS FROM TRIALS WHICH ARE NOT POSSIBLE GIVEN PAIRED TRIAL STRUCTURE. 
  ifelse(grepl(pattern = "T001", i), 
         df4 <- subset(df4, trial %in% c('T001')),
         ifelse(grepl(pattern = "T002", i), 
                df4 <- subset(df4, trial %in% c('T002','T003')),
                ifelse(grepl(pattern = "T003", i), 
                       df4 <- subset(df4, trial %in% c('T002','T003')),
                       ifelse(grepl(pattern = "T004", i), 
                              df4 <- subset(df4, trial %in% c('T004','T005')),
                              ifelse(grepl(pattern = "T005", i), 
                                     df4 <- subset(df4, trial %in% c('T004','T005')),
                                     ifelse(grepl(pattern = "T006", i), 
                                            df4 <- subset(df4, trial %in% c('T006','T007')),
                                            ifelse(grepl(pattern = "T007", i), 
                                                   df4 <- subset(df4, trial %in% c('T006','T007')),NA)))))))
  
  # change trial to the ASSIGNED TRIAL BASED ON THE METADATA!!! Which is what the trial column currently reflects. 
  names(df4)[names(df4) == 'trial'] <- "assigned_trial"
  
    # create a new column called trial which reflects where (which paddock) and when (temporally) the data were derived from. 
  # thus, if an animal from T004 somehow crossed into T005 (which ran simultaneously), assigned trial for that obs would read T004, while Trial would 
  # equal T005, which is the trial the mouse was detected in. 
  df4$trial <- trial_var
  origin <- as.POSIXct(paste(df4$scan.date[1], "12:00:00", sep =" "), format="%m/%d/%Y %H:%M:%OS")
  df4$noon_day <- ceiling(difftime(df4$field_time, origin,  units = "days"))
  df4$time_sec <- as.numeric(difftime(df4$field_time,min(df4$field_time),units="secs")+1)
  df4$weather_time <- round(as.POSIXct(df4$field_time), "hour")
  df5 <- left_join(df4,weather)
  
  all_trial_list[[i]] <- df5
  write.csv(df5, file = paste0(folders_short[flag],'_RFID_FULL_DATA.csv'))
  flag <- flag+1
  print(paste("Finished", trial_var))
}
all_trial_data = do.call(bind_rows, all_trial_list)

setwd(wd)
write.csv(all_trial_data, "ALLTRIAL_RFID_DATA_full.csv")


## create cleaner, less data heavy version. 
## Clean data to 10 days and triage relevant subjects as identified by descriptive analyses.
rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA_full.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
rfid_data <- rfid_data %>% 
  mutate(strain_sex = paste(strain, sex, sep = "-")) %>% 
  select(assigned_trial, trial, strain_sex, sex, strain, name, code, cage, family_group, 
         read_tag, scan.date, scan.time, field_time, noon_day,time_sec, zone, zone_x,zone_y) %>% 
  filter(noon_day %in% 1:10) %>%
  filter(!(read_tag == "982.126057708741" & scan.date >= "07/24/2020")) %>%  #T005: Jeeves tag falls out of body on top of antenna)
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage. 
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
gc()

write.csv(rfid_data, "data/ALLTRIAL_RFID_DATA.csv")