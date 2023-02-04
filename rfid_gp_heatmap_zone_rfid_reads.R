## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

rfid_data <- rfid_data %>% 
  filter(noon_day %in% 1:10) %>%
  filter(!(read_tag == "982.126057708741" & scan.date >= "07/24/2020")) %>%  #T005: Jeeves tag falls out of body on top of antenna)
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & noon_day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.  
         !(name == "Rae" & noon_day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
         !(name == "Hare" & noon_day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead. 
         !(name == "Isis" & noon_day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
         !(name == "Rose" & noon_day >= 10)) %>% #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  mutate(strain_sex = paste(strain, sex, sep = "_"))


# GP heat map scale should be set to the max visits for any given category. Suri = 184 visits. 
df <- rfid_data %>%
  # filter(trial == "T001", sex == "M") %>%
  # filter(trial == "T001", sex == "F") %>%
  # filter(trial == "T002", sex == "M") %>%
  # filter(trial == "T002", sex == "F") %>%
  # filter(trial == "T003", sex == "M") %>%
  # filter(trial == "T003", sex == "F") %>%
  # filter(trial == "T004", sex == "M") %>%
  # filter(trial == "T004", sex == "F") %>%
  # filter(trial == "T005", sex == "M") %>%
  # filter(trial == "T005", sex == "F") %>%
  # filter(trial == "T006", sex == "M") %>%
  # filter(trial == "T006", sex == "F") %>%
# filter(trial == "T007", sex == "M") %>%
filter(trial == "T007", sex == "F") %>%
  group_by(name, code, noon_day, zone) %>%
  tally() #number of visits

#check maximum visit number! adjust GP heatmap settings accordingly. 
df$n<-log10(df$n)

codes <- unique(df$code)
ids <- unique(df$name)#CREATE LOOP FOR ADDING NUMBER OF VISITS TO STATS FOR AN ENTIRE TRIAL. 
idlist <- list() # NOTE THAT DATA FRAMED WILL BE IN ORDER OF THIS LIST
daylist<- list()
flag=1
# aa = ids[1]
for(aa in ids[1:length(ids)]){ # LOOP THROUGH EACH INDIVIDUAL AND PULL OUT NUMBER OF VISITS PER UNIQUE zone AND PUT INTO 2X4 GRID THAT LOCALIZES TO THE 
  # CREATE STATS DATAFRAME TO MIMIC LAYOUT OF FIELD SITE
  stats <- data.frame(matrix(0, nrow = 4, ncol = 2)) # ENCLOSURE SETUP. PUT EACH NIGHT OF ACTIVITY TO THE RIGHT FOR 10 NIGHTS. COPY AND PASTE DIRECLY INTO PRISM. 
  df1 <- subset(df, name == aa)
  for(bb in 1:10){
    df2 <- subset(df1, noon_day == bb)
    # CHECK IF THERE WERE ANY DETECTED VISITS THAT NIGHT. IF YES, JUST PUT A ZERO. IF NO, ADD VISIT NUMBERS BY zone
    if(nrow(df2) == 0){
      stats <- data.frame(matrix(0, nrow = 4, ncol = 2))
    } else {
      
      for(cc in 1:nrow(df2)){
        if(df2$zone[cc] == 1){
          stats[4,1] <- print(df2$n[cc])
        } else if(df2$zone[cc] == 2){
          stats[4,2] <-print(df2$n[cc])
        } else if(df2$zone[cc] == 3){
          stats[3,1] <- print(df2$n[cc])
        } else if(df2$zone[cc] == 4){
          stats[3,2] <- print(df2$n[cc])
        } else if(df2$zone[cc] == 5){
          stats[2,1] <- print(df2$n[cc])
        } else if(df2$zone[cc] == 6){
          stats[2,2] <- print(df2$n[cc])
        } else if(df2$zone[cc] == 7){
          stats[1,1] <- print(df2$n[cc])
        } else if (df2$zone[cc] == 8){
          stats[1,2] <- print(df2$n[cc])
        } else {print("barnacles")}
        
      }
    }
    daylist[[bb]] <- stats
    stats <- data.frame(matrix(0, nrow = 4, ncol = 2)) #CHANGE FROM NA 'S TO 0 'S
  }
  master_class <- do.call("cbind",daylist) #THIS THROWS THE ERROR
  idlist[[flag]] <- master_class
  flag=flag+1
}
master_SASS <- do.call("rbind",idlist) # RBIND DATA FOR EACH INDIVIDUAL
# View(master_SASS)
# head(master_SASS)
write.table(master_SASS, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  
ids # PRINT ORDER OF THE DATA SET, COPY INTO GRAPHPAD
codes

