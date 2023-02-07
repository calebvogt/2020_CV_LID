## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# for each mouse, get groups of consecutive rfid reads in the same zone. 
# Get the interread intervals for each group, then combine all of these inter-read intervals. 
# Note that this will group consecutive rfid reads seperated by many days. 
# Then get the 99% condfidence interval across every observed interread interval for every read ever, across all mice.
## improved the older method in that it counted interread intervals purely on a zone basis, in that it didnt end bouts
## when the mice were known to have visited another zone. Thus, would have overinflated the interread interval distribution.


library(tidyverse)
library(data.table)
library(readxl)

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
gc()

df <- rfid_data
ids <- unique(df$name)
id_list <- list()
flag <- 1
aa = ids[14]
for(aa in ids[1:length(ids)]){
  print(paste("Processing ", aa, " mouse ",flag, " out of ", length(ids), sep=''))
  df1 <- df %>% 
    filter(name == aa) %>% 
    arrange(field_time) %>%
    mutate(visit_group = NA)
  df1$visit_group[1] <- 1
  
  if(nrow(df1)>1){
    # label consecutive rfid read groups (resource zone visit groups)
    bb = 2
    for(bb in 2:nrow(df1)) {
      if(df1$zone[bb] == df1$zone[bb-1]){
        df1$visit_group[bb] <- df1$visit_group[bb-1] 
      } else{
        df1$visit_group[bb] <- df1$visit_group[bb-1] + 1 
      }
    }
    
    # find interread interval for the visit group and put into list. 
    visit_group <- unique(df1$visit_group)
    visit_group_list <- list()
    cc = 3
    for(cc in visit_group[1:length(visit_group)]){
      df2 <- df1 %>% 
        filter(visit_group == cc) %>% 
        mutate(diff = field_time - lag(field_time), diff_secs = as.numeric(diff, units = 'secs')) %>% 
        filter(diff_secs > 0)
      
      if(nrow(df2) > 0){
        visit_group_list[[cc]] <- df2
      }
    }
  } else{}
  
  id_list[[aa]] <- do.call("rbind", visit_group_list)
  flag <- flag+1
}
df3 <- do.call("rbind", id_list)
df3 <- df3[!is.na(df3$diff_secs),]
summary(df3$diff_secs)
sort(df3$diff_secs)[0.95*length(df3$diff_secs)]  # 18s <- confirmed twice. 
sort(df3$diff_secs)[0.99*length(df3$diff_secs)]  #139s <- most conservative, use this. confirmed twice.  

## Graph
# png(filename = "output/rfid_data_interread_interval_threshold.png")
sdat <- summary(df3$diff_secs)
summStr <- paste(names(sdat), format(sdat, digits = 3), collapse = "; ")
op <- par(mar = c(7,4,4,2) + 0.1)
hist(df3$diff_secs, 
     xlim = c(0, 300),
     # log = "y",
     breaks = 30000,
     main = "",
     xlab = "Within zone inter-read interval (s)"
)
abline(v=c(139), col=c("blue"), lty=c(1,1,1), lwd=c(3,3,3))
title(sub = summStr, line = 5.5)
dev.off()


