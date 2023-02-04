## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)

dd <- paste(getwd(), "data", sep = "/")
filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)
social_data = lapply(filenames, fread) ## READ IN ALL FILES

# clean social data for triaged mice from social interaction bouts. 
# note that this cleaning step merely deletes columns and adds 0s where appropriate. Does not adjust GBI summary information (m_sum, f_sum, mf_sum)
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  
  df2 <- df %>% 
    #T004: George only mouse to cross between trials on Day 3. triage completely (drop column)
    dplyr::select(-one_of(c( "V1", "NYOB-M-George"))) %>% 
    #T003: Anubis visually confirmed dead by seizure on day 5.  
    mutate_at(vars(one_of(c("C57-M-Anubis"))), ~ ifelse(day >= 5, 0, .)) %>%
    #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
    mutate_at(vars(one_of(c("C57-F-Rae"))), ~ ifelse(day >= 2, 0, .)) %>%
    #T004: Hare only appears day 1. Not recovered, presumed dead. 
    mutate_at(vars(one_of(c("NYOB-M-Hare"))), ~ ifelse(day >= 2, 0, .)) %>%
    #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
    mutate_at(vars(one_of(c("NYOB-F-Isis"))), ~ ifelse(day >= 3, 0, .)) %>%
    #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
    mutate_at(vars(one_of(c("C57-F-Rose"))), ~ ifelse(day >= 10, 0, .))
  
  social_data[[aa]] <- df2
}

# gbi_cumulative_novel_mice.png -------------------------------------------
## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
males <- metadata %>% 
  filter(sex == "M") %>% 
  dplyr::select(name) %>% 
  filter(!is.na(name))
male_list <- dplyr::pull(males, name)

females <- metadata %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  dplyr::select(name) %>% 
  filter(!is.na(name))
female_list <- dplyr::pull(females, name)

trial_stats <- list()
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  colnames(df)<-gsub("C57-M-","",colnames(df))
  colnames(df)<-gsub("C57-F-","",colnames(df))
  colnames(df)<-gsub("NYOB-M-","",colnames(df))
  colnames(df)<-gsub("NYOB-F-","",colnames(df))
  
  col_ids <- colnames(df[,10:ncol(df)]) ## get mouse column names starting at col 10
  bb = col_ids[1]
  all_mouse_list <- list()
  first_flag = 1
  for(bb in col_ids[1:length(col_ids)]) {
    df2 <- df %>% 
      filter((!!as.symbol(bb)) == 1) %>% # pull all rows where bb mouse is present. 
      mutate(name = bb) %>% 
      relocate(name)
    
    non_self_ids <- col_ids[!col_ids %in% bb] # remove current mouse from the next loop to compare to other animals
    i = non_self_ids[1]
    novel_mouse_rows <- list()
    second_flag = 1
    for(i in non_self_ids[1:length(non_self_ids)]) {
      df3 <- df2 %>% 
        filter((!!as.symbol(i)) == 1) %>% 
        mutate(novel_mouse_met = i)
      
      novel_mouse_rows[[second_flag]] <- df3[1,] #save first observed meeting of the focal and novel mouse to list
      second_flag = second_flag + 1
    }
    all_mouse_list[[first_flag]] <- do.call("rbind", novel_mouse_rows)
    first_flag <- first_flag + 1
  }
  df4 <- do.call("rbind", all_mouse_list)
  
  #remove na rows which are introduced when a mouse does not ever meet a particular other mouse. 
  df4[rowSums(is.na(df4)) > 0,]
  df5 <- df4[complete.cases(df4), ] 
  
  ## ADD RELEVANT METADATA INFORMATION. 
  df6 <- merge(df5, metadata, by.x = "name", by.y = "name")
  df7 <- df6 %>%
    dplyr::rename(trial = trial.x) %>% 
    dplyr::select(trial, paddock, strain, sex, name, code, novel_mouse_met, day, zone, field_time_start, field_time_stop, m_sum, f_sum, mf_sum, duration_s)
  
  trial_stats[[aa]] <- df7
}
df8 <- do.call("rbind", trial_stats)
df9 <- df8[with(df8, order(name, day)),]

# df9 has a dataframe ordered by first meeting time with each mouse in the paddock .
df10 <- df9 %>% 
  mutate(strain_sex = paste0(strain, "-", sex)) %>% 
  group_by(trial, strain_sex,strain, sex,name, day) %>% 
  tally() %>% # 
  complete(name, day = full_seq(1:10, period = 1)) %>% #fill in missing day rows for each mouse, adds NAs. 
  replace(is.na(.), 0) %>% 
  mutate(csum = cumsum(n)) %>% #get cumulative # of novel mice met
  arrange(name, day) %>% 
  fill(csum) %>% ## fill cumulative sum data from last observed day
  dplyr::rename(novel_indivs_met = n, csum_novel_indivs_met = csum)

#data cleaning
df11 <- df10 %>% 
  filter(!(name == "George")) %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & day >= 10))  

#plot
df12 <- df11 %>% 
  group_by(strain_sex, day) %>%
  summarise(mean = mean(csum_novel_indivs_met), sd = sd(csum_novel_indivs_met), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df12, aes(x=day, y=mean, color = strain_sex)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("Cumulative novel mice met") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "none") 
)
# ggsave(p, filename = "output/gbi_cumulative_novel_mice.png", device = "png", bg = "white")
ggsave(p, filename = "output/gbi_cumulative_novel_mice.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS
df <- df11
df$trial <- as.factor(df$trial)
df$strain_sex <- as.factor(df$strain_sex)
df$strain <- as.factor(df$strain)
df$sex <- as.factor(df$sex)
df$day <- as.numeric(df$day)

m1 = lmer(csum_novel_indivs_met ~ strain*sex + (1|trial), data = subset(df, day == 10)) 
summary(m1)

m2 = lmer(csum_novel_indivs_met ~ strain*sex + (1|trial)+(1|name), data = subset(df, day == 10)) 
summary(m2)

AIC(m1,m2)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

