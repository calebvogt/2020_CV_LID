## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)

## As of 10.27.2022, these graphs are made 
# from Mike's combinedSurface.xlsx; plots made in GraphPad using proportion data 
# in the future, I would like to make this via code. He makes his figures using the GBI_summary data sheet

## import gbi_summary
gbi <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_summary.csv")
df<-gbi %>% 
  filter(!(name == "George"), #T004: George only mouse to cross between trials on Day 3. triage.
         !(name == "Anubis" & day >= 5), #T003: Anubis visually confirmed dead by seizure on day 5.
         !(name == "Rae" & day >= 2), #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded.
         !(name == "Hare" & day >= 2), #T004: Hare only appears day 1. Not recovered, presumed dead.
         !(name == "Isis" & day >= 3), #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep.
         !(name == "Rose" & day >= 10)) #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data.
  
options(scipen=4)

# social surface summary --------------------------------------------------

## range of percent mouse time spent alone across all trials. 
df1 <- aggregate(duration_s~m_sum+f_sum+strain+trial,data=df, FUN=sum) ## get sum across all trials for each type of group event. 
df2 <- df1 %>% 
  group_by(trial) %>% 
  mutate(trial_total_h = sum(duration_s/3600)) %>% 
  filter(m_sum==1 & f_sum ==0 | m_sum==0 & f_sum ==1) %>% ## pull out alone durations
  ungroup() %>% 
  group_by(trial) %>% 
  mutate(trial_alone_total_h = sum(duration_s/3600),
         trial_perc_alone=trial_alone_total_h/trial_total_h*100) 
min(df2$trial_perc_alone)
max(df2$trial_perc_alone)




## total sum of social spatiotemporal associations
df1 <- aggregate(duration_s~m_sum+f_sum+strain+trial,data=df, FUN=sum) ## get sum across all trials for each type of group event. 
df2 <- df1 %>% 
  filter(!(m_sum==1 & f_sum ==0 | m_sum==0 & f_sum ==1)) %>% ## remove alone time
  group_by(trial) %>% 
  summarize(trial_total_h = sum(duration_s/3600))
sum(df2$trial_total_h)

## percent mouse social time spent in two mouse spatiotemporal association bouts for C57 and outbred mice
df1 <- aggregate(duration_s~m_sum+f_sum+strain+trial,data=df, FUN=sum) ## get sum across all trials for each type of group event. 
df2 <- df1 %>% 
  filter(!(m_sum==1 & f_sum ==0 | m_sum==0 & f_sum ==1)) %>% ## remove alone time
  group_by(strain) %>% 
  mutate(strain_total_h = sum(duration_s/3600)) %>% 
  filter(m_sum==1 & f_sum ==1 | m_sum==2 & f_sum ==0 | m_sum==0 & f_sum ==2) %>% 
  group_by(strain) %>% 
  summarize(strain_dyadic_h = sum(duration_s/3600),
            perc_social = strain_dyadic_h/strain_total_h*100) %>% 
  unique()


## report the range per trials. 
df2 <- df1 %>% 
  filter(!(m_sum==1 & f_sum ==0 | m_sum==0 & f_sum ==1)) %>% ## remove alone time
  group_by(trial, strain) %>% 
  mutate(strain_total_h = sum(duration_s/3600)) %>% 
  filter(m_sum==1 & f_sum ==1 | m_sum==2 & f_sum ==0 | m_sum==0 & f_sum ==2) %>% 
  group_by(trial,strain) %>% 
  summarize(strain_dyadic_h = sum(duration_s/3600),
            perc_social = strain_dyadic_h/strain_total_h*100) %>% 
  unique()


  

# create social surface table for export to graphpad ----------------------------------------------
df1 <- aggregate(duration_s~m_sum+f_sum+strain,data=df, FUN=sum) ## get sum across all trials for each type of group event. 
df1$duration_h_avg <- ifelse(df1$strain=="C57", df1$duration_s/3600/4/(df1$m_sum+df1$f_sum), df1$duration_s/3600/3/(df1$m_sum+df1$f_sum)) ## take the trial average for each. 

df2<-df1 %>% 
  filter(strain=="C57") ## change for each genotype
  # filter(strain=="NYOB")

surface <- data.frame(matrix(0,nrow=11,ncol=11))
row=1
col=1
for(row in 1:11) {
  for(col in 1:11) {
    print(c(row-1,col-1))
    print(c(row,col))
    try(surface[row,col] <- df2 %>% 
      filter(m_sum==row-1&f_sum==col-1) %>% 
        select(duration_h_avg)
      
    )
     # subset(df2, m_sum==row-1&f_sum==col-1, select="duration_h_avg")         #ifelse(df2$m_sum==row-1 & df2$f_sum==bb, df2$duration_h_avg, 0)
   }
}
surface <- as.matrix(surface)
write.table(df1, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  


# deprecated --------------------------------------------------------------


### Use code to make the heatmap... 
## mikes attempt, using my old code from 2018.... 
t001.surface=read.csv("/users/michaelsheehan/Downloads/t001.secondssurface.csv",header=T)
row.names(t001.surface)=t001.surface$Row
t001.surface$Row=NULL

library(ggplot2)
t001.surface=as.matrix(t001.surface)
ggplot(data = t001.surface, aes(x=Males, y=Females, z=Seconds )) +
  geom_contour(), xlim=(0,10)+ylim(0,)

heatmap(t001.surface)

heatmap(surface, Colv = NA, Rowv = NA, scale="none")

## PLOT C: CREATE POINT_HM FOR PLOTTING TRIAL ATTACK HEATMAP.
point_hm <- point_df
point_hm <- point_hm[,c(4:9)]
colnames(point_hm) <- c("M1", "M2", "M3", "M4", "M5", "M6")
point_hm <- as.matrix(point_hm)
#point_hm <- t(point_hm)
png('T003_Male_Aggression_Heatmap.png')
plot(surface, 
     col = c("beige", "yellow", "orange", "red"),    
     breaks=c(0,75,150,225,300))
     # main = "T003: Heatmap of Male Aggression",
     # xlab= "Male Receiving Aggression", 
     # ylab = "Male Aggressor")
dev.off()