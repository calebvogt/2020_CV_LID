## Created by Lucie Michel 
## Edited/modified by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(dplyr)
library(plotrix)
library(readxl)    

metadata <- metadata <- read_excel("data/LID_2020_metadata.xlsx", sheet = 1, skip = 1)

setwd("data/video_scoring/")

## read in data.
## lucie watched EVERY purported GBI event where there were >1 males. 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, 
                                                     sheet = X, 
                                                     col_names=T,
                                                     skip=1,
                                                     col_types=c("text")))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

list <- read_excel_allsheets("LID_2020_MOVEBOUT_GBI_scoring.xlsx")
list1 <- list[-c(8,9)]
list2env(list1,envir=.GlobalEnv)
list2 <- list()
aa=1
for(aa in 1:length(list1)) {
  df <- list1[[aa]]
  list2[[aa]] <- df %>% 
    filter(!(no_video == "TRUE" | no_mice_visible == "TRUE")) %>%
    # filter(m_sum >1 & agg == "TRUE") %>%
    filter(m_sum >1) %>% 
    rename(trial = Trial, day=Day, zone=Zone, field_time_start=Field_Time_Start, field_time_stop=Field_Time_Stop,checked=checked...15) %>% 
    select(trial, day, zone,field_time_start, field_time_stop, duration_s, m_sum, f_sum,mf_sum,max_mice_obs,no_video, no_mice_visible,checked,agg,agg_type,agg_actor,agg_recipient)
  print(aa)
}
df = do.call(bind_rows, list2)
df$field_time_start <- as.POSIXct(as.numeric(df$field_time_start) * (60*60*24), origin="1899-12-30", tz="GMT")
df$field_time_stop <- as.POSIXct(as.numeric(df$field_time_stop) * (60*60*24), origin="1899-12-30", tz="GMT")

## code here. 
df$strain <- ifelse(df$trial == "T001", "C57", 
                    ifelse(df$trial == "T002", "C57", 
                           ifelse(df$trial == "T003", "C57",
                                  ifelse(df$trial == "T004", "Outbred",
                                         ifelse(df$trial == "T005", "Outbred",
                                                ifelse(df$trial == "T006", "C57",
                                                       ifelse(df$trial == "T007", "Outbred", NA)))))))

# Summary -----------------------------------------------------------------

df1 <- df %>%
  



# Number of visible fighting/chasing events per day  ---------------------------------------
df1 <- df %>%
  mutate(day=as.numeric(day)) %>% 
  group_by(strain, trial, day) %>% 
  tally() %>% 
  rename(agg_events=n) %>% 
  # group_by(trial,strain) %>% 
  complete(day = full_seq(1:10, period = 1)) %>% 
  replace(is.na(.), 0)
  
df2 <- df1 %>%
  group_by(strain, day) %>% 
  summarise(mean = mean(agg_events), 
            sd = sd(agg_events),
            count = n(), 
            se = (sd/(sqrt(count))))

(p <- ggplot(df2, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna",  "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("# fighting/chasing events") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.8,0.8))
          # legend.position = "")
)
# ggsave(p, filename = "output/rfid_data_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_data_zones_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")



# Percent of visible >1 male overlap events with fighting/chasing events per day ---------------------------------------
df1 <- df %>%
  mutate(day=as.numeric(day)) %>% 
  group_by(strain, trial, day) %>% 
  mutate(total_events = n()) %>% 
  group_by(strain, trial, day,total_events,agg) %>% 
  tally() %>% 
  filter(agg=="TRUE") %>% 
  mutate(perc_agg = n/total_events*100) %>% 
  ungroup() %>% 
  group_by(trial,strain) %>%
  complete(day = full_seq(1:10, period = 1)) %>% 
  mutate_all(~replace(., is.na(.),0))

df2 <- df1 %>%
  group_by(strain, day) %>% 
  summarise(mean = mean(perc_agg), 
            sd = sd(perc_agg),
            count = n(), 
            se = (sd/(sqrt(count))))

(p <- ggplot(df2, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna",  "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% of M-M events with aggression") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.8,0.8))
  # legend.position = "")
)
# ggsave(p, filename = "output/rfid_data_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_data_zones_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")




# Percent daily M-M events not visible (no video/no visible mice)  ---------------------------------------
df1 <- df %>%
  mutate(day=as.numeric(day)) %>% 
  group_by(strain, trial, day) %>% 
  mutate(total_events = n()) %>% 
  mutate(visible = ifelse(no_video == "TRUE" | no_mice_visible == "TRUE", 0,1)) %>% 
  group_by(strain, trial, day,total_events,visible) %>% 
  tally() %>% 
  filter(visible==1) %>% 
  mutate(perc_visible = n/total_events*100) 
  # ungroup() %>%
  # group_by(trial,strain) %>%
  # complete(day = full_seq(1:10, period = 1)) %>%
  # mutate_all(~replace(., is.na(.),0))

df2 <- df1 %>%
  group_by(strain, day) %>% 
  summarise(mean = mean(perc_visible), 
            sd = sd(perc_visible),
            count = n(), 
            se = (sd/(sqrt(count))))

(p <- ggplot(df2, aes(x = day, y = mean, color = strain)) +
    geom_line(size = 0.75) +
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    scale_color_manual(breaks = c("C57", "Outbred"),
                       values=c("sienna",  "skyblue4")) +
    theme_classic() +
    xlab("Day") +
    ylab("% of rfid detected M-M events that are visible on camera") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.8,0.8))
  # legend.position = "")
)
# ggsave(p, filename = "output/rfid_data_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_data_zones_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")







### What else? show that displacement of males by territory holders occurs under both agg and non-agg conditions? 

# Lucie's Original Code ---------------------------------------------------
setwd("data/video_scoring/")
vect=c("T001","T002","T003","T004","T005","T006","T007")

#### plotting the number of interactions and agressions in general per day for each strain ####

nbint=matrix(0,nrow=10,ncol=7) #number of interactions per trial per day
days= c("1","2","3","4","5","6","7","8","9","10")
row.names(nbint) = days
colnames(nbint) = vect
nbagg=nbint #number of agressions per trial per day
propagg=nbagg # proportion of agressions within the total of interactions between males

for (i in 1:10) {
  for (j in vect) {
    file=paste(j,".csv", sep="")
    a=read.csv(file, header=TRUE)
    b=subset(a,a$m_sum>1 & a$no_video==FALSE & a$no_mice_visible == FALSE & a$Day==i) #all interactions
    c=subset(b, b$agg == TRUE) # only the agressive ones
    nbint[as.character(i),j]=length(b$agg)
    nbagg[as.character(i),j]=length(c$agg)
    propagg[as.character(i),j]=(length(c$agg)/length(b$agg)*100)
    }
}

resume=matrix(0,nrow=10,ncol=6) #summing up the mean number of interactions per trial per day for both c57 and outbred mice
row.names(resume)=days
colnames(resume)=c("nbintc57","nbaggc57","nbintoutbred","nbaggoutbred","propaggc57","propaggoutbred")
for (i in 1:10) {
  for (j in c("T001","T002","T003","T006")) {
    resume[as.character(i),"nbintc57"]=nbint[as.character(i),j]+resume[as.character(i),"nbintc57"]
    resume[as.character(i),"nbaggc57"]=nbagg[as.character(i),j]+resume[as.character(i),"nbaggc57"]
    resume[as.character(i),"propaggc57"]=propagg[as.character(i),j]+resume[as.character(i),"propaggc57"]
  }
  for (j in c("T004","T005","T007")) {
    resume[as.character(i),"nbintoutbred"]=nbint[as.character(i),j]+resume[as.character(i),"nbintoutbred"]
    resume[as.character(i),"nbaggoutbred"]=nbagg[as.character(i),j]+resume[as.character(i),"nbaggoutbred"]
    if (is.na(propagg[as.character(i),j])){
      resume[as.character(i),"propaggoutbred"]=0+resume[as.character(i),"propaggoutbred"]  #replacing the NA values by 0 to be able to have a continuous graph
    }else{
      resume[as.character(i),"propaggoutbred"]=propagg[as.character(i),j]+resume[as.character(i),"propaggoutbred"]
    }
    }
        resume[as.character(i),"nbintc57"]=resume[as.character(i),"nbintc57"]/4
  resume[as.character(i),"nbaggc57"]=resume[as.character(i),"nbaggc57"]/4
  resume[as.character(i),"propaggc57"]=resume[as.character(i),"propaggc57"]/4
  resume[as.character(i),"nbintoutbred"]=resume[as.character(i),"nbintoutbred"]/3
  resume[as.character(i),"nbaggoutbred"]=resume[as.character(i),"nbaggoutbred"]/3
  resume[as.character(i),"propaggoutbred"]=resume[as.character(i),"propaggoutbred"]/3
}
resumedataframe=as.data.frame(resume) #resume was a matrix since then

par(mfrow=c(1,2)) #plotting everything : the interactions for each breed and the proportion of aggressive interactions for both
plot(resumedataframe$nbintc57, type='l', ylim = range(c(resumedataframe$nbintoutbred - std.error(t(nbint[,c(4,5,7)])), resumedataframe$nbintc57 + std.error(t(nbint[,c(1,2,3,6)])))), 
                                           pch=20,cex=2, col="blue", xlab="Days", ylab="Number of interactions", main="Mean number of interactions per day")
points(resumedataframe$nbintc57, ylim = range(c(resumedataframe$nbintoutbred - std.error(t(nbint[,c(4,5,7)])), resumedataframe$nbintc57 + std.error(t(nbint[,c(1,2,3,6)])))), 
     pch=20,cex=2, col="blue", xlab="Days", ylab="Number of interactions", main="Mean number of interactions per day")
lines(resumedataframe$nbintoutbred, col="darkorange", xlab="Days", ylab="Number of interactions")
points(resumedataframe$nbintoutbred, pch=20,cex=2, col="darkorange", xlab="Days", ylab="Number of interactions")
arrows(x0 = 1:10, y0 = resumedataframe$nbintc57 - std.error(t(nbint[,c(1,2,3,6)])) , x1 = 1:10,
       y1 = resumedataframe$nbintc57 + std.error(t(nbint[,c(1,2,3,6)])),       length = 0.15, code = 3, angle = 90)
arrows(x0 = 1:10, y0 = resumedataframe$nbintoutbred - std.error(t(nbint[,c(4,5,7)])) , x1 = 1:10,
       y1 = resumedataframe$nbintoutbred + std.error(t(nbint[,c(4,5,7)])),       length = 0.15, code = 3, angle = 90)
legend("topright",legend=c('C57',"Outbred"),lwd=3, col=c("blue","darkorange"))


plot(resumedataframe$nbaggc57,type= "l", ylim = range(c(resumedataframe$nbaggoutbred - std.error(t(nbagg[,c(4,5,7)])), resumedataframe$nbaggc57 + std.error(t(nbagg[,c(1,2,3,6)])))),
     pch=20,cex=2, col="blue",xlab="Days", ylab = "Number of agressions",  main="Mean number of aggressions per day")
points(resumedataframe$nbaggc57, ylim = range(c(resumedataframe$nbaggoutbred - std.error(t(nbagg[,c(4,5,7)])), resumedataframe$nbaggc57 + std.error(t(nbagg[,c(1,2,3,6)])))),
     pch=20,cex=2, col="blue",xlab="Days", ylab = "Number of agressions",  main="Mean number of aggressions per day")
lines(resumedataframe$nbaggoutbred, col="darkorange")
points(resumedataframe$nbaggoutbred, pch=20,cex=2, col="darkorange")
arrows(x0 = 1:10, y0 = resumedataframe$nbaggc57 - std.error(t(nbagg[,c(1,2,3,6)])) , x1 = 1:10,
       y1 = resumedataframe$nbaggc57 + std.error(t(nbagg[,c(1,2,3,6)])),       length = 0.15, code = 3, angle = 90)
arrows(x0 = 1:10, y0 = resumedataframe$nbaggoutbred - std.error(t(nbagg[,c(4,5,7)])) , x1 = 1:10,
       y1 = resumedataframe$nbaggoutbred + std.error(t(nbagg[,c(4,5,7)])),       length = 0.15, code = 3, angle = 90)
legend("topright",legend=c('C57',"Outbred"),lwd=3, col=c("blue","darkorange"))

plot(resumedataframe$propaggc57,type='l',cex=1, col="blue", xlab="Days", ylab="% of agressive interactions", main="% of agressive interactions")
lines(resumedataframe$propaggoutbred,pch=20,cex=1, col="darkorange")
legend("topright",legend=c('C57',"Outbred"),lwd=3, col=c("blue","darkorange"))




#### plotting the mean number of aggressive bouts per mouse per day for each strain ####


permouse=data.frame(matrix(ncol = 3, nrow = 1))
columnames=c("agg_actor",'Day',"n")
colnames(permouse)=columnames #creating a data frame with the number of agressive bouts per mouse per day

for (i in vect) {
  file=paste(i,".csv", sep="")
  df = as.data.frame(read.csv(file, header=TRUE))
  trial=as.data.frame(df %>% filter (agg==TRUE & agg_actor!="") %>% group_by(agg_actor, Day)%>%tally())
  permouse=rbind(permouse,trial)
}
permouse=na.omit(permouse) #added a blank row at the beginning bc it wouldn't work otherwise

meanpermouse=data.frame(Day=c(1:10), C57=rep(0,10), OB=rep(0,10)) #creating a data frame with the mean number of agressive bouts per mouse per day
stdepermouse=data.frame(Day=c(1:10), C57=rep(0,10), OB=rep(0,10)) #creating a data frame with the mean number of agressive bouts per mouse per day
for (i in c(1:5,7)) { #no interactions on day 6 and 8 to 10 : it would cause errors in the script
  temp = permouse %>% filter(Day==i)
  stdec57= temp %>% filter(grepl("C57",temp[,1]))
  stdeob= temp %>% filter(grepl("OB",temp[,1]))
  stdepermouse[i,2]=std.error(stdec57[,3])
  stdepermouse[i,3]=std.error(stdeob[,3])
  meanc57=0
  meanob=0
  for (j in 1:nrow(temp)) {
    if (grepl("C57",temp[j,1])){
      meanpermouse[i,2]=meanpermouse[i,2]+temp[j,3]
      meanc57=meanc57+1
    }else{
      meanpermouse[i,3]=meanpermouse[i,3]+temp[j,3]
      meanob=meanob+1
    }
  }
  meanpermouse[i,3]=meanpermouse[i,3]/meanob
  meanpermouse[i,2]=meanpermouse[i,2]/meanc57
  #meanpermouse[i,3]=meanpermouse[i,3]/30
  #meanpermouse[i,2]=meanpermouse[i,2]/40
}
meanpermouse[is.na(meanpermouse)]=0 #on some days there was no interaction for some mice so impossible to divide by 0
stdepermouse[is.na(stdepermouse)]=0

#stdepermouse10=data.frame(Day=c(1:10), C57=rep(0,10), OB=rep(0,10))
#for (i in 1:10) {
 # tempc57 = permouse %>% filter(Day==i & grepl("C57", temp[,1]))
  #tempob= permouse %>% filter(Day==i & grepl("OB", temp[,1]))
  #sumc57=0
  #sumob=0
  #for (j in 1:nrow(tempc57)) {
   #   sumc57=sumc57+(tempc57[j,1]-meanpermouse[i,2])**2+(40-nrow(tempc57))*(meanpermouse[i,2]**2)
      
    #  sumob=sumob+(temp[j,2]-meanpermouse[i,3])**2
  #}
#}


plot(meanpermouse$C57, ylim = range(c(meanpermouse$C57 - stdepermouse$C57, meanpermouse$C57 + stdepermouse$C57)),
     type="l", col="blue", xlab="Days", ylab="Number of agressive bouts per mouse",
     main="Mean number of agressive bouts per mouse per day")
lines(meanpermouse$OB, col="darkorange")
arrows(x0 = meanpermouse$Day, y0 = meanpermouse$C57 - stdepermouse$C57 , x1 = meanpermouse$Day,
       y1 = meanpermouse$C57 + stdepermouse$C57,       length = 0.15, code = 3, angle = 90)
arrows(x0 = meanpermouse$Day, y0 = meanpermouse$OB - stdepermouse$OB , x1 = meanpermouse$Day,
       y1 = meanpermouse$OB + stdepermouse$OB,       length = 0.15, code = 3, angle = 90)
legend("topright",legend=c('C57',"Outbred"),lwd=3, col=c("blue","darkorange"))
