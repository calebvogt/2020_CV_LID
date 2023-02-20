## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(magrittr)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

library(MASS)
library(dplyr)
# library(rattle)
library(heplots) # for BoxM test
# library(caret) # can be used for discriminant analysis
library(klaR)



# Figure 2X and Figure S2X: PCA of space use and movement patterns  --------------------------------------------------------------
## To generate PCA variable data from raw data, see Section 2
## Section 1: Create PCA graphs
df <- read_excel("EXTRA_Data_Stats_v2.xlsx", sheet = "FigX_PCA_data", col_names = TRUE) # note, only includes mice that made it to day 10. 
df1 <- df %>% 
  column_to_rownames("name")

#extra schtuff
df2 <- merge(df, meta_short, by.x = "name", by.y = "name")
cols <-colnames(df)

df1 <- df2 %>% 
  filter(strain == "NYOB") %>%
  # filter(strain == "C57") %>% 
  select(all_of(cols)) %>% 
  column_to_rownames("name")



pca <- prcomp(df1, scale = TRUE) # always scale the data. 
plot(pca$x[,1], pca$x[,2])
biplot(pca)
scores <- pca$x
loadings <- pca$rotation
sqrt(1/ncol(df1)) # cutoff for significant loadings
write.table(loadings, "clipboard", sep="\t", row.names=FALSE)
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per[1:3]
sum(pca.var.per[1:3])

# PC scree plot
barplot(pca.var.per, 
        xlab = "PC", 
        ylab = "Percent Variation", 
        names.arg = 1:length(pca.var.per), las = 1)
abline(h=1/ncol(df)*100, col = "red")

# Quick Biplot
biplot(scores[, 1:2], 
       loadings[, 1:2], 
       col = c("red", "blue"),
       cex=c(0.8, 0.8),
       xlabs = rep("o", nrow(scores)),
       xlim = c(-10, 10))


## ggplot
loading_scores <- pca$rotation[,1]
move_scores <- abs(loading_scores)
move_scores_ranked <- sort(move_scores, decreasing = TRUE)
top_moves <- names(move_scores_ranked)
pca$rotation[top_moves,1]
df2 <- data.frame(name = rownames(pca$x), X=pca$x[,1], Y=pca$x[,2])
df3 <- left_join(df2, meta) 
df4 <- df3 %>% 
  mutate(strain_sex = paste(strain, sex, sep= "-")) %>% 
  select(trial, strain_sex, strain,sex, name, X, Y)

(p <- ggplot(df4, aes(x = X, y = Y, color = strain_sex, shape = strain_sex)) +
    geom_point(aes(shape = factor(strain_sex))) +
    # scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M")) + 
    # values=c("red", "blue", "green", "darkgrey")) +
    theme_classic() +
    xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
    theme(axis.text.x = element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.key.size = unit(0.2, 'cm'), 
          legend.text = element_text(size=7),
          legend.position = "top"))
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent")

## PCA Stats
# 1. Take PCA scores for each mouse and associate with other (per mouse) data
# 2. Then do something like the follow. Note that Comp.1 here is equivalent to PC1. 
tracking60s.hometub=read.csv("/Users/michaelsheehan/Downloads/tracking60s.hometub.csv", header=T)

pca.mod1=lmer(Comp.1~night+(1|mouse), data=tracking60s.hometub)
anova(pca.mod1)

pca.mod2=lmer(Comp.2~night+(1|mouse), data=tracking60s.hometub)
anova(pca.mod2)



# Section 2: Get all data that you want to put into the PCA. 
# total_zone_time_h
df <- move_data %>% 
  mutate(duration_hours = duration_s / 60 / 60) %>%
  group_by(name) %>% 
  tally(sum(duration_hours)) 
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# mean_bout_visit_s
df <- move_data %>% 
  group_by(name) %>% 
  tally(mean(duration_s)) 
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_unique_zones_per_trial
df <- move_data %>% 
  group_by(name,antenna) %>% 
  tally() %>% 
  dplyr::count(name)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# mean_dist_zones_per_night
df <- move_data %>% # 
  group_by(name, noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each zone
  group_by(name, noon_day) %>%
  tally() %>% # get number of distinct zone visits
  select(name, n) %>% 
  summarize_all(funs(mean))
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# %_time_rank_1_zone
options(scipen=999)
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) %>% 
  group_by(name) %>% 
  slice(which.max(percent_time)) %>% 
  select(name, percent_time)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# %_time_rank_2_zone
options(scipen=999)
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  mutate(percent_time = n / sum(n)) %>% 
  arrange(desc(percent_time)) %>% 
  group_by(name) %>% 
  slice(-which.max(percent_time)) %>% 
  slice(which.max(percent_time)) %>% 
  select(name, percent_time)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_time_rank_1_zone
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  arrange(desc(n)) %>% 
  group_by(name) %>% 
  slice(which.max(n)) %>% 
  select(name, n)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

# total_time_rank_2_zone
df <- move_data %>% 
  mutate(duration_min = duration_s / 60) %>% 
  group_by(strain, sex, name, antenna) %>% 
  tally(sum(duration_min)) %>% 
  arrange(desc(n)) %>% 
  group_by(name) %>% 
  slice(-which.max(n)) %>% 
  slice(which.max(n)) %>% 
  select(name, n)
write.table(df, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD  

## Create df to get avg and min daily distance travelled. 
df <- move_data %>% 
  select(trial, strain_sex,sex,strain,name, noon_day, zone, zone_x, zone_y, antenna)

df$zone_x[df$zone_x =="A"] <- 3.75
df$zone_x[df$zone_x =="B"] <- 11.25
df$zone_y[df$zone_y =="A"] <- 7.6
df$zone_y[df$zone_y =="B"] <- 15.2
df$zone_y[df$zone_y =="C"] <- 22.8
df$zone_y[df$zone_y =="D"] <- 30.4

df$zone_x <- as.numeric(df$zone_x)
df$zone_y <- as.numeric(df$zone_y)
summary(df)

ids <- unique(df$name)
data_list <- list()
aa = ids[60]
for(aa in ids[1:length(ids)]){
  ## create df of percent time
  df1 <- df %>% 
    filter(name == aa) 
  # delete consecutive repeat antenna hits, only keep rows where antennas change. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$zone)]
  df2$dist <- NA
  n <- nrow(df2)
  
  if(n==1){
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  } else{
    df2$dist[2:n] <- sqrt((df2$zone_x[2:n] - df2$zone_x[1:n-1]) ^ 2 + (df2$zone_y[2:n] - df2$zone_y[1:n-1]) ^ 2)
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  }
}

df3 <- do.call("rbind", data_list)
summary(df3)
df3$noon_day <- as.numeric(df3$noon_day)

df4 <- df3 %>% 
  group_by(trial, strain_sex,strain, sex, name, noon_day) %>% 
  tally(sum(dist)) %>% 
  complete(noon_day = 1:10, fill = list(n = 0)) %>%  # fill in days where mouse doesnt appear with 0s.
  dplyr::rename(dist = n)

# data cleaning. 
df5 <- df4 %>% 
  #T003: Anubis visually confirmed dead by seizure on day 5.  
  filter(!(name == "Anubis" & noon_day >= 5)) %>% 
  #T003: Rae appears once on the first day, but she is captured at the end of the trial. Only female to do this, so excluded. 
  filter(!(name == "Rae" & noon_day >= 2)) %>%  
  #T004: Hare only appears day 1. Not recovered, presumed dead. 
  filter(!(name == "Hare" & noon_day >= 2)) %>% 
  #T004: Isis lost after day 2. Not recovered, presumed dead. #T004: Gilmour lost on day 10 only but recovered/trapped. Keep. 
  filter(!(name == "Isis" & noon_day >= 3)) %>% 
  #T003: Rose lost on Day 10, but trapped WITHOUT RFID tag. triage day 10 data. 
  filter(!(name == "Rose" & noon_day >= 10))  

# total_min_distance_traveled
df4 <- df3 %>% 
  group_by(strain_sex, name, noon_day) %>% 
  tally(sum(dist_moved)) %>% 
  group_by(name) %>%
  tally(sum(n))
write.table(df4, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD 

# avg_min_dist_traveled_per_night
df4 <- df3 %>% 
  group_by(name, noon_day) %>% 
  tally(mean(dist_moved)) %>% 
  complete(noon_day = 1:10, fill = list(n = 0)) # fill in days where mouse doesnt appear with 0s
df4 <- subset(df4, n > 0 | (name != "Anubis" & name != "Hare" & name != "Gilmore" & name != "Isis" & name != "Ray" & name != "Rose"))## if mouse has all 0s consecutively till the end, delete 0 rows for that animal. 
df5 <- df4 %>% 
  group_by(name) %>%
  tally(mean(n))
write.table(df5, "clipboard", sep="\t", row.names=FALSE) # COPY THE OUTPUT TO THE CLIPBOARD 



# Figure S2C: Quadratic Discriminant Analysis --------------------------------------------


df <- read_excel("Data_Stats_v6.xlsx", sheet = "FigX_QDA_data", col_names = TRUE) # note, only includes mice that made it to day 10. 
df1 <- merge(df, meta_short, by = "name")
df2 <- df1 %>% 
  mutate(strain_sex = paste(strain, sex, sep= "-")) %>% 
  relocate(strain_sex) %>% 
  dplyr::select(-c(name, paddock, sex, strain, family_group, trial, code))
df2$strain_sex <- as.factor(df2$strain_sex)

# LDA model
#check assumptions for LDA
boxM(df2[,2:12], df2$strain_sex) # technically, a significant value suggests that LDA assumptions are not met, and therefore quadratic discriminant analysis should be used... 


model <- lda(strain_sex ~ ., data = df2)
model.p <- predict(model)$class
model.p_table <- table(model.p, df2$strain_sex) # model.p are rows truth are columns
mean(model.p == df2$strain_sex)
write.table(model.p_table, "clipboard-16384", sep="\t", row.names=T, col.names=T) #col.names = false for second time. 

# cross validate the model. leave one out approach. 
model2 <- lda(strain_sex ~ ., data = df2, CV = TRUE)
model2.p <- model2$class
model2.p_table <- table(model2.p, df2$strain_sex)
mean(model2.p == df2$strain_sex)
table(model2$class, df2$strain_sex) # performs worse. 
plot(model, dimen = 1, type = "b")

#partimat plots
library(klaR)
# plot the linear discriminant functions for highest PC1 variables? 
partimat(strain_sex ~ ., data = df2, col.correct = 'green', col.wrong = 'red', method="lda")

## get the LD1 (x) and LD2 (y) coordinates for the LDA plot
data.lda.values <- predict(model)
## create a dataframe that has all the info we need to draw a graph
plot.data <- data.frame(X=data.lda.values$x[,1], Y=data.lda.values$x[,2], strain_sex = df2$strain_sex)
head(plot.data)

## draw a graph using ggplot2
(p <- ggplot(data=plot.data, aes(x=X, y=Y, color = strain_sex, shape = strain_sex)) +
    geom_point(aes(shape = factor(strain_sex))) +
    scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       values=c("red", "blue", "green", "darkgrey")) +
    xlab("LD1") +
    ylab("LD2") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8, face = "bold"),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8, face = "bold"),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.key.size = unit(0.2, 'cm'), 
          legend.text = element_text(size=7),
          legend.position = "top")
) 
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.5, bg = "transparent")


## QDA model
set.seed(123)
train <- sample(c(T,F),nrow(df2), replace = T, prob = c(2/3, 1/3))
trainset <- df2[train,]
test <- df2[!train,]

model <- qda(strain_sex ~ ., data = trainset)
model.p <- predict(model, test)
mean(model.p$class == test$strain_sex) # 75% accuracy. 
table(model.p$class, test$strain_sex, dnn = c('Predicted Group', 'Actual Group'))


## QDA modelv2. C57-F vs OTHER TBD. 

# Partition plots. Creates a shit ton of plots
library(klaR)
partimat(strain_sex ~ ., data = df2, method = "qda", plot.matrix = T, col.correct = 'green', col.wrong = 'red')

## get the x,y coordinates for the LDA plot
data.qda.values <- predict(model)$class
## create a dataframe that has all the info we need to draw a graph
plot.data <- data.frame(X=data.qda.values$x[,1], Y=data.lda.values$x[,2], strain_sex = df2$strain_sex)
head(plot.data)


## QDA model
set.seed(1234)
train <- sample(c(T,F),nrow(df2), replace = T, prob = c(0.6, 0.4))
trainset <- df2[train,]
test <- df2[!train,]

model <- lda(strain_sex ~ ., data = trainset)
model.p <- predict(model, test)
mean(model.p$class == test$strain_sex) # 75% accuracy. 
table(model.p$class, test$strain_sex, dnn = c('Predicted Group', 'Actual Group'))




