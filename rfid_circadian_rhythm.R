## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

wd <- getwd()
output_fp <- paste(getwd(), "output", sep = "/")
df <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

df$hours <- as.numeric(format(df$field_time, format="%H"))

(p <- ggplot(df, aes(x=hours, fill = strain_sex)) +
    geom_histogram(binwidth = 1) +
    xlab("Hours") +
    ylab("RFID Reads") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 10, face = "bold"),
          axis.title.x = element_text(color = "black", size = 10, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 10, face = "bold"),
          axis.title.y = element_text(color = "black", size = 10, face = "bold"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.position = "top")
)
ggsave(p, filename = "output/rfid_data_circadian_rhythm.png", device = "png", bg = "white")
ggsave(p, filename = "output/rfid_data_circadian_rhythm.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

