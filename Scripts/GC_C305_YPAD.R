##########################
# Analysis and figure generation for overexpression of C3_05110W in all four genetic backgrounds
###################

# Preparing environment
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
library(data.table)
library(ggplot2)
library(growthcurver)
library(dplyr)
library(cowplot)
library(openxlsx)
#Make sure there is nothing in the environment
rm(list = ls())

#Read in data file
data.dir <- "Data/"
fig.dir <- "Figures/"
Plate1 <- read.table(paste0(data.dir,"2026_04_17_GC_CRISPRA_C305_YPAD.txt"), sep = "\t", header = TRUE)

#Create a time column - measurements were taken every 15min
x <- 0.25
for (i in 1:nrow(Plate1)) {
  Plate1[i,"Time"] <- x
  x <- x + 0.25
}
#Limit this Plate1 to 48hrs
Plate1 <- Plate1[Plate1$Time <= 48,]
#Take a quick look at the raw data for all wells
matplot(Plate1[,1:96],type = "l")

#Read in the Plate1 setup information
Plate1key <- read.xlsx(paste0(data.dir,"2026_04_17_Platekey.xlsx")) #NEed to get actual correct key off of the google drive!

#Subtract the blanks from all wells
Blank <- mean(as.matrix(Plate1[,Plate1key[Plate1key$Strain == "blank","Well"]]))
Plate1[,1:96] <- Plate1[,1:96] - Blank

# Plotting
Plate1melt <- data.table(Plate1[,1:96])
Plate1melt <- melt(Plate1melt)
colnames(Plate1melt)[1] <- "Well"
Plate1melt <- left_join(Plate1melt, Plate1key, by = "Well")
#Need to add time back in 
Plate1melt$Time <- rep(Plate1$Time, 96)

Platemelt <- Plate1melt

# L26 background
ggplot(data = Platemelt[Platemelt$Strain %in% c(grep("2868", unique(Platemelt$Strain), value = TRUE)),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  stat_summary(
    fun = mean,
    geom='line',
    aes(color=Strain)) +
  stat_summary(
    fun=mean,
    geom='point',
    size = 0.8) +
  stat_summary(
    geom='errorbar',
    width=0.01) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("black","red3","red3")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

# P75016 background
ggplot(data = Platemelt[Platemelt$Strain %in% c(grep("2875", unique(Platemelt$Strain), value = TRUE)),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  stat_summary(
    fun = mean,
    geom='line',
    aes(color=Strain)) +
  stat_summary(
    fun=mean,
    geom='point',
    size = 0.8) +
  stat_summary(
    geom='errorbar',
    width=0.01) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("black","red3","red3")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

# AMS5192 background
ggplot(data = Platemelt[Platemelt$Strain %in% c(grep("5192", unique(Platemelt$Strain), value = TRUE)),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  stat_summary(
    fun = mean,
    geom='line',
    aes(color=Strain)) +
  stat_summary(
    fun=mean,
    geom='point',
    size = 0.8) +
  stat_summary(
    geom='errorbar',
    width=0.01) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black")) +
  scale_color_manual(values = c("firebrick", "firebrick","black")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))
ggsave("AMS5192C305YPADGC.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)
