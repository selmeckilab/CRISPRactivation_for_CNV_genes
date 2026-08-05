#########################
# Code for analysis and figure generation for the overexpression of C3_03370C in micafungin
#######

# Prepare environment
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
Plate1 <- read.table(paste0(data.dir,"2026_01_30_GC_CRISPRA_MCF1_C3.txt"), header = 1)
Plate2 <- read.xlsx(paste0(data.dir,"2026_01_30_GC_CRISPRA_MCF2_C3.xlsx"), sheet = 1, rows = c(40:233), cols = c(4:100))
Plate3 <- read.xlsx(paste0(data.dir,"2026_01_30_GC_CRISPRA_MCF3_C3.xlsx"), sheet = 1, rows = c(40:233), cols = c(4:100))

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
Plate1key <- read.xlsx(paste0(data.dir,"2026_01_30_GCPlatekey_MCF.xlsx"), colNames = TRUE) 
#Subtract the blanks from all wells
Blank <- mean(as.matrix(Plate1[,Plate1key[Plate1key$Strain == "blank","Well"]]))
Blank
Plate1[,1:96] <- Plate1[,1:96] - Blank

Plate1melt <- data.table(Plate1[,1:96])
Plate1melt <- melt(Plate1melt)
colnames(Plate1melt)[1] <- "Well"
Plate1melt <- left_join(Plate1melt, Plate1key, by = "Well")
#Need to add time back in 
Plate1melt$Time <- rep(Plate1$Time, 96)

#Create a time column - measurements were taken every 15min
x <- 0.25
for (i in 1:nrow(Plate2)) {
  Plate2[i,"Time"] <- x
  x <- x + 0.25
}
#Limit this Plate2 to 48hrs
Plate2 <- Plate2[Plate2$Time <= 48,]
#Take a quick look at the raw data for all wells
matplot(Plate2[,1:96],type = "l")

#Read in the Plate2 setup information
Plate2key <- read.xlsx(paste0(data.dir,"2026_01_30_GCPlatekey_MCF.xlsx"), colNames = TRUE) 
#Subtract the blanks from all wells
Blank <- mean(as.matrix(Plate2[,Plate2key[Plate2key$Strain == "blank","Well"]]))
Blank
Plate2[,1:96] <- Plate2[,1:96] - Blank

Plate2melt <- data.table(Plate2[,1:96])
Plate2melt <- melt(Plate2melt)
colnames(Plate2melt)[1] <- "Well"
Plate2melt <- left_join(Plate2melt, Plate2key, by = "Well")
#Need to add time back in 
Plate2melt$Time <- rep(Plate2$Time, 96)

#Create a time column - measurements were taken every 15min
x <- 0.25
for (i in 1:nrow(Plate3)) {
  Plate3[i,"Time"] <- x
  x <- x + 0.25
}
#Limit this Plate3 to 48hrs
Plate3 <- Plate3[Plate3$Time <= 48,]
#Take a quick look at the raw data for all wells
matplot(Plate3[,1:96],type = "l")

#Read in the Plate3 setup information
Plate3key <- read.xlsx(paste0(data.dir,"2026_01_30_GCPlatekey_MCF.xlsx"), colNames = TRUE) 
#Subtract the blanks from all wells
Blank <- mean(as.matrix(Plate3[,Plate3key[Plate3key$Strain == "blank","Well"]]))
Blank
Plate3[,1:96] <- Plate3[,1:96] - Blank

Plate3melt <- data.table(Plate3[,1:96])
Plate3melt <- melt(Plate3melt)
colnames(Plate3melt)[1] <- "Well"
Plate3melt <- left_join(Plate3melt, Plate3key, by = "Well")
#Need to add time back in 
Plate3melt$Time <- rep(Plate3$Time, 96)

#Joining each plate together
Platemelt <- rbind(Plate1melt, Plate2melt, Plate3melt)

#Factorize the concentrations
Platemelt$Condition <- factor(Platemelt$Condition, levels = c("YPAD","0.004µg/ml","0.008µg/ml","0.016µg/ml","0.032µg/ml","0.064µg/ml","0.128µg/ml"))

#Plotting
a <- ggplot(data = Platemelt[Platemelt$Strain %in% c("5192.0","5192_C3_03370C_C31","5192_C3_03370C_C32"),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  facet_grid(~Condition) +
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
  scale_color_manual(values = c("cyan3","cyan3","black")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

b <- ggplot(data = Platemelt[Platemelt$Strain %in% c("P75016","P75016_C3_03370C_C11","P75016_C3_03370C_C12"),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  facet_grid(~Condition) +
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
  scale_color_manual(values = c("black","mediumvioletred","mediumvioletred")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

c <- ggplot(data = Platemelt[Platemelt$Strain %in% c("P75063","P75063_C3_03370C_C24","P75063_C3_03370C_C26"),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  facet_grid(~Condition) +
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
  scale_color_manual(values = c("black","goldenrod3","goldenrod3")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

d <- ggplot(data = Platemelt[Platemelt$Strain %in% c("L26","L26_C3_03370C_C1","L26_C3_03370C_C2"),], aes(x = Time, y = value, color=Strain, group = Strain)) +
  facet_grid(~Condition) +
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
  scale_color_manual(values = c("black","forestgreen","forestgreen")) +
  xlab("Time (hours)") +
  ylab("Growth (OD600)") +
  scale_y_continuous(breaks=c(0,0.5,0.75,1.0,1.25)) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

plot_grid(a,d,c,b, nrow = 4)
ggsave("C3MCFGC.tiff", path = fig.dir, plot = last_plot(), width = 12, height = 10)

