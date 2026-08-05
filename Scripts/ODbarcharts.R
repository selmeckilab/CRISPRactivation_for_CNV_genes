###############
## Coding for plotting the relative OD values as bar charts. 
## Raw data and analysis generating the relative OD values are present in the associated spreadsheets drawn from by this script.
##############

#Clean environment
rm(list = ls())
#Load libraries and set working directory
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
library(ggplot2)
library(openxlsx)
library(tidyr)
library(dplyr)
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
#Set directories to read in processed excel files and directory to write figures
data.dir <- "Data/"
fig.dir <- "Figures/Updated/"

#SET3 MCF
SET3MCF <- read.xlsx(paste0(data.dir,"2026_04_01_MIC_SET3_MCF.xlsx"), sheet = 1, rows = c(43,47,56,65), cols = c(3:15), colNames = TRUE)
SET3MCF <- pivot_longer(SET3MCF, cols = colnames(SET3MCF), values_to = "Relative_OD", names_to = "Strain")
SET3MCF$Mutant <- rep(c("L26","L26 SET3 OE 1","L26 SET3 OE 2","P75016","P75016 SET3 OE 1","P75016 SET3 OE 2","P75063","P75063 SET3 OE 1","P75063 SET3 OE 2","AMS5192","AMS5192 SET3 OE 1","AMS5192 SET3 OE 2"), 3)
t.test(SET3MCF[SET3MCF$Mutant == "L26","Relative_OD"], SET3MCF[SET3MCF$Mutant == "L26 SET3 OE 1","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "L26","Relative_OD"], SET3MCF[SET3MCF$Mutant == "L26 SET3 OE 2","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "P75016","Relative_OD"], SET3MCF[SET3MCF$Mutant == "P75016 SET3 OE 1","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "P75016","Relative_OD"], SET3MCF[SET3MCF$Mutant == "P75016 SET3 OE 2","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "P75063","Relative_OD"], SET3MCF[SET3MCF$Mutant == "P75063 SET3 OE 1","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "P75063","Relative_OD"], SET3MCF[SET3MCF$Mutant == "P75063 SET3 OE 2","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "AMS5192","Relative_OD"], SET3MCF[SET3MCF$Mutant == "AMS5192 SET3 OE 1","Relative_OD"])
t.test(SET3MCF[SET3MCF$Mutant == "AMS5192","Relative_OD"], SET3MCF[SET3MCF$Mutant == "AMS5192 SET3 OE 2","Relative_OD"])

SET3MCF <- SET3MCF %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
SET3MCFSUM <- unique(SET3MCF[,c("Mutant","Mean","SD")])
SET3MCFSUM$Mutant <- factor(SET3MCFSUM$Mutant, levels = c("AMS5192","AMS5192 SET3 OE 1","AMS5192 SET3 OE 2", "L26","L26 SET3 OE 1","L26 SET3 OE 2","P75063","P75063 SET3 OE 1","P75063 SET3 OE 2","P75016","P75016 SET3 OE 1","P75016 SET3 OE 2"))

ggplot(data = SET3MCFSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = SET3MCFSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("SET3"), " OE in MCF")))
ggsave("SET3MCFBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#SET3 SDS
SET3SDS <- read.xlsx(paste0(data.dir,"2026_04_09_MIC_SET3_SDS.xlsx"), sheet = 1, rows = c(43,48,57,66), cols = c(3:15), colNames = TRUE)
SET3SDS <- pivot_longer(SET3SDS, cols = colnames(SET3SDS), values_to = "Relative_OD", names_to = "Strain")
SET3SDS$Mutant <- rep(c("L26","L26 SET3 OE 1","L26 SET3 OE 2","P75016","P75016 SET3 OE 1","P75016 SET3 OE 2","P75063","P75063 SET3 OE 1","P75063 SET3 OE 2","AMS5192","AMS5192 SET3 OE 1","AMS5192 SET3 OE 2"), 3)
t.test(SET3SDS[SET3SDS$Mutant == "L26","Relative_OD"], SET3SDS[SET3SDS$Mutant == "L26 SET3 OE 1","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "L26","Relative_OD"], SET3SDS[SET3SDS$Mutant == "L26 SET3 OE 2","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "P75016","Relative_OD"], SET3SDS[SET3SDS$Mutant == "P75016 SET3 OE 1","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "P75016","Relative_OD"], SET3SDS[SET3SDS$Mutant == "P75016 SET3 OE 2","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "P75063","Relative_OD"], SET3SDS[SET3SDS$Mutant == "P75063 SET3 OE 1","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "P75063","Relative_OD"], SET3SDS[SET3SDS$Mutant == "P75063 SET3 OE 2","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "AMS5192","Relative_OD"], SET3SDS[SET3SDS$Mutant == "AMS5192 SET3 OE 1","Relative_OD"])
t.test(SET3SDS[SET3SDS$Mutant == "AMS5192","Relative_OD"], SET3SDS[SET3SDS$Mutant == "AMS5192 SET3 OE 2","Relative_OD"])
SET3SDS <- SET3SDS %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
SET3SDSSUM <- unique(SET3SDS[,c("Mutant","Mean","SD")])
SET3SDSSUM$Mutant <- factor(SET3SDSSUM$Mutant, levels = c("AMS5192","AMS5192 SET3 OE 1","AMS5192 SET3 OE 2", "L26","L26 SET3 OE 1","L26 SET3 OE 2","P75063","P75063 SET3 OE 1","P75063 SET3 OE 2","P75016","P75016 SET3 OE 1","P75016 SET3 OE 2"))

ggplot(data = SET3SDSSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = SET3SDSSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("SET3"), " OE in SDS")))
ggsave("SET3SDSBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#SET3 FLC - this one is a bit trickier because I need different rows for L26 than for the other three do to a different progenitor MIC50.
SET3FLC <- read.xlsx(paste0(data.dir,"2026_04_09_MIC_SET3_FLC.xlsx"), sheet = 1, rows = c(43,49,58,67), cols = c(3:5), colNames = TRUE)
SET3FLC2 <- read.xlsx(paste0(data.dir,"2026_04_09_MIC_SET3_FLC.xlsx"), sheet = 1, rows = c(43,46,55,64), cols = c(6:15), colNames = TRUE)
SET3FLC <- cbind(SET3FLC, SET3FLC2)
SET3FLC <- pivot_longer(SET3FLC, cols = colnames(SET3FLC), values_to = "Relative_OD", names_to = "Strain")
SET3FLC$Mutant <- rep(c("L26","L26 SET3 OE 1","L26 SET3 OE 2","P75016","P75016 SET3 OE 1","P75016 SET3 OE 2","P75063","P75063 SET3 OE 1","P75063 SET3 OE 2","AMS5192","AMS5192 SET3 OE 1","AMS5192 SET3 OE 2"), 3)
t.test(SET3FLC[SET3FLC$Mutant == "L26","Relative_OD"], SET3FLC[SET3FLC$Mutant == "L26 SET3 OE 1","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "L26","Relative_OD"], SET3FLC[SET3FLC$Mutant == "L26 SET3 OE 2","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "P75016","Relative_OD"], SET3FLC[SET3FLC$Mutant == "P75016 SET3 OE 1","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "P75016","Relative_OD"], SET3FLC[SET3FLC$Mutant == "P75016 SET3 OE 2","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "P75063","Relative_OD"], SET3FLC[SET3FLC$Mutant == "P75063 SET3 OE 1","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "P75063","Relative_OD"], SET3FLC[SET3FLC$Mutant == "P75063 SET3 OE 2","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "AMS5192","Relative_OD"], SET3FLC[SET3FLC$Mutant == "AMS5192 SET3 OE 1","Relative_OD"])
t.test(SET3FLC[SET3FLC$Mutant == "AMS5192","Relative_OD"], SET3FLC[SET3FLC$Mutant == "AMS5192 SET3 OE 2","Relative_OD"])

SET3FLC <- SET3FLC %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
SET3FLCSUM <- unique(SET3FLC[,c("Mutant","Mean","SD")])
SET3FLCSUM$Mutant <- factor(SET3FLCSUM$Mutant, levels = c("AMS5192","AMS5192 SET3 OE 1","AMS5192 SET3 OE 2", "L26","L26 SET3 OE 1","L26 SET3 OE 2","P75063","P75063 SET3 OE 1","P75063 SET3 OE 2","P75016","P75016 SET3 OE 1","P75016 SET3 OE 2"))

ggplot(data = SET3FLCSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = SET3FLCSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("SET3"), " OE in FLC")))
ggsave("SET3FLCBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#C303460C in MCF - We had to redo one background, so I have to pull from a two different places. 
C303MCF <- read.xlsx(paste0(data.dir,"2025_10_08_MIC_C303_MCF.xlsx"), sheet = 1, rows = c(55,59,68,77), cols = c(3:11), colNames = TRUE)
C303MCFRedo <- read.xlsx(paste0(data.dir,"2025_12_18_TROUBLECRISPRA_MCF_SDS.xlsx"), sheet = 2, rows = c(43,47,55,63), cols = c(16:18), colNames = TRUE)
C303MCF <- cbind(C303MCF, C303MCFRedo)
C303MCF <- pivot_longer(C303MCF, cols = colnames(C303MCF), values_to = "Relative_OD", names_to = "Strain")
C303MCF$Mutant <- rep(c("L26","L26 C303 OE 1","L26 C303 OE 2","P75016","P75016 C303 OE 1","P75016 C303 OE 2","P75063","P75063 C303 OE 1","P75063 C303 OE 2","AMS5192","AMS5192 C303 OE 1","AMS5192 C303 OE 2"), 3)
t.test(C303MCF[C303MCF$Mutant == "L26","Relative_OD"], C303MCF[C303MCF$Mutant == "L26 C303 OE 1","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "L26","Relative_OD"], C303MCF[C303MCF$Mutant == "L26 C303 OE 2","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "P75016","Relative_OD"], C303MCF[C303MCF$Mutant == "P75016 C303 OE 1","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "P75016","Relative_OD"], C303MCF[C303MCF$Mutant == "P75016 C303 OE 2","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "P75063","Relative_OD"], C303MCF[C303MCF$Mutant == "P75063 C303 OE 1","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "P75063","Relative_OD"], C303MCF[C303MCF$Mutant == "P75063 C303 OE 2","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "AMS5192","Relative_OD"], C303MCF[C303MCF$Mutant == "AMS5192 C303 OE 1","Relative_OD"])
t.test(C303MCF[C303MCF$Mutant == "AMS5192","Relative_OD"], C303MCF[C303MCF$Mutant == "AMS5192 C303 OE 2","Relative_OD"])


C303MCF <- C303MCF %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
C303MCFSUM <- unique(C303MCF[,c("Mutant","Mean","SD")])
C303MCFSUM$Mutant <- factor(C303MCFSUM$Mutant, levels = c("AMS5192","AMS5192 C303 OE 1","AMS5192 C303 OE 2", "L26","L26 C303 OE 1","L26 C303 OE 2","P75063","P75063 C303 OE 1","P75063 C303 OE 2","P75016","P75016 C303 OE 1","P75016 C303 OE 2"))

ggplot(data = C303MCFSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = C303MCFSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("C3_03460C"), " OE in MCF")))
ggsave("C303MCFBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#C303 in SDS - also taking different concentrations because of different progenitor MICs
C303SDS <- read.xlsx(paste0(data.dir,"2026_04_23_MIC_C303_SDS.xlsx"), sheet = 1, rows = c(34,44,54), cols = c(2:7), colNames = FALSE)
C303SDS2 <- read.xlsx(paste0(data.dir,"2026_04_23_MIC_C303_SDS.xlsx"), sheet = 1, rows = c(35,45,55), cols = c(8:13), colNames = FALSE)
C303SDS <- cbind(C303SDS, C303SDS2)
C303SDS <- pivot_longer(C303SDS, cols = colnames(C303SDS), values_to = "Relative_OD", names_to = "Strain")
C303SDS$Mutant <- rep(c("L26","L26 C303 OE 1","L26 C303 OE 2","P75016","P75016 C303 OE 1","P75016 C303 OE 2","P75063","P75063 C303 OE 1","P75063 C303 OE 2","AMS5192","AMS5192 C303 OE 1","AMS5192 C303 OE 2"), 3)
t.test(C303SDS[C303SDS$Mutant == "L26","Relative_OD"], C303SDS[C303SDS$Mutant == "L26 C303 OE 1","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "L26","Relative_OD"], C303SDS[C303SDS$Mutant == "L26 C303 OE 2","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "P75016","Relative_OD"], C303SDS[C303SDS$Mutant == "P75016 C303 OE 1","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "P75016","Relative_OD"], C303SDS[C303SDS$Mutant == "P75016 C303 OE 2","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "P75063","Relative_OD"], C303SDS[C303SDS$Mutant == "P75063 C303 OE 1","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "P75063","Relative_OD"], C303SDS[C303SDS$Mutant == "P75063 C303 OE 2","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "AMS5192","Relative_OD"], C303SDS[C303SDS$Mutant == "AMS5192 C303 OE 1","Relative_OD"])
t.test(C303SDS[C303SDS$Mutant == "AMS5192","Relative_OD"], C303SDS[C303SDS$Mutant == "AMS5192 C303 OE 2","Relative_OD"])

C303SDS <- C303SDS %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
C303SDSSUM <- unique(C303SDS[,c("Mutant","Mean","SD")])
C303SDSSUM$Mutant <- factor(C303SDSSUM$Mutant, levels = c("AMS5192","AMS5192 C303 OE 1","AMS5192 C303 OE 2", "L26","L26 C303 OE 1","L26 C303 OE 2","P75063","P75063 C303 OE 1","P75063 C303 OE 2","P75016","P75016 C303 OE 1","P75016 C303 OE 2"))

ggplot(data = C303SDSSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = C303SDSSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("C3_03460C"), " OE in SDS")))
ggsave("C303SDSBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#C303 FLC 
C303FLC <- read.xlsx(paste0(data.dir,"2026_03_25_MIC_C303_FLC.xlsx"), sheet = 1, rows = c(43,49,58,67), cols = c(3:5), colNames = TRUE)
C303FLC2 <- read.xlsx(paste0(data.dir,"2026_03_25_MIC_C303_FLC.xlsx"), sheet = 1, rows = c(43,47,56,65), cols = c(6:15), colNames = TRUE)
C303FLC <- cbind(C303FLC, C303FLC2)
C303FLC <- pivot_longer(C303FLC, cols = colnames(C303FLC), values_to = "Relative_OD", names_to = "Strain")
C303FLC$Mutant <- rep(c("L26","L26 C303 OE 1","L26 C303 OE 2","P75016","P75016 C303 OE 1","P75016 C303 OE 2","P75063","P75063 C303 OE 1","P75063 C303 OE 2","AMS5192","AMS5192 C303 OE 1","AMS5192 C303 OE 2"), 3)
t.test(C303FLC[C303FLC$Mutant == "L26","Relative_OD"], C303FLC[C303FLC$Mutant == "L26 C303 OE 1","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "L26","Relative_OD"], C303FLC[C303FLC$Mutant == "L26 C303 OE 2","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "P75016","Relative_OD"], C303FLC[C303FLC$Mutant == "P75016 C303 OE 1","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "P75016","Relative_OD"], C303FLC[C303FLC$Mutant == "P75016 C303 OE 2","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "P75063","Relative_OD"], C303FLC[C303FLC$Mutant == "P75063 C303 OE 1","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "P75063","Relative_OD"], C303FLC[C303FLC$Mutant == "P75063 C303 OE 2","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "AMS5192","Relative_OD"], C303FLC[C303FLC$Mutant == "AMS5192 C303 OE 1","Relative_OD"])
t.test(C303FLC[C303FLC$Mutant == "AMS5192","Relative_OD"], C303FLC[C303FLC$Mutant == "AMS5192 C303 OE 2","Relative_OD"])
C303FLC <- C303FLC %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
C303FLCSUM <- unique(C303FLC[,c("Mutant","Mean","SD")])
C303FLCSUM$Mutant <- factor(C303FLCSUM$Mutant, levels = c("AMS5192","AMS5192 C303 OE 1","AMS5192 C303 OE 2", "L26","L26 C303 OE 1","L26 C303 OE 2","P75063","P75063 C303 OE 1","P75063 C303 OE 2","P75016","P75016 C303 OE 1","P75016 C303 OE 2"))

ggplot(data = C303FLCSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = C303FLCSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("C3_03460C"), " OE in FLC")))
ggsave("C303FLCBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  


#SKN1 in MCF
SKN1MCF <- read.xlsx(paste0(data.dir,"2026_04_01_MIC_SKN1_MCF.xlsx"), sheet = 1, rows = c(43,47,56,65), cols = c(3:15), colNames = TRUE)
SKN1MCF <- pivot_longer(SKN1MCF, cols = colnames(SKN1MCF), values_to = "Relative_OD", names_to = "Strain")
SKN1MCF$Mutant <- rep(c("L26","L26 SKN1 OE 1","L26 SKN1 OE 2","P75016","P75016 SKN1 OE 1","P75016 SKN1 OE 2","P75063","P75063 SKN1 OE 1","P75063 SKN1 OE 2","AMS5192","AMS5192 SKN1 OE 1","AMS5192 SKN1 OE 2"), 3)
t.test(SKN1MCF[SKN1MCF$Mutant == "L26","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "L26 SKN1 OE 1","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "L26","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "L26 SKN1 OE 2","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "P75016","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "P75016 SKN1 OE 1","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "P75016","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "P75016 SKN1 OE 2","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "P75063","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "P75063 SKN1 OE 1","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "P75063","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "P75063 SKN1 OE 2","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "AMS5192","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "AMS5192 SKN1 OE 1","Relative_OD"])
t.test(SKN1MCF[SKN1MCF$Mutant == "AMS5192","Relative_OD"], SKN1MCF[SKN1MCF$Mutant == "AMS5192 SKN1 OE 2","Relative_OD"])
SKN1MCF <- SKN1MCF %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
SKN1MCFSUM <- unique(SKN1MCF[,c("Mutant","Mean","SD")])
SKN1MCFSUM$Mutant <- factor(SKN1MCFSUM$Mutant, levels = c("AMS5192","AMS5192 SKN1 OE 1","AMS5192 SKN1 OE 2", "L26","L26 SKN1 OE 1","L26 SKN1 OE 2","P75063","P75063 SKN1 OE 1","P75063 SKN1 OE 2","P75016","P75016 SKN1 OE 1","P75016 SKN1 OE 2"))

ggplot(data = SKN1MCFSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = SKN1MCFSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("SKN1"), " OE in MCF")))
ggsave("SKN1MCFBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  



#C3_03370C in FLC 
C3_03370CFLC <- read.xlsx(paste0(data.dir,"2025_11_06_MIC_SKN12_C3.xlsx"), sheet = 1, rows = c(55,62,71,80), cols = c(3:5), colNames = TRUE)
C3_03370CFLC2 <- read.xlsx(paste0(data.dir,"2025_11_06_MIC_SKN12_C3.xlsx"), sheet = 1, rows = c(55,59,68,77), cols = c(6:15), colNames = TRUE)
C3_03370CFLC <- cbind(C3_03370CFLC, C3_03370CFLC2)
C3_03370CFLC <- pivot_longer(C3_03370CFLC, cols = colnames(C3_03370CFLC), values_to = "Relative_OD", names_to = "Strain")
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "2868.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C1","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "2868.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C5","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "2875.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C9","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "2875.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C13","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "2876.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C18","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "2876.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C19","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "5192.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C25","Relative_OD"])
t.test(C3_03370CFLC[C3_03370CFLC$Strain == "5192.0","Relative_OD"], C3_03370CFLC[C3_03370CFLC$Strain == "C27","Relative_OD"])
C3_03370CFLC$Mutant <- rep(c("L26","L26 C3_03370C OE 1","L26 C3_03370C OE 2","P75016","P75016 C3_03370C OE 1","P75016 C3_03370C OE 2","P75063","P75063 C3_03370C OE 1","P75063 C3_03370C OE 2","AMS5192","AMS5192 C3_03370C OE 1","AMS5192 C3_03370C OE 2"), 3)
C3_03370CFLC <- C3_03370CFLC %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
C3_03370CFLCSUM <- unique(C3_03370CFLC[,c("Mutant","Mean","SD")])
C3_03370CFLCSUM$Mutant <- factor(C3_03370CFLCSUM$Mutant, levels = c("AMS5192","AMS5192 C3_03370C OE 1","AMS5192 C3_03370C OE 2", "L26","L26 C3_03370C OE 1","L26 C3_03370C OE 2","P75063","P75063 C3_03370C OE 1","P75063 C3_03370C OE 2","P75016","P75016 C3_03370C OE 1","P75016 C3_03370C OE 2"))

ggplot(data = C3_03370CFLCSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = C3_03370CFLCSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media")
  #ggtitle(expression(paste(italic("C3_03370C"), " OE in FLC")))
ggsave("C3_03370CFLCBar.tiff", path = fig.dir, width = 3, height = 3.5, plot = last_plot())  


#TYE7 in FLC
TYE7FLC <- read.xlsx(paste0(data.dir,"2025_10_02_MIC_TYE7_FLC.xlsx"), sheet = 1, rows = c(55,62,71,80), cols = c(3:5), colNames = TRUE)
TYE7FLC2 <- read.xlsx(paste0(data.dir,"2025_10_02_MIC_TYE7_FLC.xlsx"), sheet = 1, rows = c(55,60,69,78), cols = c(6:15), colNames = TRUE)
TYE7FLC <- cbind(TYE7FLC, TYE7FLC2)
TYE7FLC <- pivot_longer(TYE7FLC, cols = colnames(TYE7FLC), values_to = "Relative_OD", names_to = "Strain")
TYE7FLC$Mutant <- rep(c("L26","L26 TYE7 OE 1","L26 TYE7 OE 2","P75016","P75016 TYE7 OE 1","P75016 TYE7 OE 2","P75063","P75063 TYE7 OE 1","P75063 TYE7 OE 2","AMS5192","AMS5192 TYE7 OE 1","AMS5192 TYE7 OE 2"), 3)
t.test(TYE7FLC[TYE7FLC$Mutant == "L26","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "L26 TYE7 OE 1","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "L26","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "L26 TYE7 OE 2","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "P75016","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "P75016 TYE7 OE 1","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "P75016","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "P75016 TYE7 OE 2","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "P75063","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "P75063 TYE7 OE 1","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "P75063","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "P75063 TYE7 OE 2","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "AMS5192","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "AMS5192 TYE7 OE 1","Relative_OD"])
t.test(TYE7FLC[TYE7FLC$Mutant == "AMS5192","Relative_OD"], TYE7FLC[TYE7FLC$Mutant == "AMS5192 TYE7 OE 2","Relative_OD"])
TYE7FLC <- TYE7FLC %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
TYE7FLCSUM <- unique(TYE7FLC[,c("Mutant","Mean","SD")])
TYE7FLCSUM$Mutant <- factor(TYE7FLCSUM$Mutant, levels = c("AMS5192","AMS5192 TYE7 OE 1","AMS5192 TYE7 OE 2", "L26","L26 TYE7 OE 1","L26 TYE7 OE 2","P75063","P75063 TYE7 OE 1","P75063 TYE7 OE 2","P75016","P75016 TYE7 OE 1","P75016 TYE7 OE 2"))

ggplot(data = TYE7FLCSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = TYE7FLCSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media")
  #ggtitle(expression(paste(italic("TYE7"), " OE at progenitor MIC")))
ggsave("TYE7FLCBar.tiff", path = fig.dir, width = 3, height = 3.5, plot = last_plot())  

#C111 in MCF - one mutant was not in agreement and was redone
C111MCF <- read.xlsx(paste0(data.dir,"2025_09_25_MIC_C111_MCF.xlsx"), sheet = 1, rows = c(45, 54, 63), cols = c(3:8), colNames = FALSE)
C111MCF2 <- read.xlsx(paste0(data.dir,"2026_05_14_MIC_C111_MCF.xlsx"), sheet = 1, rows = c(34,44,54), cols = c(8:10), colNames = FALSE)
C111MCF3 <- read.xlsx(paste0(data.dir,"2025_09_25_MIC_C111_MCF.xlsx"), sheet = 1, rows = c(45, 54, 63), cols = c(12:14), colNames = FALSE)
C111MCF <- cbind(C111MCF, C111MCF2, C111MCF3)
C111MCF <- pivot_longer(C111MCF, cols = colnames(C111MCF), values_to = "Relative_OD", names_to = "Strain")
C111MCF$Mutant <- rep(c("L26","L26 C111 OE 1","L26 C111 OE 2","P75016","P75016 C111 OE 1","P75016 C111 OE 2","P75063","P75063 C111 OE 1","P75063 C111 OE 2","AMS5192","AMS5192 C111 OE 1","AMS5192 C111 OE 2"), 3)
t.test(C111MCF[C111MCF$Mutant == "L26","Relative_OD"], C111MCF[C111MCF$Mutant == "L26 C111 OE 1","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "L26","Relative_OD"], C111MCF[C111MCF$Mutant == "L26 C111 OE 2","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "P75016","Relative_OD"], C111MCF[C111MCF$Mutant == "P75016 C111 OE 1","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "P75016","Relative_OD"], C111MCF[C111MCF$Mutant == "P75016 C111 OE 2","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "P75063","Relative_OD"], C111MCF[C111MCF$Mutant == "P75063 C111 OE 1","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "P75063","Relative_OD"], C111MCF[C111MCF$Mutant == "P75063 C111 OE 2","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "AMS5192","Relative_OD"], C111MCF[C111MCF$Mutant == "AMS5192 C111 OE 1","Relative_OD"])
t.test(C111MCF[C111MCF$Mutant == "AMS5192","Relative_OD"], C111MCF[C111MCF$Mutant == "AMS5192 C111 OE 2","Relative_OD"])
C111MCF <- C111MCF %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
C111MCFSUM <- unique(C111MCF[,c("Mutant","Mean","SD")])
C111MCFSUM$Mutant <- factor(C111MCFSUM$Mutant, levels = c("AMS5192","AMS5192 C111 OE 1","AMS5192 C111 OE 2", "L26","L26 C111 OE 1","L26 C111 OE 2","P75063","P75063 C111 OE 1","P75063 C111 OE 2","P75016","P75016 C111 OE 1","P75016 C111 OE 2"))

ggplot(data = C111MCFSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = C111MCFSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("C1_11720W"), " OE in MCF")))
ggsave("C111MCFBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#YCK2 in MCF 
YCK2MCF <- read.xlsx(paste0(data.dir,"2026_04_22_MIC_YCK2_MCF.xlsx"), sheet = 1, rows = c(34,44,54), cols = c(2:14), colNames = FALSE)
YCK2MCF <- pivot_longer(YCK2MCF, cols = colnames(YCK2MCF), values_to = "Relative_OD", names_to = "Strain")
YCK2MCF$Mutant <- rep(c("L26","L26 YCK2 OE 1","L26 YCK2 OE 2","P75016","P75016 YCK2 OE 1","P75016 YCK2 OE 2","P75063","P75063 YCK2 OE 1","P75063 YCK2 OE 2","AMS5192","AMS5192 YCK2 OE 1","AMS5192 YCK2 OE 2"), 3)
t.test(YCK2MCF[YCK2MCF$Mutant == "L26","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "L26 YCK2 OE 1","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "L26","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "L26 YCK2 OE 2","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "P75016","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "P75016 YCK2 OE 1","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "P75016","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "P75016 YCK2 OE 2","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "P75063","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "P75063 YCK2 OE 1","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "P75063","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "P75063 YCK2 OE 2","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "AMS5192","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "AMS5192 YCK2 OE 1","Relative_OD"])
t.test(YCK2MCF[YCK2MCF$Mutant == "AMS5192","Relative_OD"], YCK2MCF[YCK2MCF$Mutant == "AMS5192 YCK2 OE 2","Relative_OD"])
YCK2MCF <- YCK2MCF %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
YCK2MCFSUM <- unique(YCK2MCF[,c("Mutant","Mean","SD")])
YCK2MCFSUM$Mutant <- factor(YCK2MCFSUM$Mutant, levels = c("AMS5192","AMS5192 YCK2 OE 1","AMS5192 YCK2 OE 2", "L26","L26 YCK2 OE 1","L26 YCK2 OE 2","P75063","P75063 YCK2 OE 1","P75063 YCK2 OE 2","P75016","P75016 YCK2 OE 1","P75016 YCK2 OE 2"))

ggplot(data = YCK2MCFSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = YCK2MCFSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("YCK2"), " OE in MCF")))
ggsave("YCK2MCFBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#YCK2 in SDS
YCK2SDS <- read.xlsx(paste0(data.dir,"2026_04_15_MIC_YCK2_SDS.xlsx"), sheet = 1, rows = c(22, 27, 36, 45), cols = c(3:5), colNames = TRUE)
YCK2SDS2 <- read.xlsx(paste0(data.dir,"2026_04_15_MIC_YCK2_SDS.xlsx"), sheet = 1, rows = c(22, 26, 35, 44), cols = c(6:8), colNames = TRUE)
YCK2SDS3 <- read.xlsx(paste0(data.dir,"2026_04_15_MIC_YCK2_SDS.xlsx"), sheet = 1, rows = c(22, 27, 36, 45), cols = c(9:15), colNames = TRUE)
YCK2SDS <- cbind(YCK2SDS, YCK2SDS2, YCK2SDS3)
YCK2SDS <- pivot_longer(YCK2SDS, cols = colnames(YCK2SDS), values_to = "Relative_OD", names_to = "Strain")
YCK2SDS$Mutant <- rep(c("L26","L26 YCK2 OE 1","L26 YCK2 OE 2","P75016","P75016 YCK2 OE 1","P75016 YCK2 OE 2","P75063","P75063 YCK2 OE 1","P75063 YCK2 OE 2","AMS5192","AMS5192 YCK2 OE 1","AMS5192 YCK2 OE 2"), 3)
t.test(YCK2SDS[YCK2SDS$Mutant == "L26","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "L26 YCK2 OE 1","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "L26","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "L26 YCK2 OE 2","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "P75016","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "P75016 YCK2 OE 1","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "P75016","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "P75016 YCK2 OE 2","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "P75063","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "P75063 YCK2 OE 1","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "P75063","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "P75063 YCK2 OE 2","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "AMS5192","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "AMS5192 YCK2 OE 1","Relative_OD"])
t.test(YCK2SDS[YCK2SDS$Mutant == "AMS5192","Relative_OD"], YCK2SDS[YCK2SDS$Mutant == "AMS5192 YCK2 OE 2","Relative_OD"])
YCK2SDS <- YCK2SDS %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
YCK2SDSSUM <- unique(YCK2SDS[,c("Mutant","Mean","SD")])
YCK2SDSSUM$Mutant <- factor(YCK2SDSSUM$Mutant, levels = c("AMS5192","AMS5192 YCK2 OE 1","AMS5192 YCK2 OE 2", "L26","L26 YCK2 OE 1","L26 YCK2 OE 2","P75063","P75063 YCK2 OE 1","P75063 YCK2 OE 2","P75016","P75016 YCK2 OE 1","P75016 YCK2 OE 2"))

ggplot(data = YCK2SDSSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = YCK2SDSSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  ylab("Growth (OD600) \nRelative to rich media") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ggtitle(expression(paste(italic("YCK2"), " OE in SDS")))
ggsave("YCK2SDSBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  

#YCK2 FLC 
YCK2FLC <- read.xlsx(paste0(data.dir,"2026_04_22_MIC_YCK2_FLC.xlsx"), sheet = 1, rows = c(37, 47, 57), cols = c(2:4), colNames = FALSE)
YCK2FLC2 <- read.xlsx(paste0(data.dir,"2026_04_22_MIC_YCK2_FLC.xlsx"), sheet = 1, rows = c(34, 44, 54), cols = c(5:14), colNames = FALSE)
YCK2FLC <- cbind(YCK2FLC, YCK2FLC2)
YCK2FLC <- pivot_longer(YCK2FLC, cols = colnames(YCK2FLC), values_to = "Relative_OD", names_to = "Strain")
YCK2FLC$Mutant <- rep(c("L26","L26 YCK2 OE 1","L26 YCK2 OE 2","P75016","P75016 YCK2 OE 1","P75016 YCK2 OE 2","P75063","P75063 YCK2 OE 1","P75063 YCK2 OE 2","AMS5192","AMS5192 YCK2 OE 1","AMS5192 YCK2 OE 2"), 3)
t.test(YCK2FLC[YCK2FLC$Mutant == "L26","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "L26 YCK2 OE 1","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "L26","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "L26 YCK2 OE 2","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "P75016","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "P75016 YCK2 OE 1","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "P75016","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "P75016 YCK2 OE 2","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "P75063","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "P75063 YCK2 OE 1","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "P75063","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "P75063 YCK2 OE 2","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "AMS5192","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "AMS5192 YCK2 OE 1","Relative_OD"])
t.test(YCK2FLC[YCK2FLC$Mutant == "AMS5192","Relative_OD"], YCK2FLC[YCK2FLC$Mutant == "AMS5192 YCK2 OE 2","Relative_OD"])
YCK2FLC <- YCK2FLC %>%
  group_by(Mutant) %>%
  mutate(Mean = mean(Relative_OD), SD = sd(Relative_OD))
YCK2FLCSUM <- unique(YCK2FLC[,c("Mutant","Mean","SD")])
YCK2FLCSUM$Mutant <- factor(YCK2FLCSUM$Mutant, levels = c("AMS5192","AMS5192 YCK2 OE 1","AMS5192 YCK2 OE 2", "L26","L26 YCK2 OE 1","L26 YCK2 OE 2","P75063","P75063 YCK2 OE 1","P75063 YCK2 OE 2","P75016","P75016 YCK2 OE 1","P75016 YCK2 OE 2"))

ggplot(data = YCK2FLCSUM, aes(x = Mutant, y = Mean)) + 
  geom_col(aes(fill = Mutant)) +
  scale_fill_manual(values = c("grey","cyan3","cyan3","grey", "forestgreen", "forestgreen","grey","goldenrod3","goldenrod3","grey","mediumorchid","mediumorchid")) +
  geom_linerange(data = YCK2FLCSUM, aes(ymin = Mean - SD, ymax = Mean + SD)) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab("Growth (OD600) \nRelative to rich media") +
  ggtitle(expression(paste(italic("YCK2"), " OE in FLC")))
ggsave("YCK2FLCBar.tiff", path = fig.dir, width = 3, height = 3.75, plot = last_plot())  
