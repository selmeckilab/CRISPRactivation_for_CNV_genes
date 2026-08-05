#########################
## Analysis and figure generation for RT-qPCR and MICs of control strains
## overexpressing CDR1 via CRISPR-activation
########################

#Setting up the environment
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
library(dplyr)
library(ggplot2)

#Make sure there is nothing lurking in the environment
rm(list = ls())

#Read in data file
data.dir <- "Data/"
fig.dir <- "Figures/"
Platedata <- read.table(paste0(data.dir,"20250610run -  Quantification Cq Results_0.txt"), sep = "\t", header = 1)
Platedata <- Platedata[,c("Well","Cq")]
Platedata <- Platedata[1:96,]
#Read in sample information
Platekey <- read.table(paste0(data.dir,"Key_2025_06_06.txt"), sep = "\t", header = 1) #Same key as last time

Platedata <- left_join(Platedata, Platekey, by = "Well")

#Calculating average and standard deviation for each strain/primer pair
Platesum <- Platedata %>%
  group_by(Background, Strain, Primer) %>%
  summarize(Mean = mean(Cq), StDev = sd(Cq))

#Calculating deltaCt for each primer relative to TEF1 (Gene - TEF1) and error (sqrt(summed squares))
for(i in 1:nrow(Platesum)) {
  Platesum[i,"DeltaCt"] <- Platesum[i,"Mean"] - Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Strain == Platesum$Strain[i] & Platesum$Primer == "TEF1","Mean"]
  Platesum[i,"DeltaCtsd"] <- sqrt((Platesum[i,"StDev"])^2 + (Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Strain == Platesum$Strain[i] & Platesum$Primer == "TEF1","StDev"])^2)
}

#Calculating deltadeltaCt for each Strain relative to the WT (Strain - WT) and error(sqrt(summed squares))
for(i in 1:nrow(Platesum)) {
  Platesum[i,"DeltaDeltaCt"] <- Platesum[i,"DeltaCt"] - Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Primer == Platesum$Primer[i] & Platesum$Strain == "CTRL","DeltaCt"]
  Platesum[i,"DeltaDeltaCtsd"] <- sqrt((Platesum[i,"DeltaCtsd"])^2 + (Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Primer == Platesum$Primer[i] & Platesum$Strain == "CTRL","DeltaCtsd"])^2)
}

#Calculating the fold change 2^-deltadeltaCt and error (plus and minus, then 2^-)
for(i in 1:nrow(Platesum)) {
  Platesum[i,"FoldChange"] <- 2^-Platesum[i,"DeltaDeltaCt"]
  Platesum[i,"Lower"] <- Platesum[i,"DeltaDeltaCt"] + Platesum[i,"DeltaDeltaCtsd"]
  Platesum[i,"Upper"] <- Platesum[i,"DeltaDeltaCt"] - Platesum[i,"DeltaDeltaCtsd"]
  Platesum[i,"UpperFC"] <- 2^-Platesum[i,"Upper"]
  Platesum[i,"LowerFC"] <- 2^-Platesum[i,"Lower"]
}

# Plotting
# Setting the order correctly and making nicer names
Platesum$Strain <- factor(Platesum$Strain, levels = c("CTRL","BF.33", "L26.33", "63.33","16.33","16.34","SC.33","SC.34"))
Labs <- c("AMS5192\nguide 1", "L26\nguide 1", "P75063\nguide 1", "P75016\nguide 1", "P75016\nguide 2")
ggplot(data = Platesum[Platesum$Primer == "CDR1" & !Platesum$Strain %in% c("CTRL","SC.33","SC.34"),], aes(x = Strain, y = FoldChange)) +
  geom_col() +
  geom_linerange(data = Platesum[Platesum$Primer == "CDR1" & !Platesum$Strain %in% c("CTRL","SC.33","SC.34"),], ymax = Platesum[Platesum$Primer == "CDR1" & !Platesum$Strain %in% c("CTRL","SC.33","SC.34"),]$UpperFC, ymin = Platesum[Platesum$Primer == "CDR1" & !Platesum$Strain %in% c("CTRL","SC.33","SC.34"),]$LowerFC) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, color = "black"), axis.text.x = element_text(color = "black")) +
  scale_x_discrete(labels = Labs) +
  ylim(0, 14.5) +
  ylab("Fold Change \nin CDR1 Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
ggsave(paste0(fig.dir,"AllCDR1RT.tiff"), plot = last_plot(), width = 2.5, height = 3)

####################
# A quick plot of the change in MIC values for these CRISPRa strains - Taken from raw data in XXXXX
MICchangevec <- c(4, 2, 2, 2, 2)
Namevec <- c("AMS5192\nguide 1","L26\nguide 1","P75063\nguide 1","P75016\nguide 1", "P75016\nguide 2")
Plotdf <- data.frame(MIC = MICchangevec, Strain = Namevec)

Plotdf$Strain <- factor(Plotdf$Strain, levels = Namevec)

ggplot(data = Plotdf, aes(x = Strain, y = MIC)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size = 12, color = "black"), axis.text.x = element_text(color = "black")) +
  ylab(expr(Fold~Change~"in"~MIC[50]))
ggsave(paste0(fig.dir,"CDR1MICchange.tiff"), plot = last_plot(), width = 2.5, height = 3)
