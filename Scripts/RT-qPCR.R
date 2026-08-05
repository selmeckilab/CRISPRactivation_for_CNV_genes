#########################
## RT-qPCR analysis and figure generation for SKN1 and C3_03370C CRISPRa Mutants
########################
# 1. SKN1 --> line 17
# 2. C3_03370C --> line 

#Setting up the environment
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
library(dplyr)
library(ggplot2)
#Clear any lurking environmental variables
rm(list = ls())
data.dir <- "Data/"
fig.dir <- "Figures/"

#####################
#SKN1 data
# Read in the data
Platedata <- read.csv(paste0(data.dir,"20260327run -  Quantification Cq Results_0.csv"), sep = ",", header = 1)
Platedata <- Platedata[,c("Well","Cq")]
#Read in sample information
Platekey <- read.table(paste0(data.dir,"Key_2026_03_27.txt"), sep = "\t", header = 1)

Platedata <- left_join(Platedata, Platekey, by = "Well")

#Calculating average and standard deviation for each strain/primer pair
Platesum <- Platedata %>%
  group_by(Background, Strain, Primer) %>%
  summarize(Mean = mean(Cq), StDev = sd(Cq))

#Calculating deltaCt for each primer relative to ACT1 (Gene - ACT1) and error (sqrt(summed squares))
for(i in 1:nrow(Platesum)) {
  Platesum[i,"DeltaCt"] <- Platesum[i,"Mean"] - Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Strain == Platesum$Strain[i] & Platesum$Primer == "ACT1","Mean"]
  Platesum[i,"DeltaCtsd"] <- sqrt((Platesum[i,"StDev"])^2 + (Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Strain == Platesum$Strain[i] & Platesum$Primer == "ACT1","StDev"])^2)
}

#Calculating deltadeltaCt for each Strain relative to the WT (Strain - WT) and error(sqrt(summed squares))
for(i in 1:nrow(Platesum)) {
  Platesum[i,"DeltaDeltaCt"] <- Platesum[i,"DeltaCt"] - Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Primer == Platesum$Primer[i] & Platesum$Strain == "CTRL","DeltaCt"]
  Platesum[i,"DeltaDeltaCtsd"] <- sqrt((Platesum[i,"DeltaCtsd"])^2 + (Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Primer == Platesum$Primer[i] & Platesum$Strain == "CTRL","DeltaCtsd"])^2)
} #Just errors on the NAs, which were no template controls or empty

#Calculating the fold change 2^-deltadeltaCt and error (plus and minus, then 2^-)
for(i in 1:nrow(Platesum)) {
  Platesum[i,"FoldChange"] <- 2^-Platesum[i,"DeltaDeltaCt"]
  Platesum[i,"Lower"] <- Platesum[i,"DeltaDeltaCt"] + Platesum[i,"DeltaDeltaCtsd"]
  Platesum[i,"Upper"] <- Platesum[i,"DeltaDeltaCt"] - Platesum[i,"DeltaDeltaCtsd"]
  Platesum[i,"UpperFC"] <- 2^-Platesum[i,"Upper"]
  Platesum[i,"LowerFC"] <- 2^-Platesum[i,"Lower"]
}

#Getting rid of the NA strains (empty columns)
Platesum <- Platesum[!is.na(Platesum$Strain),]
Platesum <- Platesum[!is.na(Platesum$Background),]

#Factoring the Strain names so that they will be in the order I would like
Platesum$Strain <- factor(Platesum$Strain, levels = c("CTRL","C25","C26","C1","C18","C20","C9"))
# Plotting
Labs <- c("AMS5192 OE1","AMS5192 OE2", "L26 OE1", "P75063 OE1","P75063 OE2", "P75016 OE1")
ggplot(data = Platesum[Platesum$Primer == "SKN1" & Platesum$Strain != "CTRL",], aes(x = Strain, y = FoldChange)) +
  geom_col() +
  geom_linerange(data = Platesum[Platesum$Primer == "SKN1"& Platesum$Strain != "CTRL",], ymax = Platesum[Platesum$Primer == "SKN1"& Platesum$Strain != "CTRL",]$UpperFC, ymin = Platesum[Platesum$Primer == "SKN1"& Platesum$Strain != "CTRL",]$LowerFC) +
  theme_bw() +
  scale_x_discrete(labels = Labs) +
  ylim(0, 35) +
  ylab("Fold Change in SKN1\n Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(paste0(fig.dir,"AllSKN1RT.tiff"), plot = last_plot(), width = 2.5, height = 3)

#########################
# C3_03370C
#Clearing environment
rm(list = ls())
data.dir <- "Data/"
fig.dir <- "Figures/"

# Load data
Platedata <- read.csv(paste0(data.dir,"20260313run -  Quantification Cq Results_0.csv"), sep = ",", header = 1)
Platedata <- Platedata[,c("Well","Cq")]
#Read in sample information
Platekey <- read.table(paste0(data.dir,"Key_2026_03_13.txt"), sep = "\t", header = 1)

Platedata <- left_join(Platedata, Platekey, by = "Well")
#There is one crazy replicate in 5192, so I am removing that sample.
Platedata <- Platedata[Platedata$Strain != "C25",]
#Calculating average and standard deviation for each strain/primer pair
Platesum <- Platedata %>%
  group_by(Background, Strain, Primer) %>%
  summarize(Mean = mean(Cq), StDev = sd(Cq))

#Calculating deltaCt for each primer relative to ACT1 (Gene - ACT1) and error (sqrt(summed squares))
for(i in 1:nrow(Platesum)) {
  Platesum[i,"DeltaCt"] <- Platesum[i,"Mean"] - Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Strain == Platesum$Strain[i] & Platesum$Primer == "ACT1","Mean"]
  Platesum[i,"DeltaCtsd"] <- sqrt((Platesum[i,"StDev"])^2 + (Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Strain == Platesum$Strain[i] & Platesum$Primer == "ACT1","StDev"])^2)
}

#Calculating deltadeltaCt for each Strain relative to the WT (Strain - WT) and error(sqrt(summed squares))
for(i in 1:nrow(Platesum)) {
  Platesum[i,"DeltaDeltaCt"] <- Platesum[i,"DeltaCt"] - Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Primer == Platesum$Primer[i] & Platesum$Strain == "CTRL","DeltaCt"]
  Platesum[i,"DeltaDeltaCtsd"] <- sqrt((Platesum[i,"DeltaCtsd"])^2 + (Platesum[Platesum$Background == Platesum$Background[i] & Platesum$Primer == Platesum$Primer[i] & Platesum$Strain == "CTRL","DeltaCtsd"])^2)
} #Just errors on the NAs, which were no template controls or empty

#Calculating the fold change 2^-deltadeltaCt and error (plus and minus, then 2^-)
for(i in 1:nrow(Platesum)) {
  Platesum[i,"FoldChange"] <- 2^-Platesum[i,"DeltaDeltaCt"]
  Platesum[i,"Lower"] <- Platesum[i,"DeltaDeltaCt"] + Platesum[i,"DeltaDeltaCtsd"]
  Platesum[i,"Upper"] <- Platesum[i,"DeltaDeltaCt"] - Platesum[i,"DeltaDeltaCtsd"]
  Platesum[i,"UpperFC"] <- 2^-Platesum[i,"Upper"]
  Platesum[i,"LowerFC"] <- 2^-Platesum[i,"Lower"]
}

#Removing the empty well spots
Platesum <- Platesum[Platesum$Strain != "",]

#Factoring the Strain names so that they will be in the order I would like
Platesum$Strain <- factor(Platesum$Strain, levels = c("CTRL","C27","C1","C5","C18","C9"))
# Plotting
Labs <- c("AMS5192 OE1", "L26 OE1","L26 OE2", "P75063 OE1","P75016 OE1")
ggplot(data = Platesum[Platesum$Primer == "C3" & Platesum$Strain != "CTRL",], aes(x = Strain, y = FoldChange)) +
  geom_col() +
  geom_linerange(data = Platesum[Platesum$Primer == "C3"& Platesum$Strain != "CTRL",], ymax = Platesum[Platesum$Primer == "C3"& Platesum$Strain != "CTRL",]$UpperFC, ymin = Platesum[Platesum$Primer == "C3"& Platesum$Strain != "CTRL",]$LowerFC) +
  theme_bw() +
  scale_x_discrete(labels = Labs) +
  ylim(0, 15) +
  ylab("Fold Change in C3_03370C\n Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(paste0(fig.dir,"AllC3RT.tiff"), plot = last_plot(), width = 2.5, height = 3)
