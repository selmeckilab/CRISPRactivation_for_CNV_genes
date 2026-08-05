#######################
# Analysis of growth curve data for clean CNV strains. Output (CNVGC.txt) will be used for comparison
# to additive single gene overexpression effects in 'AnalysisandFigs.R'
#######################

# Setting up the environment
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
#Load libraries
library(tidyr)
library(ggplot2)
library(growthcurver)
library(dplyr)
library(cowplot)
library(openxlsx)
#Make sure there is nothing in the environment
rm(list = ls())

#Read in data files (Plate 3 was a plate of additional MCF concentrations. They were not used in the final analysis, as 0.016ug/mL was the most informative concentration and was already included in Plates 1/2)
data.dir <- "Data"
Plate1 <- read.xlsx(paste0(data.dir,"/2026_07_02_GC_PLATE1.xlsx"), sheet = 1, rows = c(40:233), cols = c(4:100)) 
Plate2 <- read.xlsx(paste0(data.dir,"/2026_07_02_GC_PLATE2.xlsx"), sheet = 1, rows = c(40:233), cols = c(4:100)) 
Plate4 <- read.xlsx(paste0(data.dir,"/2026_07_02_GC_37C.xlsx"), sheet = 1, rows = c(29:222), cols = c(4:100)) 
Platekey <- read.xlsx(paste0(data.dir,"/2026_07_02_Platekey.xlsx")) 

#Processing and blanking each plate separately and then joining them together. 
#Plate1
x <- 0
for (i in 1:nrow(Plate1)) {
  Plate1[i,"Time"] <- x
  x <- x + 0.25
}

#A few quick looks at the data to make sure it all read in fine. 
matplot(Plate1[,1:96],type = "l")

Plate1melt <- pivot_longer(Plate1, cols = colnames(Plate1)[1:96], names_to = "Well", values_to = "OD")
ggplot(data = Plate1melt, aes(x = Time, y = OD)) +
  geom_point() +
  facet_wrap(~ Well, ncol = 12)

#Modeling growth curves with growthcurveR
gc_outPlate1 <- SummarizeGrowthByPlate(Plate1, plot_fit = FALSE)
#Checking a few
gc_fitPlate1 <- SummarizeGrowth(Plate1$Time, Plate1$D1) # Not bad even for the weird drug ones. 
plot(gc_fitPlate1)

#Adding on sample info for the growthcurveR stats
colnames(gc_outPlate1)[1] <- "Well"
gc_out_infoPlate1 <- left_join(gc_outPlate1, Platekey[Platekey$Plate == "Plate1",], by = "Well")

#Plate2
x <- 0
for (i in 1:nrow(Plate2)) {
  Plate2[i,"Time"] <- x
  x <- x + 0.25
}

#A few quick looks at the data to make sure it all read in fine. 
matplot(Plate2[,1:96],type = "l") #Something a little odd here where cells get an extra 'bump' about half way through.

Plate2melt <- pivot_longer(Plate2, cols = colnames(Plate2)[1:96], names_to = "Well", values_to = "OD")
ggplot(data = Plate2melt, aes(x = Time, y = OD)) + #It is the bottom row that has that bump - pH8, I think it might have dried out a bit. 
  geom_point() +
  facet_wrap(~ Well, ncol = 12)

#Modeling growth curves with growthcurveR
gc_outPlate2 <- SummarizeGrowthByPlate(Plate2, plot_fit = FALSE)
#Checking a few
gc_fitPlate2 <- SummarizeGrowth(Plate2$Time, Plate2$H1) #Get pretty well smoothed out in the model 
plot(gc_fitPlate2)

#Adding on sample info for the growthcurveR stats
colnames(gc_outPlate2)[1] <- "Well"
gc_out_infoPlate2 <- left_join(gc_outPlate2, Platekey[Platekey$Plate == "Plate2",], by = "Well")

#Plate4
x <- 0
for (i in 1:nrow(Plate4)) {
  Plate4[i,"Time"] <- x
  x <- x + 0.25
}

#A few quick looks at the data to make sure it all read in fine. 
matplot(Plate4[,1:96],type = "l")

Plate4melt <- pivot_longer(Plate4, cols = colnames(Plate4)[1:96], names_to = "Well", values_to = "OD")
ggplot(data = Plate4melt, aes(x = Time, y = OD)) +
  geom_point() +
  facet_wrap(~ Well, ncol = 12)

#Modeling growth curves with growthcurveR
gc_outPlate4 <- SummarizeGrowthByPlate(Plate4, plot_fit = FALSE)
#Checking a few
gc_fitPlate4 <- SummarizeGrowth(Plate4$Time, Plate4$B1)
plot(gc_fitPlate4)

#Adding on sample info for the growthcurveR stats
colnames(gc_outPlate4)[1] <- "Well"
gc_out_infoPlate4 <- left_join(gc_outPlate4, Platekey[Platekey$Plate == "Plate4",], by = "Well")
#And now getting rid of the empty wells in this plate
gc_out_infoPlate4 <- gc_out_infoPlate4[gc_out_infoPlate4$Well %in% Platekey[Platekey$Plate == "Plate4",]$Well,]

###############################
#Make everything relative to YPAD by subtracting
#Combining all other plate dataframes to look at all conditions 
gc_out_info <- rbind(gc_out_infoPlate1, gc_out_infoPlate2, gc_out_infoPlate4)
#Getting rid of the blanks
gc_out_info <- gc_out_info[gc_out_info$Strain != "blank",]
for (i in 1:nrow(gc_out_info)) {
  gc_out_info[i,"YPADREL"] <- gc_out_info[i,"auc_e"] - gc_out_info[gc_out_info$Rep == gc_out_info$Rep[i] & gc_out_info$Strain == gc_out_info$Strain[i] & gc_out_info$Condition == "YPAD","auc_e"]
}

#And now making each strain relative to the mean progenitor values
for (i in 1:nrow(gc_out_info)) {
  gc_out_info[i,"PROGREL"] <- gc_out_info[i,"YPADREL"] - mean(gc_out_info[gc_out_info$Condition == gc_out_info$Condition[i] & gc_out_info$Strain == gc_out_info$Progenitor[i],"YPADREL"])
}

#Keeping replicates relative to each other
for (i in 1:nrow(gc_out_info)) {
  gc_out_info[i,"PROGRELREP"] <- gc_out_info[i,"YPADREL"] - gc_out_info[gc_out_info$Condition == gc_out_info$Condition[i] & gc_out_info$Replicate == gc_out_info$Replicate[i] & gc_out_info$Strain == gc_out_info$Progenitor[i],"YPADREL"]
}

#Reading out the processed data, to be used in "AnalysisandFigs.R"
output.dir <- "Output/"
write.table(gc_out_info, paste0(output.dir,"CNVGC.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
