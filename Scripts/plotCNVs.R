############
# Code for generating figures of read depth plots for all CNV strains
############
# Contents
# 1. Setting up the environment --> line 9
# 2. Calculating rolling median of rolling means --> line 29
# 3. Figure generation --> line 130

####################
# Environment set-up
# clear environment
rm(list = ls())
# set working directory
setwd('/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/')
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
#Load required libraries
library(RcppRoll)
library(ggplot2)
library(ape)
library(openxlsx)
#Set data directory containing all needed depth files, generated in CNVstraindepth.sh
data.dir <- "Data/"
#Set directory for Figure output
fig.dir <- "Figures/Updated/"
#Read in file names
depth_filenames <- list.files(path = data.dir, pattern = "_depth.txt", full.names = TRUE)

##########################
# Set up function for calculating rolling read depth
n <- 500
#Function to calculate relative read depth
RELATIVEDEPTH <- function(x, n) {
  if(missing(n) == TRUE) {
    n <- 500
  }
  x1 <- read.delim(x, header=FALSE)
  x_genome<-x1[x1$V1 != "Ca19-mtDNA",]
  x_mediandepth<-median(x_genome$V3) 
  x_rollmean <- roll_mean(x_genome$V3,n)[seq_len(length(x_genome$V3)-n+1)%%n==1]
  x_relativedepth <- x_rollmean/x_mediandepth
  x_relativedepth <- data.frame(x_relativedepth)
  x_relativedepth$Copy_Number <- x_relativedepth$x_relativedepth*2 #For a diploid organism
  position <- c(x_genome$V2[seq_len(length(x_genome$V3)-n+1)%%n==1])
  x_relativedepth$position <- position[1:(length(position) -1)]
  x_relativedepth$pos_Mb <- x_relativedepth$position/1000000
  chromosome <- c(x_genome$V1[seq_len(length(x_genome$V3)-n+1)%%n==1])
  x_relativedepth$Chromosome <- chromosome[1:(length(chromosome)-1)]
  return(x_relativedepth)
}
#Apply function across list of filenames - this may take a little time.
Relativedepth <- lapply(depth_filenames, RELATIVEDEPTH)
nopath <- list.files(path = data.dir, pattern = "_depth.txt", full.names = FALSE)
names(Relativedepth) <- nopath

# Now specifying the function for calculating a rolling median of the rolling means for further smoothing
ROLLMEDOVERLAP <- function(x, n) {
  if(missing(n) == TRUE) {
    n <- 25
  }
  x_sorted <- x
  x_split_1 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr1_C_albicans_SC5314"),]
  x_split_2 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr2_C_albicans_SC5314"),]
  x_split_3 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr3_C_albicans_SC5314"),]
  x_split_4 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr4_C_albicans_SC5314"),]
  x_split_5 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr5_C_albicans_SC5314"),]
  x_split_6 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr6_C_albicans_SC5314"),]
  x_split_7 <- x_sorted[which(x_sorted$Chromosome == "Ca21chr7_C_albicans_SC5314"),]
  x_split_R <- x_sorted[which(x_sorted$Chromosome == "Ca21chrR_C_albicans_SC5314"),]
  x_rollmedian_1 <- roll_median(x_split_1$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_2 <- roll_median(x_split_2$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_3 <- roll_median(x_split_3$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_4 <- roll_median(x_split_4$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_5 <- roll_median(x_split_5$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_6 <- roll_median(x_split_6$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_7 <- roll_median(x_split_7$Copy_Number,n, na.rm = TRUE)
  x_rollmedian_R <- roll_median(x_split_R$Copy_Number,n, na.rm = TRUE)
  x_split_1 <- cbind(x_split_1[1:(nrow(x_split_1) - 24),], x_rollmedian_1)
  x_split_2 <- cbind(x_split_2[1:(nrow(x_split_2) - 24),], x_rollmedian_2)
  x_split_3 <- cbind(x_split_3[1:(nrow(x_split_3) - 24),], x_rollmedian_3)
  x_split_4 <- cbind(x_split_4[1:(nrow(x_split_4) - 24),], x_rollmedian_4)
  x_split_5 <- cbind(x_split_5[1:(nrow(x_split_5) - 24),], x_rollmedian_5)
  x_split_6 <- cbind(x_split_6[1:(nrow(x_split_6) - 24),], x_rollmedian_6)
  x_split_7 <- cbind(x_split_7[1:(nrow(x_split_7) - 24),], x_rollmedian_7)
  x_split_R <- cbind(x_split_R[1:(nrow(x_split_R) - 24),], x_rollmedian_R)
  colnames(x_split_1)[6] <- "rollmedian"
  colnames(x_split_2)[6] <- "rollmedian"
  colnames(x_split_3)[6] <- "rollmedian"
  colnames(x_split_4)[6] <- "rollmedian"
  colnames(x_split_5)[6] <- "rollmedian"
  colnames(x_split_6)[6] <- "rollmedian"
  colnames(x_split_7)[6] <- "rollmedian"
  colnames(x_split_R)[6] <- "rollmedian"
  x_backtogether <- rbind(x_split_1, x_split_2, x_split_3, x_split_4, x_split_5, x_split_6, x_split_7, x_split_R)
  return(x_backtogether)
}

Rollmedlist <- lapply(Relativedepth, ROLLMEDOVERLAP)
names(Rollmedlist) <- names(Relativedepth)

#Also reading in files for the positions of all genes in the CNV regions, the genes targeted by our CRISPRa guides, and positions of long inverted repeats
Allfeatures22 <- read.gff(paste0(data.dir,"C_albicans_SC5314_version_A22-s07-m01-r213_features.gff"))
Allgenes22 <- Allfeatures22[Allfeatures22$type == "gene",]
getAttributeField <- function (x, field, attrsep = ";") { #From the 'ballgown' package, which appears to be depricated
  s = strsplit(x, split = attrsep, fixed = TRUE)          
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}
Allgenes22$SysName <- getAttributeField(Allgenes22$attributes, field = "ID") 
#This is in a different format than the names in the guide list, so I will edit them quickly
Allgenes22 <- Allgenes22[grep("_A",Allgenes22$SysName),]
Allgenes22$SysName <- gsub("_","", Allgenes22$SysName)
Allgenes22$Copy <- rep(0, nrow(Allgenes22)) #To have them listed at the bottom of the plots
# Reading in the bed file for all genes targeted by guides
Guidebed <- read.table(paste0(data.dir,"Chr134guides.bed"), sep = "\t")
Guidebed$SysName <- sapply(strsplit(Guidebed$V4, split = "_", fixed = TRUE), `[`, 1)
Allgenes22$Targeted <- ifelse(Allgenes22$SysName %in% Guidebed$SysName, "TARG", "NO")

#Also adding in the positions of the long inverted repeats from Todd et al, 2019
LIRS <- read.xlsx(paste0(data.dir,"elife-45954-supp2-v2.xlsx"))

###########################
# Figure generation
# Original CNV strains with positions of long inverted repeats and ticks for genes, colored by whether targeted with CRISPRa

#AMS4105 chromosome 1 section
ggplot(data = Rollmedlist[["AMS4105_depth.txt"]][Rollmedlist[["AMS4105_depth.txt"]]$Chromosome == "Ca21chr1_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  geom_segment(data = Allgenes22[Allgenes22$seqid == "Ca22chr1A_C_albicans_SC5314",], aes(x = start/1000000, xend = end/1000000, y = Copy, color = Targeted), linewidth = 2) +
  scale_color_manual(values = c("#064259","#C50043")) +
  geom_vline(xintercept = LIRS[LIRS$Chromosome.Sequence.1 == "1" & LIRS$Chromosome.Sequence.2 == "1" & LIRS$Start.Sequence.1 >= 2000000 & LIRS$Start.Sequence.2 >= 2000000,"Start.Sequence.1"]/1000000, linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,7) +
  xlim(2,3.177) +
  xlab("Position on Chromosome 1 (Mb)") +
  ylab("Copy Number \n(Relative Sequencing Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS4105Chr1black.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)

#AMS3092 with Chr3L CNV
ggplot(data = Rollmedlist[["AMS3092_depth.txt"]][Rollmedlist[["AMS3092_depth.txt"]]$Chromosome == "Ca21chr3_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  geom_segment(data = Allgenes22[Allgenes22$seqid == "Ca22chr3A_C_albicans_SC5314",], aes(x = start/1000000, xend = end/1000000, y = Copy, color = Targeted), linewidth = 2) +
  scale_color_manual(values = c("#064259","#C50043")) +
  geom_vline(xintercept = LIRS[LIRS$Chromosome.Sequence.1 == "3" & LIRS$Chromosome.Sequence.2 == "3" & LIRS$Start.Sequence.1 <= 823240 & LIRS$Start.Sequence.1 >= 600000 & LIRS$Start.Sequence.2 <= 823240 & LIRS$Start.Sequence.2 >= 600000,"Start.Sequence.1"]/1000000, linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,20) +
  xlim(0.5,0.823240) +
  xlab("Chromosome 3L (Mb)") +
  ylab("Copy Number \n(Relative Sequencing Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS3092Chr3black.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 4)

#AMS4397 with Chr3L CNV
ggplot(data = Rollmedlist[["AMS4397_depth.txt"]][Rollmedlist[["AMS4397_depth.txt"]]$Chromosome == "Ca21chr3_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  geom_segment(data = Allgenes22[Allgenes22$seqid == "Ca22chr3A_C_albicans_SC5314",], aes(x = start/1000000, xend = end/1000000, y = Copy, color = Targeted), linewidth = 2) +
  scale_color_manual(values = c("#064259","#C50043")) +
  geom_vline(xintercept = LIRS[LIRS$Chromosome.Sequence.1 == "3" & LIRS$Chromosome.Sequence.2 == "3" & LIRS$Start.Sequence.1 <= 1500000 & LIRS$Start.Sequence.1 >= 1000000 & LIRS$Start.Sequence.2 <= 1500000 & LIRS$Start.Sequence.2 >= 1000000,"Start.Sequence.1"]/1000000, linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,12) +
  xlim(0.823240,1.6) +
  xlab("Chromosome 3R (Mb)") +
  ylab("Copy Number \n(Relative Sequencing Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS4379Chr3black.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 4)

#AMS5778 with Chr4 CNV
ggplot(data = Rollmedlist[["AMS5778_depth.txt"]][Rollmedlist[["AMS5778_depth.txt"]]$Chromosome == "Ca21chr4_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  geom_segment(data = Allgenes22[Allgenes22$seqid == "Ca22chr4A_C_albicans_SC5314",], aes(x = start/1000000, xend = end/1000000, y = Copy, color = Targeted), linewidth = 2) +
  scale_color_manual(values = c("#064259","#C50043")) +
  geom_vline(xintercept = LIRS[LIRS$Chromosome.Sequence.1 == "4" & LIRS$Chromosome.Sequence.2 == "4" & LIRS$Start.Sequence.1 <= 900000 & LIRS$Start.Sequence.1 >= 500000 & LIRS$Start.Sequence.2 <= 900000 & LIRS$Start.Sequence.2 >= 500000,"Start.Sequence.1"]/1000000, linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,7) +
  xlim(0.3,1.0) +
  xlab("Chromosome 4L (Mb)") +
  ylab("Copy Number \n(Relative Sequencing Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS5778Chr4black.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)

#AMS4702 with Chr4 CNV
ggplot(data = Rollmedlist[["AMS4702_depth.txt"]][Rollmedlist[["AMS4702_depth.txt"]]$Chromosome == "Ca21chr4_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  geom_segment(data = Allgenes22[Allgenes22$seqid == "Ca22chr4A_C_albicans_SC5314",], aes(x = start/1000000, xend = end/1000000, y = Copy, color = Targeted), linewidth = 2) +
  scale_color_manual(values = c("#064259","#C50043")) +
  geom_vline(xintercept = LIRS[LIRS$Chromosome.Sequence.1 == "4" & LIRS$Chromosome.Sequence.2 == "4" & LIRS$Start.Sequence.1 <= 1000000 & LIRS$Start.Sequence.1 >= 500000 & LIRS$Start.Sequence.2 <= 1000000 & LIRS$Start.Sequence.2 >= 500000,"Start.Sequence.1"]/1000000, linetype = "dashed") +
  geom_vline(xintercept = 0.85, linetype = "dashed") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,15) +
  xlim(0.3,1.0) +
  xlab("Chromosome 4L (Mb)") +
  ylab("Copy Number \n(Relative Sequencing Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS4702Chr4black.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)

#Getting the number of genes in each CNV region.
nrow(Allgenes22[Allgenes22$seqid == "Ca22chr1A_C_albicans_SC5314" & Allgenes22$Targeted == "TARG",])
nrow(Allgenes22[Allgenes22$seqid == "Ca22chr3A_C_albicans_SC5314" & Allgenes22$start <= 800000 & Allgenes22$Targeted == "TARG",])
nrow(Allgenes22[Allgenes22$seqid == "Ca22chr3A_C_albicans_SC5314" & Allgenes22$start >= 800000 & Allgenes22$Targeted == "TARG",])
nrow(Allgenes22[Allgenes22$seqid == "Ca22chr4A_C_albicans_SC5314" & Allgenes22$start <= 703000 & Allgenes22$Targeted == "TARG",])
nrow(Allgenes22[Allgenes22$seqid == "Ca22chr4A_C_albicans_SC5314" & Allgenes22$start <= 673000 & Allgenes22$Targeted == "TARG",])


# CNV regions with specific genes of interest highlighted
# AMS7083 with a few key genes highlighted
ggplot(data = Rollmedlist[["AMS7083_depth.txt"]][Rollmedlist[["AMS7083_depth.txt"]]$Chromosome == "Ca21chr1_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,8) +
  xlim(2,3.177) +
  geom_vline(xintercept = 2869627/1000000, color = "purple4") + #TYE7
  geom_vline(xintercept = 2568085/1000000, color = "forestgreen") + #C1_11720W
  xlab("Chromosome 1 (Mb)") +
  ylab("Copy Number \n(Relative WGS Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS7083Chr1TYE7Genes.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)

#AMS7084 Chr3 CNV with CDR1, SKN1, and YCK2 highlighted.
ggplot(data = Rollmedlist[["AMS7084_depth.txt"]][Rollmedlist[["AMS7084_depth.txt"]]$Chromosome == "Ca21chr3_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,8) +
  xlim(1,1.5) +
  geom_vline(xintercept = 1146022/1000000, color = "purple4") + #CDR1
  geom_vline(xintercept = 1302477/1000000, color = "forestgreen") + #SKN1
  geom_vline(xintercept = 1263719/1000000, color = "forestgreen") + #YCK2
  xlab("Chromosome 3 (Mb)") +
  ylab("Copy Number \n(Relative WGS Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS7084Chr3Genes.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)

#AMS5778 with NCP1 highlighted
ggplot(data = Rollmedlist[["AMS5778_depth.txt"]][Rollmedlist[["AMS5778_depth.txt"]]$Chromosome == "Ca21chr4_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylim(0,7) +
  xlim(0.45,0.80) +  
  geom_vline(xintercept = 668228/1000000, color = "purple4") + #NCP1
  xlab("Chromosome 4L (Mb)") +
  ylab("Copy Number \n(Relative Sequencing Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS5778Chr4Genes.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)

#AMS4397 with specific genes of interest highlighted
ggplot(data = Rollmedlist[["AMS4397_depth.txt"]][Rollmedlist[["AMS4397_depth.txt"]]$Chromosome == "Ca21chr3_C_albicans_SC5314",], aes(x=pos_Mb, y=Copy_Number)) +
  #geom_point(alpha=0.5, size=1, color = "#4C4281") + 
  geom_point(aes(x = pos_Mb, y = rollmedian), color = "black", alpha = 0.5) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  geom_vline(xintercept = 1146022/1000000, color = "purple4") + #CDR1
  geom_vline(xintercept = 1302477/1000000, color = "forestgreen") + #SKN1
  geom_vline(xintercept = 1263719/1000000, color = "forestgreen") + #YCK2
  ylim(0,12) +
  xlim(0.823240,1.60) +
  xlab("Chromosome 3R (Mb)") +
  ylab("Copy Number \n(Relative WGS Depth)") +
  theme(legend.position = "none", axis.text = element_text(size = 15), axis.title = element_text(size = 20))
ggsave("AMS4397Chr3Genes.tiff", path = fig.dir, plot = last_plot(), height = 5, width = 6)


