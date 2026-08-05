###################
# Code to analyze a pilot experiment to compare our two P75063 pools growth in YPAD.
##################
# Contents:
# 1. Environment set-up
# 2. 

##################
# Clear environment, load libraries, and specify file paths
rm(list = ls())
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
library(tidyr)
library(dplyr)
library(ggplot2)
library(openxlsx)

# location of pipeline functions
source("Scripts/ForGithub/pipeline_functions_edited.r") 
# Location of count data for all samples
data.dir <- "Data/"
# Location where figures will be stored
fig.dir <- "Figures/Updated/"
# Location where processed data will be read out
output.dir <- "Output/"

#########################
# Begin data input

# Reading in the raw data
filenames <- list.files(path = paste0(data.dir,"countfiles/pilot/"), pattern = "counts.csv", full.names = FALSE)
# Function for grabbing all count files from each sample and creating a list
GETCOUNTS <- function(x) {
  countdf <- read.csv(paste0(data.dir,"countfiles/pilot/",x), header = 1)
  return(countdf)
}
# Applying the function over the list of files
countlist <- lapply(filenames, GETCOUNTS)
names(countlist) <- filenames

# Generating a count matrix from the list of count files
for (i in 1:length(countlist)) {
  countlist[[i]]$Sample <- c(rep(names(countlist[i]), nrow(countlist[[i]])))
  countlist[[i]]$Sample <- gsub("counts.csv","",countlist[[i]]$Sample)
}
Alltog <- do.call("rbind",countlist)

#Now removing the 8nt anchors 
Alltog$Barcode <- substring(Alltog$Barcode, 9, 28)
Allrawcounts <- pivot_wider(Alltog, names_from = Sample, values_from = Count, id_cols = Barcode)

write.table(Allrawcounts, "Data/Pilotcounts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# List of neutral barcodes
ref=c("CTGCAGAGAGATCATAACTT","TAACATTGGAGTCGCTACTC","AAGGGAGTCGTAATAGGAAG","GCATGACTGTGGATAGAGGT","TCCACTGGACTCGCCCCCTC","TACCACGTCTTCAGATCTCG","ATCAGTCCGAAGCCGTGCGG","GTCACAGCTCTACTACAGCA","TCTCGCGACTGCCACATAAC","TAAAGCTTAGGGACCGGAGT","CCGTGTGCGATCCCTTCTTT","AGACTCCAGCGTTGGCCAAA","AGCTCTACTGTTACAAGGGG","GCTAGCTTCACTAGGGGTAA","TCGGCCCTAGACTCAATATG","GATCTGGCTCGGGACTCGTT","TCCGCTGCGCAAACTTCACT","TGGACCAGGTTCACCGAGTG","TCTTCCCATTCGGACCCAAG","TGTTTAGGCAGTGATCAGCC","GGTATCGTCGTGTCTCATTA","TGTTGGATCGTCCCTAGGAA","TCTGACGATCTGTTGTGAGA","TAATCGGTAACTCGTAATCT","AATACGTTTGCCGTGAAGAT","GGCCGGATAGATAAGGGAGT","TACTCAAGGTGGCTGACTAT","TACCATAAGGCGGAGTTCGT","TCGAGACCAGCTCCTAGTAA","GGTGTATGAGGTTGCGGCCA","GCTGAGTTGTGGCTAGGAGC","GTTTACGGGCATGAAGTGCA","CAGCACAGCTAGGACCCAAC","CGTCAGACTTAGTACGCGTA","GGATTAGCAGGGTAATGCAA","GAAGCTAATCAGGCCAAGAT","GTCCTGATCCTCGCCACAAC","GCAATCTCGCGGGTTGGCTA","TGTTTGGGGCCTTCCCGGAC","ACGAGGCGCCGCGTCGTAAT","TTACCTGCAACGGTCGATTG","AATCTCACGTTCCGCCTGAC","ACCATGACGCGTCTATGTAC","TCCGACGAGGAAGTCTGCAT","GCATATAGGCCTGCCCGCGG","CTTCGCCGAGTAACGTCCGA","ATGAGGGAAGAGGCCGAACG","GGGCCGTAATCTAGACTGGC","TGTCATACGGTCGTTTCACT","ATACTCTGTGTCGCGATGAG","TAGCCTTGATCCGCGTGGGC","ACGATTTTGACCGAGCCCCA","TTATGCCGAAACCTGACATG","GTTGTCAGACTGCCAAGTTA","AGTGTCCCGTCATGCTCCAG","TCTTTCGGGAGCGCCTGATT","TCCCCAAGGCTGTCCCCCAA","GACCTCTTTCGACTAGGCCA","GCACTAGCGTCCCACGAATG","AGACTGAAACGGCAGCGCGA") # here you write the neutral barcode references
# Specify parameters
poolcount_loc = "Pilotcounts.txt" 
pseudocount = 0.1 
T0_thresh = 50.1  # minimum count for a T0 barcode

#####read files
raw_counts = read.table(paste0(data.dir,poolcount_loc), sep = "\t", stringsAsFactors = FALSE, header = 1) 
#To make this more easily compatible with the functions - I should just change this in the input file
colnames(raw_counts)[1] <- 'barcode'
#generating the metadata file from the sample column names directly
metadata <- data.frame(Samples = colnames(raw_counts[2:ncol(raw_counts)]))
metadata$Timepoint <- rep("T1", nrow(metadata))
metadata$Condition <- substring(metadata$Samples, 3, 3) #Warning, this needs to be fixed for the overnights
metadata$Replicate <- substring(metadata$Samples, 4, 4) #Warning, this needs to be fixed for the overnights
T0s = c("R1_S53","R2_S54","R3_S55")
metadata[metadata$Samples %in% T0s, "Timepoint"] <- "T0"
metadata[metadata$Samples %in% T0s, "Condition"] <- "ALL"
metadata[metadata$Samples %in% T0s, "Replicate"] <- substring(metadata[metadata$Samples %in% T0s, "Samples"], 2, 2)

#########
# Data filtering
######## 
counts = raw_counts
# 1. remove barcodes which never have more than 50 reads in the T0 samples
low_barcodes = get_barcodes_with_any_below_threshold(counts, T0s, 51)
print(paste("removed", 
            100 * ((length(low_barcodes) / nrow(counts))),
            " percent of remaining barcodes because the read count is too low in at least one T0"))
counts = counts %>%
  filter(!(barcode %in% low_barcodes))

# 2. Identify samples with less than a median of 50 reads per barcode
low_count_reps = get_replicates_with_low_median_gene_reads(counts, T0_thresh)
print(paste("removing reps: ", low_count_reps, ", due to low median barcode counts [if empty, all passed threshold]"))
counts = counts %>%
  select(-low_count_reps)

# 3. Remove barcodes that are the only barcode targeting a gene. Uses Supplementary File S4 - sgRNA sequences, sheet 2
barconvert <- read.xlsx(paste0(data.dir,"Supplementary Table S4 - sgRNA sequences.xlsx"), sheet = 2, colNames = TRUE)
colnames(barconvert)[3] <- "barcode"
colnames(barconvert)[1] <- "gRNA"
few_barcodes <- names(which(table(barconvert$gRNA) == 1))
few_barcodes_sequence <- barconvert[barconvert$gRNA %in% few_barcodes, "barcode"]
counts = counts %>%
  filter(!(barcode %in% few_barcodes_sequence)) #Removes 26

######################
### the main loop; it uses a different T0 for each strain fitness calculation,
### calculating fitness per-sample, each replicate individually. 

strain_fitness = data.frame()
for (rep in unique(metadata$Replicate)){
  print(rep)
  # subsample this rep's data, with a different T0 per loop
  samples = metadata$Samples[metadata$Replicate == rep]
  T0_sample = metadata$Samples[metadata$Replicate == rep & metadata$Timepoint == "T0"]
  this_rep = counts %>%
    select(barcode, all_of(samples)) %>%
    rename(T0 = T0_sample)
  
  # perform the rest of the fitness calculation for this replicate
  this_rep = add_pseudocount(this_rep, 0.1)
  this_rep = normalize_using_reference_genes(this_rep, ref)
  this_rep = calculate_strain_fitness_from_T0(this_rep) %>% #This says it uses depth normalized counts, but i am not sure about that - normalized by the reference barcodes, not total
    mutate(replicate = ifelse(replicate == "T0", T0_sample, replicate))
  
  strain_fitness = rbind(strain_fitness, this_rep)
  rm(this_rep)
}
fitness = strain_fitness %>%
  filter(!(replicate %in% T0s))


########################
# Adding information about which gene each barcode is targeting and the distances to the gene to make gene-level weighted averages
barcodeinfo <- read.table(paste0(data.dir,"barcodeinfo.txt"), sep = "\t", header = 1)
# there are 91 barcodes that are located far away from their target genes and are also low quality in 
# guide scores. They should be removed. 
farbarcodes <- barcodeinfo[abs(barcodeinfo$DIST) > 1000 & !is.na(barcodeinfo$DIST),"barcode"] 
fitness <- fitness[!fitness$barcode %in% farbarcodes,]
#Joining on the barcode distance information
fitness_id <- left_join(fitness, barcodeinfo, by = 'barcode')
#Tagging the nontargeting guides as such, their distance will be NA
fitness_id[fitness_id$barcode %in% ref,"gRNA"] <- "NonTargeting"

#Calculating weighted gene average fitness
fitness_id = calculate_weighted_gene_fitness_edit(fitness_id, count_cap = 50)

###########
# Processing replicates to get mean and standard deviation values

colnames(fitness_id)[2] <- "Samples"
fitness_id <- left_join(fitness_id, metadata, by = 'Samples')

#After filtering, how many genes total are we left with?
length(unique(fitness_id$gRNA)) #789 for this data as well as the bigbatch experiment

#Replacing 'fitness' with 'strain_fitness' for the NonTargeting controls, because we do not want to average those
fitness_id$fitness <- ifelse(fitness_id$gRNA == "NonTargeting", fitness_id$strain_fitness, fitness_id$fitness)

#Normalize the fitnesses to the fitness of the NonTargeting controls in the same sample - first remove individual Nontargeting barcodes that are suspiciously low fitness
Barpiv <- pivot_wider(fitness_id[fitness_id$gRNA == "NonTargeting",], id_cols = c("barcode"), names_from = Samples, values_from = strain_fitness)
fitness_id[fitness_id$gRNA == "NonTargeting","fitness"] <- ifelse(fitness_id[fitness_id$gRNA == "NonTargeting",]$strain_fitness_var >= 0.04, NA, fitness_id[fitness_id$gRNA == "NonTargeting",]$fitness)

NonTargetFit <- fitness_id %>%
  filter(gRNA == "NonTargeting") %>%
  group_by(Samples) %>%
  summarize(MedianFit = median(fitness, na.rm = TRUE), MeanFit = mean(fitness, na.rm = TRUE))

fitness_id <- left_join(fitness_id, NonTargetFit[,c("Samples","MedianFit")], by = "Samples")
fitness_id$AdjFitness <- fitness_id$fitness - fitness_id$MedianFit

#Eliminating redundency of gene-wise fitness estimates by taking one line per gene
gene_only <- fitness_id[,c(c("Samples","gRNA","genestart","geneend","genestrand","fitness","fitness_left","fitness_right","Timepoint","Condition","Replicate","AdjFitness"))]
gene_only <- unique(gene_only)

gene_only <- gene_only %>%
  group_by(Condition, Timepoint, gRNA) %>%
  mutate(Mean = mean(AdjFitness), SE = (sd(AdjFitness))/3)

#Removing Non-Targeting controls, since those have different Means with the same gRNA
gene_only <- gene_only[gene_only$gRNA != "NonTargeting",]

#Removing the now-redundant replicate measures
gene_only_norep <- gene_only[gene_only$Replicate == 2,]

#Calculating means and sd values
CompmultirepsMean <- pivot_wider(gene_only_norep, id_cols = c('gRNA'), names_from = c('Condition','Timepoint'), values_from = 'Mean')
CompmultirepsSD <- pivot_wider(gene_only_norep, id_cols = c('gRNA'), names_from = c('Condition','Timepoint'), values_from = 'SE')
colnames(CompmultirepsMean)[2:length(colnames(CompmultirepsMean))] <- paste0(colnames(CompmultirepsMean)[2:length(colnames(CompmultirepsMean))],"_Mean")
colnames(CompmultirepsSD)[2:length(colnames(CompmultirepsSD))] <- paste0(colnames(CompmultirepsSD)[2:length(colnames(CompmultirepsSD))],"_SD")

CompmultirepsMean <- left_join(CompmultirepsMean, CompmultirepsSD, by = "gRNA")

# Reading out these fitness values to compare to the fitness values for the other P75063 pool in the bulk competition
Pilotfitness <- fitness
Pilotfitness_id <- fitness_id
PilotMean <- CompmultirepsMean

#Now loading in the environment from the full competition experiment analysis (generated in "AnalysisandFigs.R")
load(paste0(output.dir,"/ManuscriptFull20260206.R"))

#Now joining the CompmultirepsMean files
Comppools <- left_join(CompmultirepsMean, PilotMean, by = "gRNA")

# Comparing the fitness of the two P75063 pools

ggplot(gene_only_norep[gene_only_norep$Condition == "YPAD" & gene_only_norep$Timepoint == "T1" & gene_only_norep$gRNA != "NonTargeting" & gene_only_norep$Strain == "63",], aes(x = Mean)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(-2,1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Fitness \n (Weighted Mean Gene Fitness)") +
  ylab("Density") +
  ggtitle("P75063 Pool 1")
ggsave("Batch63FitDist.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3)

ggplot(PilotMean, aes(x = Y_T1_Mean)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(-2,1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Fitness \n (Weighted Mean Gene Fitness)") +
  ylab("Density") +
  ggtitle("P75063 Pool 2")
ggsave("Pilot63FitDist.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3)


ggplot(data = Comppools, aes(y = Y_T1_Mean, x = `63_YPAD_T1_Mean`)) +
  geom_point() +
  geom_errorbar(aes(ymin = `Y_T1_Mean` - `Y_T1_SD`, ymax = `Y_T1_Mean` + `Y_T1_SD`)) +
  geom_errorbarh(aes(xmin = `63_YPAD_T1_Mean` - `63_YPAD_T1_SD`, xmax = `63_YPAD_T1_Mean` + `63_YPAD_T1_SD`)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(-1,0.85) +
  ylim(-1,0.85) +
  ylab("P75063 Pool 2 Fitness") +
  xlab("P75063 Pool 1 Fitness") +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("Pool63Comp.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3)

