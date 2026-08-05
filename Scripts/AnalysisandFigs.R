###################
# Code for analysis of bar-seq competition data and generation of main text figures.
# The bar-seq analysis is adapted from code first published in Martinson et al, 2023.
##################

# Contents:
# 1. Environment set-up and data input --> line 13
# 2. Data filtering --> line 101
# 3. Fitness estimates --> line 128
# 4. Comparison of individual barcodes and weighted averages --> line 193
# 5. Analysis and figure generation --> line 254

##################
# Clear environment, load libraries, and specify file paths
rm(list = ls())
setwd("/projects/standard/selmecki/pvzande/CRISPRa/CRISPRaPaper/")
.libPaths('/projects/standard/selmecki/pvzande/MyRlibs_pvzande/')
# Loading libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(UpSetR)
library(openxlsx)
library(RColorBrewer)
library(limma)
library(MKmisc)
library(cowplot)
library(VennDiagram)
library(purrr)

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
filenames <- list.files(path = paste0(data.dir,"countfiles/"), pattern = "counts.csv", full.names = FALSE)
# Function for grabbing all count files from each sample and creating a list
GETCOUNTS <- function(x) {
  countdf <- read.csv(paste0(data.dir,"countfiles/",x), header = 1)
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

# Removing the 8nt anchors from the sgRNA sequences to associate them with the genes they target
Alltog$Barcode <- substring(Alltog$Barcode, 9, 28)
Allrawcounts <- pivot_wider(Alltog, names_from = Sample, values_from = Count, id_cols = Barcode)

# Storing this raw count matrix - this is supplementary file S6
write.table(Allrawcounts, "Data/countfiles/Allrawcounts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# List of neutral barcodes
ref=c("CTGCAGAGAGATCATAACTT","TAACATTGGAGTCGCTACTC","AAGGGAGTCGTAATAGGAAG","GCATGACTGTGGATAGAGGT","TCCACTGGACTCGCCCCCTC","TACCACGTCTTCAGATCTCG","ATCAGTCCGAAGCCGTGCGG","GTCACAGCTCTACTACAGCA","TCTCGCGACTGCCACATAAC","TAAAGCTTAGGGACCGGAGT","CCGTGTGCGATCCCTTCTTT","AGACTCCAGCGTTGGCCAAA","AGCTCTACTGTTACAAGGGG","GCTAGCTTCACTAGGGGTAA","TCGGCCCTAGACTCAATATG","GATCTGGCTCGGGACTCGTT","TCCGCTGCGCAAACTTCACT","TGGACCAGGTTCACCGAGTG","TCTTCCCATTCGGACCCAAG","TGTTTAGGCAGTGATCAGCC","GGTATCGTCGTGTCTCATTA","TGTTGGATCGTCCCTAGGAA","TCTGACGATCTGTTGTGAGA","TAATCGGTAACTCGTAATCT","AATACGTTTGCCGTGAAGAT","GGCCGGATAGATAAGGGAGT","TACTCAAGGTGGCTGACTAT","TACCATAAGGCGGAGTTCGT","TCGAGACCAGCTCCTAGTAA","GGTGTATGAGGTTGCGGCCA","GCTGAGTTGTGGCTAGGAGC","GTTTACGGGCATGAAGTGCA","CAGCACAGCTAGGACCCAAC","CGTCAGACTTAGTACGCGTA","GGATTAGCAGGGTAATGCAA","GAAGCTAATCAGGCCAAGAT","GTCCTGATCCTCGCCACAAC","GCAATCTCGCGGGTTGGCTA","TGTTTGGGGCCTTCCCGGAC","ACGAGGCGCCGCGTCGTAAT","TTACCTGCAACGGTCGATTG","AATCTCACGTTCCGCCTGAC","ACCATGACGCGTCTATGTAC","TCCGACGAGGAAGTCTGCAT","GCATATAGGCCTGCCCGCGG","CTTCGCCGAGTAACGTCCGA","ATGAGGGAAGAGGCCGAACG","GGGCCGTAATCTAGACTGGC","TGTCATACGGTCGTTTCACT","ATACTCTGTGTCGCGATGAG","TAGCCTTGATCCGCGTGGGC","ACGATTTTGACCGAGCCCCA","TTATGCCGAAACCTGACATG","GTTGTCAGACTGCCAAGTTA","AGTGTCCCGTCATGCTCCAG","TCTTTCGGGAGCGCCTGATT","TCCCCAAGGCTGTCCCCCAA","GACCTCTTTCGACTAGGCCA","GCACTAGCGTCCCACGAATG","AGACTGAAACGGCAGCGCGA") # here you write the neutral barcode references
# Specify parameters
poolcount_loc = "countfiles/Allrawcounts.txt" 
pseudocount = 0.1 
T0_thresh = 50.1  # minimum count of 50 for a T0 barcode

# Re-read in the count matrix
raw_counts = read.table(paste0(data.dir,poolcount_loc), sep = "\t", stringsAsFactors = FALSE, header = 1) 
# Re-naming first column for compatibility with pipeline functions
colnames(raw_counts)[1] <- 'barcode'

# Generating the metadata file from the sample column names directly
metadata <- data.frame(Samples = colnames(raw_counts[2:ncol(raw_counts)]))
metadata$Timepoint <- sapply(strsplit(metadata$Samples, split = "_", fixed = TRUE), `[`, 1)
Initials <- grep("T0",metadata$Samples)
T0s = metadata$Samples[Initials]
for (i in 1:length(Initials)) {
    metadata[i,"Condition"] <- "ALL"
    metadata[i,"Hold"] <- sapply(strsplit(metadata$Samples[i], split = "_", fixed = TRUE), `[`, 2)
    metadata[i,"Strain"] <- sapply(strsplit(metadata$Hold[i], split = ".", fixed = TRUE), `[`, 1)
    metadata[i,"Replicate"] <- sapply(strsplit(metadata$Hold[i], split = ".", fixed = TRUE), `[`, 2)
}
#Careful of this step - Initials must be the first 14 samples
for (i in 15:nrow(metadata)) { #Note, this is hard coded!!! Initials must be the first 18 samples
    metadata$Timepoint[i] <- sapply(strsplit(metadata$Samples[i], split = "_", fixed = TRUE), `[`, 1)
    metadata$Condition[i] <- sapply(strsplit(metadata$Samples[i], split = "_", fixed = TRUE), `[`, 2)
    metadata$Hold[i] <- sapply(strsplit(metadata$Samples[i], split = "_", fixed = TRUE), `[`, 3)
    metadata$Strain[i] <- sapply(strsplit(metadata$Hold[i], split = ".", fixed = TRUE), `[`, 1)
    metadata$Replicate[i] <- sapply(strsplit(metadata$Hold[i], split = ".", fixed = TRUE), `[`, 2)
}

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
### The main loop; it uses a different T0 for each strain fitness calculation,
### calculating fitness per-sample, each replicate individually. 

strain_fitness = data.frame()
for (rep in unique(metadata$Hold)){
    print(rep)
    # subsample this rep's data, with a different T0 per loop
    samples = metadata$Samples[metadata$Hold == rep]
    T0_sample = metadata$Samples[metadata$Hold == rep & metadata$Timepoint == "T0"]
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

# Adjusting fitness according to the fitness effect estimates of the NonTargeting control strains in each sample
colnames(fitness_id)[2] <- "Samples"
fitness_id <- left_join(fitness_id, metadata, by = 'Samples')

# After filtering, how many genes total are we left with?
length(unique(fitness_id$gRNA)) #789

# Replacing 'fitness' with 'strain_fitness' for the NonTargeting controls, because we do not want averages of those
fitness_id$fitness <- ifelse(fitness_id$gRNA == "NonTargeting", fitness_id$strain_fitness, fitness_id$fitness)

# Normalize the fitness to the fitness of the NonTargeting controls in the same sample
# First remove individual Nontargeting barcodes that are suspiciously low fitness and high variance

fitness_id[fitness_id$gRNA == "NonTargeting","fitness"] <- ifelse(fitness_id[fitness_id$gRNA == "NonTargeting",]$strain_fitness_var >= 0.1, NA, fitness_id[fitness_id$gRNA == "NonTargeting",]$fitness)

NonTargetFit <- fitness_id %>%
  filter(gRNA == "NonTargeting") %>%
  group_by(Samples) %>%
  summarize(MedianFit = median(fitness, na.rm = TRUE), MeanFit = mean(fitness, na.rm = TRUE))

fitness_id <- left_join(fitness_id, NonTargetFit[,c("Samples","MedianFit")], by = "Samples")
fitness_id$AdjFitness <- fitness_id$fitness - fitness_id$MedianFit
# Adjusted Fitness is very similar to fitness, as the non-targeting median fitness effects are close to zero, as shown by the following plot.
# plot(fitness_id$AdjFitness ~ fitness_id$fitness)

######################
# Comparing individual barcodes to weighted averages

ggplot(data = fitness_id[fitness_id$Timepoint == "T1" & fitness_id$Condition == "YPAD",], aes(x = fitdiff, y = strain_fitness_var)) +
  geom_point(aes(color = Strain), alpha = 0.5) +
  theme_bw() +
  xlab("Difference in average and individual\nbarcode fitness") +
  ylab("Individual barcode variance") +
  scale_color_manual(values = c("cyan3","forestgreen", "goldenrod3", "mediumvioletred"), labels = c("AMS5192","L26","P75063","P75016"))
ggsave("VarAveCompYPAD.tiff", plot = last_plot(), path = fig.dir, height = 3, width = 4)

ggplot(data = fitness_id[fitness_id$Timepoint == "T1" & fitness_id$Condition == "FLC",], aes(x = fitdiff, y = strain_fitness_var)) +
  geom_point(aes(color = Strain), alpha = 0.5) +
  theme_bw() +
  xlab("Difference in average and individual\nbarcode fitness") +
  ylab("Individual barcode variance") +
  scale_color_manual(values = c("cyan3","forestgreen", "goldenrod3", "mediumvioletred"), labels = c("AMS5192","L26","P75063","P75016"))
ggsave("VarAveCompFLC.tiff", plot = last_plot(), path = fig.dir, height = 3, width = 4)

ggplot(data = fitness_id[fitness_id$Timepoint == "T1" & fitness_id$Condition == "MCF",], aes(x = fitdiff, y = strain_fitness_var)) +
  geom_point(aes(color = Strain), alpha = 0.5) +
  theme_bw() +
  xlab("Difference in average and individual\nbarcode fitness") +
  ylab("Individual barcode variance") +
  scale_color_manual(values = c("cyan3","forestgreen", "goldenrod3", "mediumvioletred"), labels = c("AMS5192","L26","P75063","P75016"))
ggsave("VarAveCompMCF.tiff", plot = last_plot(), path = fig.dir, height = 3, width = 4)

ggplot(data = BFsamp[BFsamp$gRNA == "C305110WA",], aes(x = DIST, y = strain_fitness)) +
  geom_point(aes(color = BFYPADSIG)) +
  theme_bw() +
  xlab("Distance from Start Codon") +
  ylab("Individual Barcode Fitness") +
  scale_color_manual(labels = c("Not \nSignificant","Significant"), values = c("grey", "darkred")) +
  labs(color = "Fitness Effect")
ggsave("C305110WBarComp.tiff", plot = last_plot(), path = fig.dir, height = 3, width = 4)

ggplot(data = BFsamp[BFsamp$gRNA == "C305070WA",], aes(x = DIST, y = strain_fitness)) +
  geom_point(aes(color = BFYPADSIG)) +
  theme_bw() +
  xlab("Distance from Start Codon") +
  ylab("Individual Barcode Fitness") +
  scale_color_manual(labels = c("Not \nSignificant","Significant"), values = c("grey", "darkred")) +
  labs(color = "Fitness Effect")
ggsave("C305070WBarComp.tiff", plot = last_plot(), path = fig.dir, height = 3, width = 4)

ggplot(data = BFsamp[BFsamp$gRNA == "C114070WA",], aes(x = DIST, y = strain_fitness)) +
  geom_point(aes(color = BFYPADSIG)) +
  theme_bw() +
  xlab("Distance from Start Codon") +
  ylab("Individual Barcode Fitness") +
  scale_color_manual(labels = c("Not \nSignificant","Significant"), values = c("grey", "darkred")) +
  labs(color = "Fitness Effect")
ggsave("C114070WBarComp.tiff", plot = last_plot(), path = fig.dir, height = 3, width = 4)

#Grab the weighted average fitness for these 
mean(BFsamp[BFsamp$gRNA == "C114070WA",]$fitness)
mean(BFsamp[BFsamp$gRNA == "C305070WA",]$fitness)
mean(BFsamp[BFsamp$gRNA == "C305110WA",]$fitness)



######################
# Main text and supplementary figure generation

#Plotting distributions of fitness effects in rich media (YPAD)
Facetlab <- c('P75016','P75063',"AMS5192",'L26')
names(Facetlab) <- c('16','63','BF','L26')
fitness_id$Strain <- factor(fitness_id$Strain, levels = c("BF","L26","63","16"))

#Plotting non-targetin fitness distributions separately so that I can inset them in the full distributions
ggplot(fitness_id[fitness_id$Condition == "YPAD" & fitness_id$Timepoint == "T1" & fitness_id$gRNA == "NonTargeting" & fitness_id$Strain == 'BF',], aes(x = AdjFitness)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(-2,1) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("") +
  ylab("")
ggsave("AMS5192NTFitDist.tiff", path = fig.dir, plot = last_plot(), width = 2, height = 2)

ggplot(fitness_id[fitness_id$Condition == "YPAD" & fitness_id$Timepoint == "T1" & fitness_id$gRNA == "NonTargeting" & fitness_id$Strain == '63',], aes(x = AdjFitness)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(-2,1) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("") +
  ylab("")
ggsave("P75063NTFitDist.tiff", path = fig.dir, plot = last_plot(), width = 2, height = 2)

ggplot(fitness_id[fitness_id$Condition == "YPAD" & fitness_id$Timepoint == "T1" & fitness_id$gRNA == "NonTargeting" & fitness_id$Strain == '16',], aes(x = AdjFitness)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(-2,1) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("") +
  ylab("")
ggsave("P75016NTFitDist.tiff", path = fig.dir, plot = last_plot(), width = 2, height = 2)

ggplot(fitness_id[fitness_id$Condition == "YPAD" & fitness_id$Timepoint == "T1" & fitness_id$gRNA == "NonTargeting" & fitness_id$Strain == 'L26',], aes(x = AdjFitness)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  xlim(-2,1) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("") +
  ylab("")
ggsave("L26NTFitDist.tiff", path = fig.dir, plot = last_plot(), width = 2, height = 2)

#Eliminating redundancy for gene-wise fitness estimates by taking one line per gene
gene_only <- fitness_id[,c(c("Samples","gRNA","genestart","geneend","genestrand","fitness","fitness_left","fitness_right","Timepoint","Condition","Hold","Strain","Replicate","AdjFitness"))]
gene_only <- unique(gene_only)

############################
## Calculating means and standard deviations across all three replicates
gene_only_allsamp <- gene_only #Preserving this so I can use it later for calculating residuals for each environment paired to its specific YPAD sample

#A few samples have a 5th replicate that I will remove before taking the mean - they were done because we were worried about quality, but the original replicates look fine.
gene_only <- gene_only[gene_only$Replicate != 5,]

#Also have to remove two low quality samples that do not cluster with the other replicates and had very poor coverage.
gene_only <- gene_only[gene_only$Samples != "T1_FLC_L26.4_S59",]
gene_only <- gene_only[gene_only$Samples != "T1_SDS_16.4_S111",]

gene_only <- gene_only %>%
  group_by(Strain, Condition, Timepoint, gRNA) %>%
  mutate(Mean = mean(AdjFitness), SE = (sd(AdjFitness))/3)

#Removing Non-Targeting controls since fitness has already been adjusted
gene_only_noNT <- gene_only[gene_only$gRNA != "NonTargeting",]

#Removing the now-redundant replicate measures
gene_only_norep <- gene_only_noNT[gene_only_noNT$Replicate == 2,]

# Plotting the overall distributions of fitness effects between backgrounds in rich media
Facetlab <- c('P75016','P75063',"AMS5192",'L26')
names(Facetlab) <- c('16','63','BF','L26')
gene_only_norep$Strain <- factor(gene_only_norep$Strain, levels = c("BF","L26","63","16"))
ggplot(gene_only_norep[gene_only_norep$Condition == "YPAD" & gene_only_norep$Timepoint == "T1" & gene_only_norep$gRNA != "NonTargeting",], aes(x = Mean)) +
  geom_density(fill = "grey") +
  geom_vline(xintercept = 0) +
  facet_wrap(~ Strain, nrow = 2, labeller = labeller(Strain = Facetlab)) +
  xlim(-2,1) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.title = element_text(size = 10, color = "black"), axis.text = element_text(size = 10, color = "black")) +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Fitness \n (Weighted Mean Gene Fitness)") +
  ylab("Density")
ggsave("FitDist.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 4)

# Pivoting the dataframe to better compare samples to each other.
CompmultirepsMean <- pivot_wider(gene_only_norep, id_cols = c('gRNA'), names_from = c('Strain','Condition','Timepoint'), values_from = 'Mean')
CompmultirepsSD <- pivot_wider(gene_only_norep, id_cols = c('gRNA'), names_from = c('Strain','Condition','Timepoint'), values_from = 'SE')
colnames(CompmultirepsMean)[2:length(colnames(CompmultirepsMean))] <- paste0(colnames(CompmultirepsMean)[2:length(colnames(CompmultirepsMean))],"_Mean")
colnames(CompmultirepsSD)[2:length(colnames(CompmultirepsSD))] <- paste0(colnames(CompmultirepsSD)[2:length(colnames(CompmultirepsSD))],"_SD")

CompmultirepsMean <- left_join(CompmultirepsMean, CompmultirepsSD, by = "gRNA")

######################
#Calling significantly different fitness in YPAD from the inoculum for each background
Strainvec <- c("BF","63","16","L26")
Conditionvec <- c("37C", "4NQO", "FLC", "MCF", "SDS", "pH4", "pH8")
Testpiv <- pivot_wider(gene_only_noNT, id_cols = c("gRNA"), names_from = Samples, values_from = AdjFitness)

YPADGenes <- list()
for (j in 1:length(Strainvec)) {
  Testpivselect <- Testpiv[,c(grep(paste0("T1_YPAD_",Strainvec[j]), colnames(Testpiv), value = TRUE))]
  Testpivselect$gRNA <- Testpiv$gRNA
  Testmat <- as.matrix(Testpivselect[,!colnames(Testpivselect) %in% c("barcode","gRNA")])
  rownames(Testmat) <- Testpivselect$gRNA
  NTVec <- gene_only_allsamp[gene_only_allsamp$Samples %in% colnames(Testpivselect) & gene_only_allsamp$gRNA == "NonTargeting" & !is.na(gene_only_allsamp$AdjFitness),"AdjFitness"]
  groups <- factor(c(rep("group 1", ncol(Testmat)), rep("group 2", nrow(NTVec))))
  n <- 1
  NTVecMat <- t(NTVec)
  while(n < nrow(Testmat)) {
    NTVecMat <- rbind(NTVecMat, t(NTVec))
    n <- n + 1
  }
  Testmat <- cbind(Testmat, NTVecMat)
  Testtest <- mod.t.test(Testmat, group = groups, adjust.method = "BH",
                         sort.by = "p")
  Testtest$gRNA <- rownames(Testtest)
  HITS <- unique(Testtest[Testtest$adj.p.value <= 0.05,"gRNA"])
  YPADGenes[[j]] <- HITS
}
names(YPADGenes) <- Strainvec

# Venn diagram
tiff(paste0(fig.dir,"VennYPAD.tiff"), width = 2000, height = 1500, units = "px", res = 300)
venn.diagram(x = list(YPADGenes[["BF"]], YPADGenes[["L26"]],YPADGenes[["63"]],YPADGenes[["16"]]), filename = NULL, 
             category.names = c("AMS5192","L26","P75063","P75016"), 
             fill = c("cyan3","forestgreen", "goldenrod3", "mediumvioletred"), # Colors for the circles
             alpha = 0.5, # Opacity
             cex = 1.5, # Font size of the numbers
             cat.cex = 1.5) # Font size of the category labels
dev.off()

#Comparing fitness values between backgrounds in rich media
CompmultirepsMean %>%
  select(gRNA, BF_YPAD_T1_Mean, BF_YPAD_T1_SD, `16_YPAD_T1_Mean`, `16_YPAD_T1_SD`) %>%
  mutate(Categ = if_else(gRNA %in% YPADGenes$`16` & !gRNA %in% YPADGenes$BF, "In16", "NDE"),
         Categ = if_else(!gRNA %in% YPADGenes$`16` & gRNA %in% YPADGenes$BF, "InBF", Categ),
         Categ = if_else(gRNA %in% YPADGenes$`16` & gRNA %in% YPADGenes$BF, "Both", Categ)) %>%
  ggplot(aes(x = BF_YPAD_T1_Mean, y = `16_YPAD_T1_Mean`)) +
  geom_point(aes(color = Categ)) +
  geom_errorbarh(aes(xmin = BF_YPAD_T1_Mean - BF_YPAD_T1_SD, xmax = BF_YPAD_T1_Mean + BF_YPAD_T1_SD, color = Categ)) +
  geom_errorbar(aes(ymin = `16_YPAD_T1_Mean` - `16_YPAD_T1_SD`, ymax = `16_YPAD_T1_Mean` + `16_YPAD_T1_SD`, color = Categ)) +
  scale_color_manual(values = c("firebrick","cyan3","grey")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  #theme(legend.position = "none") +
  ylab("P75016 Mean Fitness") +
  xlab("AMS5192 Mean Fitness") +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("BF16CompYPAD.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3) 

#What does a rank order correlation look like for these strains?
cor.test(CompmultirepsMean$BF_YPAD_T1_Mean, CompmultirepsMean$`16_YPAD_T1_Mean`, method = 'spearman')
cor.test(CompmultirepsMean$`63_YPAD_T1_Mean`, CompmultirepsMean$`16_YPAD_T1_Mean`, method = 'spearman') #These two are significant: p-value = 2.846e-15, rho = 0.2773747 
cor.test(CompmultirepsMean$BF_YPAD_T1_Mean, CompmultirepsMean$L26_YPAD_T1_Mean, method = 'spearman')
cor.test(CompmultirepsMean$BF_YPAD_T1_Mean, CompmultirepsMean$`63_YPAD_T1_Mean`, method = 'spearman')

CompmultirepsMean %>%
  select(gRNA, BF_YPAD_T1_Mean, BF_YPAD_T1_SD, `63_YPAD_T1_Mean`, `63_YPAD_T1_SD`) %>%
  mutate(Categ = if_else(gRNA %in% YPADGenes$`63` & !gRNA %in% YPADGenes$BF, "In63", "NDE"),
         Categ = if_else(!gRNA %in% YPADGenes$`63` & gRNA %in% YPADGenes$BF, "InBF", Categ),
         Categ = if_else(gRNA %in% YPADGenes$`63` & gRNA %in% YPADGenes$BF, "Both", Categ)) %>%
  ggplot(aes(x = BF_YPAD_T1_Mean, y = `63_YPAD_T1_Mean`, label = gRNA)) +
  geom_point(aes(color = Categ)) +
  #geom_text() +
  geom_errorbarh(aes(xmin = BF_YPAD_T1_Mean - BF_YPAD_T1_SD, xmax = BF_YPAD_T1_Mean + BF_YPAD_T1_SD, color = Categ)) +
  geom_errorbar(aes(ymin = `63_YPAD_T1_Mean` - `63_YPAD_T1_SD`, ymax = `63_YPAD_T1_Mean` + `63_YPAD_T1_SD`, color = Categ)) +
  scale_color_manual(values = c("firebrick","goldenrod3","cyan3","grey")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  #theme(legend.position = "none") +
  ylab("P75063 Mean Fitness") +
  xlab("AMS5192 Mean Fitness") +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("BF63CompYPAD.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3) 

CompmultirepsMean %>%
  select(gRNA, BF_YPAD_T1_Mean, BF_YPAD_T1_SD, `L26_YPAD_T1_Mean`, `L26_YPAD_T1_SD`) %>%
  mutate(Categ = if_else(gRNA %in% YPADGenes$`L26` & !gRNA %in% YPADGenes$BF, "InL26", "NDE"),
         Categ = if_else(!gRNA %in% YPADGenes$`L26` & gRNA %in% YPADGenes$BF, "InBF", Categ),
         Categ = if_else(gRNA %in% YPADGenes$`L26` & gRNA %in% YPADGenes$BF, "Both", Categ)) %>%
  ggplot(aes(x = BF_YPAD_T1_Mean, y = `L26_YPAD_T1_Mean`)) +
  geom_point(aes(color = Categ)) +
  geom_errorbarh(aes(xmin = BF_YPAD_T1_Mean - BF_YPAD_T1_SD, xmax = BF_YPAD_T1_Mean + BF_YPAD_T1_SD, color = Categ)) +
  geom_errorbar(aes(ymin = `L26_YPAD_T1_Mean` - `L26_YPAD_T1_SD`, ymax = `L26_YPAD_T1_Mean` + `L26_YPAD_T1_SD`, color = Categ)) +
  scale_color_manual(values = c("firebrick","cyan3","forestgreen","grey")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  #theme(legend.position = "none") +
  ylab("L26 Mean Fitness") +
  xlab("AMS5192 Mean Fitness") +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("BFL26CompYPAD.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3) 

CompmultirepsMean %>%
  select(gRNA, `16_YPAD_T1_Mean`, `16_YPAD_T1_SD`, `63_YPAD_T1_Mean`, `63_YPAD_T1_SD`) %>%
  mutate(Categ = if_else(gRNA %in% YPADGenes$`63` & !gRNA %in% YPADGenes$`16`, "In63", "NDE"),
         Categ = if_else(!gRNA %in% YPADGenes$`63` & gRNA %in% YPADGenes$`16`, "In16", Categ),
         Categ = if_else(gRNA %in% YPADGenes$`63` & gRNA %in% YPADGenes$`16`, "Both", Categ)) %>%
  ggplot(aes(x = `16_YPAD_T1_Mean`, y = `63_YPAD_T1_Mean`)) +
  geom_point(aes(color = Categ)) +
  geom_errorbarh(aes(xmin = `16_YPAD_T1_Mean` - `16_YPAD_T1_SD`, xmax = `16_YPAD_T1_Mean` + `16_YPAD_T1_SD`, color = Categ)) +
  geom_errorbar(aes(ymin = `63_YPAD_T1_Mean` - `63_YPAD_T1_SD`, ymax = `63_YPAD_T1_Mean` + `63_YPAD_T1_SD`, color = Categ)) +
  scale_color_manual(values = c("firebrick","goldenrod3","grey")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  #theme(legend.position = "none") +
  ylab("P75063 Mean Fitness") +
  xlab("P75016 Mean Fitness") +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave("1663CompYPAD.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3) 

###########################
# Comparing all environments to YPAD using residuals
# Grabbing the residuals for each linear model - each replicate is done separately because replicates are paired environments

gene_only_allsamp <- gene_only_allsamp[gene_only_allsamp$gRNA != "NonTargeting",]

CompReps <- pivot_wider(gene_only_allsamp, id_cols = c('gRNA'), names_from = c('Strain','Condition','Timepoint','Replicate'), values_from = 'AdjFitness')

for (i in 1:length(Strainvec)) {
  x <- CompReps[,c(grep(Strainvec[i],colnames(CompReps)))]
  for (j in 1:length(Conditionvec)) {
    for (h in 2:4) { 
      if(!all(is.na(CompReps[,c(grep(paste0(Conditionvec[j],"_T1_",h), colnames(x), value = TRUE))]))) {
        xsub <- x[,c(grep(paste0(Conditionvec[j],"_T1_",h), colnames(x), value = TRUE), grep(paste0("YPAD_T1_",h), colnames(x), value = TRUE))]
        colnames(xsub) <- c("V1","V2")
        lmobj <- lm(xsub$V1 ~ xsub$V2)
        CompReps[,paste0(Strainvec[i],"_",Conditionvec[j],"_Resid", h)] <- lmobj$residuals
      }
    }
  }
}

TopResids <- list()
for (i in 1:length(Conditionvec)) {
  TopResids[[i]] <- list()
  for (j in 1:length(Strainvec)) {
    CompRepsselect <- CompReps[,c(grep(paste0(Strainvec[j], "_", Conditionvec[i], "_Resid"), colnames(CompReps), value = TRUE))]
    CompRepsselect$gRNA <- CompReps$gRNA
    Testmat <- as.matrix(CompRepsselect[,!colnames(CompRepsselect) %in% c("barcode","gRNA")])
    rownames(Testmat) <- CompRepsselect$gRNA
    Testtest <- mod.t.test(Testmat, adjust.method = "BH",
                           sort.by = "p")
    Testtest$gRNA <- rownames(Testtest)
    HITS <- unique(Testtest[Testtest$adj.p.value <= 0.05,"gRNA"])
    TopResids[[i]][[j]] <- HITS
  }
}
names(TopResids) <- Conditionvec
for (i in 1:length(TopResids)) {
  names(TopResids[[i]]) <- Strainvec
}
TopResidVec <- unlist(TopResids)
TopResidVec <- unique(TopResidVec)

# For plotting, getting the residuals between the mean values
for (i in 1:length(Strainvec)) {
  x <- CompmultirepsMean[,c(grep(Strainvec[i],colnames(CompmultirepsMean)))]
  for (j in 1:length(Conditionvec)) {
    xsub <- x[,c(grep(paste0(Conditionvec[j],"_T1_Mean"), colnames(x)), grep("YPAD_T1_Mean", colnames(x)))]
    colnames(xsub) <- c("V1","V2")
    lmobj <- lm(xsub$V1 ~ xsub$V2)
    CompmultirepsMean[,paste0(Strainvec[i],"_",Conditionvec[j],"_Resid")] <- lmobj$residuals
  }
}

Resids <- CompmultirepsMean[,grep("Resid", colnames(CompmultirepsMean))]
# Getting the columns into a logical order
Resids <- Resids[,c(grep("FLC", colnames(Resids)),grep("MCF", colnames(Resids)),grep("SDS", colnames(Resids)), grep("37C", colnames(Resids)),grep("4NQO", colnames(Resids)),grep("pH4", colnames(Resids)),grep("pH8", colnames(Resids)))]
Resids <- as.data.frame(Resids)
rownames(Resids) <- CompmultirepsMean$gRNA

#######################
## Visualizing residual fitness values

TopResidsdf <- Resids[rownames(Resids) %in% TopResidVec,]

# I want to also impose a fold change cut-off for the heatmap to visualize it better. 
FCPassGenes <- apply(TopResidsdf, 1, function(x) { any(abs(x) >= 0.3) }) 
TopResidsdflim <- TopResidsdf[which(FCPassGenes),]
TopResidslimtoclust <- TopResidsdflim[!rownames(TopResidsdflim) %in% c("C305220WA","C403180WA"),] #Excluding these because they dominate the clustering too much
row_distance <- dist(TopResidslimtoclust, method = "euclidean") #Calculating clustering before inserting NAs for nonsignificant genes
row_cluster <- hclust(row_distance)
TopResidslimNA <- TopResidslimtoclust
for (i in 1:ncol(TopResidslimtoclust)) {
  TopResidslimNA[,i] <- ifelse(rownames(TopResidslimtoclust) %in% TopResids[[paste0(strsplit(colnames(TopResidslimtoclust[i]), "_")[[1]][2])]][[paste0(strsplit(colnames(TopResidslimtoclust[i]), "_")[[1]][1])]], TopResidslimtoclust[,i], NA)
}
# Defining the heatmap color pallette
my_palette_up <- (colorRampPalette(brewer.pal(9, "Reds"))(25))
my_palette_down <- rev(colorRampPalette(brewer.pal(9, "Blues"))(25))
my_palette <- c(my_palette_down[3:25], my_palette_up[1:25])
a <- pheatmap(TopResidslimNA, scale = "none", cluster_cols = FALSE, cluster_rows = row_cluster, show_rownames = FALSE, show_colnames = TRUE, color = my_palette, na_col = "grey")
ggsave("Residmappoint3NA.tiff", path = fig.dir, plot = a, width = 5, height = 7)
ggsave("Residmappoint3NA.pdf", path = fig.dir, plot = a, width = 5, height = 7)

b <- pheatmap(TopResidslimtoclust, scale = "none", cluster_cols = FALSE, cluster_rows = row_cluster, show_rownames = FALSE, show_colnames = TRUE, color = my_palette, na_col = "grey")
ggsave("Residmappoint3all.tiff", path = fig.dir, plot = b, width = 5, height = 7)
ggsave("Residmappoint3all.pdf", path = fig.dir, plot = b, width = 5, height = 7)

CDRNCPNA <- TopResidsdf[rownames(TopResidsdf) %in% c("C305220WA","C403180WA"),]
for (i in 1:ncol(CDRNCPNA)) {
    CDRNCPNA[,i] <- ifelse(rownames(CDRNCPNA) %in% TopResids[[paste0(strsplit(colnames(CDRNCPNA[i]), "_")[[1]][2])]][[paste0(strsplit(colnames(CDRNCPNA[i]), "_")[[1]][1])]], CDRNCPNA[,i], NA)
}

c <- pheatmap(CDRNCPNA, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = my_palette_up, na_col = "grey")
ggsave("ResidmapNCPCDR.tiff", path = fig.dir, plot = c, width = 6, height = 0.5)
ggsave("ResidmapNCPCDR.pdf", path = fig.dir, plot = c, width = 6, height = 0.5)

# Heatmap of numbers of significant genes in each background by each environment
NumberComp <- matrix(nrow = length(Strainvec), ncol = length(Conditionvec))
rownames(NumberComp) <- Strainvec
colnames(NumberComp) <- Conditionvec
for (i in 1:length(TopResids)) {
  for (j in 1:length(TopResids[[i]])) {
    NumberComp[j,i] <- length(TopResids[[i]][[j]])
  }
}
NumberComp <- NumberComp[,c(grep("FLC", colnames(NumberComp)),grep("MCF", colnames(NumberComp)),grep("SDS", colnames(NumberComp)),grep("37C", colnames(NumberComp)),grep("4NQO", colnames(NumberComp)),grep("pH4", colnames(NumberComp)),grep("pH8", colnames(NumberComp)))]

d <- pheatmap(NumberComp, scale = "none", color = brewer.pal(name = "Purples", n = 9), cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_color = "black", fontsize_number = 12, number_format = "%.0f")
ggsave("NumberCompHeatmap.tiff", path = fig.dir, width = 4, height = 4, plot = d)

# How does this compare to each background's sensitivity to the environment it is being grown in?
Doublingsdat <- read.table(paste0(data.dir,"bigbatchT1Doublings.txt"), sep = "\t", header = 1)
#Oops, have to get rid of Xs before the numbers
colnames(Doublingsdat) <- gsub("X","", colnames(Doublingsdat))
#Calculating the difference in doublings in each environment and in YPAD
for (i in 4:ncol(Doublingsdat)) {
  Doublingsdat[,i] <- Doublingsdat[,"YPAD"] - Doublingsdat[,i]
}
Doublingsdatmelt <- pivot_longer(Doublingsdat, cols = c("YPAD","37C","pH4","pH8","SDS","FLC","4NQO","MCF"))

#Calculating summary stats:
Doublingsdatsum <- Doublingsdatmelt %>%
  group_by(Strain, name) %>%
  summarise(Mean = mean(value), SD = sd(value))
colnames(Doublingsdatsum)[2] <- "Condition"

Doublingscomp <- pivot_longer(as.data.frame(NumberComp), cols = colnames(NumberComp), values_to = "NumSig", names_to = "Condition")
Doublingscomp$Strain <- rep(rownames(NumberComp), each = ncol(NumberComp))
Doublingscomp <- left_join(Doublingscomp, Doublingsdatsum, by = c("Strain","Condition"))

ggplot(Doublingscomp, aes(x = Mean, y = NumSig)) +
  geom_point(aes(color = Condition, shape = Strain), size = 4) +
  geom_errorbarh(aes(xmin = Mean + SD, xmax = Mean - SD, color = Condition)) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  xlab("Sensitivity\n Doublings in stress - doublings in rich medium") +
  ylab("Number of genes\n with significant fitness effects") +
  scale_color_manual(values = c("firebrick","darkorange","goldenrod3","forestgreen","cyan3",'violetred',"purple4"))

# Interesting, how about the number that have positive effects separate from those with negative effects?
PosResids <- list()
for (i in 1:length(Conditionvec)) {
  PosResids[[i]] <- list()
  for (j in 1:length(Strainvec)) {
    CompRepsselect <- CompReps[,c(grep(paste0(Strainvec[j], "_", Conditionvec[i], "_Resid"), colnames(CompReps), value = TRUE))]
    CompRepsselect$gRNA <- CompReps$gRNA
    Testmat <- as.matrix(CompRepsselect[,!colnames(CompRepsselect) %in% c("barcode","gRNA")])
    rownames(Testmat) <- CompRepsselect$gRNA
    Testtest <- mod.t.test(Testmat, adjust.method = "BH",
                           sort.by = "p")
    Testtest$gRNA <- rownames(Testtest)
    HITS <- unique(Testtest[Testtest$adj.p.value <= 0.05 & Testtest$mean > 0,"gRNA"])
    PosResids[[i]][[j]] <- HITS
  }
}
names(PosResids) <- Conditionvec
for (i in 1:length(PosResids)) {
  names(PosResids[[i]]) <- Strainvec
}
PosResidVec <- unlist(PosResids)
PosResidVec <- unique(PosResidVec)

NumUpComp <- matrix(nrow = length(Strainvec), ncol = length(Conditionvec))
rownames(NumUpComp) <- Strainvec
colnames(NumUpComp) <- Conditionvec
for (i in 1:length(PosResids)) {
  for (j in 1:length(PosResids[[i]])) {
    NumUpComp[j,i] <- length(PosResids[[i]][[j]])
  }
}
NumUpComp <- NumUpComp[,c(grep("FLC", colnames(NumUpComp)),grep("MCF", colnames(NumUpComp)),grep("SDS", colnames(NumUpComp)),grep("37C", colnames(NumUpComp)),grep("4NQO", colnames(NumUpComp)),grep("pH4", colnames(NumUpComp)),grep("pH8", colnames(NumUpComp)))]

Doublingscomp2 <- pivot_longer(as.data.frame(NumUpComp), cols = colnames(NumUpComp), values_to = "NumSig", names_to = "Condition")
Doublingscomp2$Strain <- rep(rownames(NumUpComp), each = ncol(NumUpComp))
Doublingscomp <- left_join(Doublingscomp, Doublingscomp2, by = c("Strain","Condition"))

ggplot(Doublingscomp, aes(x = Mean, y = NumSig.y)) +
  geom_point(aes(color = Condition, shape = Strain), size = 4) +
  geom_errorbarh(aes(xmin = Mean + SD, xmax = Mean - SD, color = Condition)) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  xlab("Difference in cell doublings") +
  ylab("Number of genes\n with positive fitness effects") +
  scale_color_manual(values = c("firebrick","darkorange","goldenrod3","forestgreen","cyan3",'violetred',"purple4"))
ggsave("DoubvsNum.tiff", path = fig.dir, plot = last_plot(), width = 4.5, height = 4)

#Getting the linear model stats
doublm <- lm(Doublingscomp$NumSig.y ~ Doublingscomp$Mean)
summary(doublm)

# Looking at significant fitness genes shared across environments and strains
SHAREDGENES <- function(x) {
  EnvComp <- fromList(TopResids[[x]])
  EnvComp$Gene <- rownames(EnvComp)
  CompDF <- pivot_longer(EnvComp, cols = colnames(EnvComp)[1:4])
  for (i in 1:nrow(CompDF)) {
      CompDF[i, "Resid"] <- TopResidsdf[rownames(TopResidsdf) == CompDF$Gene[i],colnames(TopResidsdf) == paste0(CompDF$name[i],"_",x,"_Resid")]
  }
  CompDFdrop <- CompDF[CompDF$value != 0,]
  CompDFdrop <- CompDFdrop %>%
    group_by(Gene) %>%
    mutate("Sum" = sum(value))
  CompDFdrop$Environment <- rep(x, nrow(CompDFdrop))
  return(CompDFdrop)
}

StrainComplist <- lapply(Conditionvec, SHAREDGENES)
StrainComp <- do.call(rbind, StrainComplist)

ggplot(StrainComp, aes(x = Sum, y = abs(Resid))) +
  geom_jitter(height = 0, width = 0.2) +
  #geom_smooth(method = "glm") +
  stat_summary(fun = median, color = "deeppink1") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  #scale_color_manual(values = c("firebrick","darkorange","goldenrod3","forestgreen","cyan3",'violetred',"purple4")) +
  ylab("Absolute Fitness Effect") +
  xlab("Number of Backgrounds") +
  ylim(0, 2)
ggsave("NoStrainComp.tiff", plot = last_plot(), width = 3, height = 3, path = fig.dir)

#Testing if these bins are significantly different
aov <- anova(lm(abs(Resid) ~ Sum, data = StrainComp))

#Comparing all residuals between MCF and SDS
Resids$Gene <- rownames(Resids)
Residsmelt <- pivot_longer(Resids, cols = colnames(Resids)[1:ncol(Resids)-1], names_to = "Sample", values_to = "Fitness")
Residsmelt$Environment <- sapply(strsplit(Residsmelt$Sample, "_"),`[`,2)
Residsmelt$Strain <- sapply(strsplit(Residsmelt$Sample, "_"),`[`,1)

CompmultirepsSD$Gene <- CompmultirepsSD$gRNA
Resids2melt <- pivot_longer(CompmultirepsSD[,c(grep("T1_SD", colnames(CompmultirepsSD), value = TRUE),"Gene")], cols = colnames(CompmultirepsSD)[grep("T1_SD", colnames(CompmultirepsSD))], names_to = "Sample", values_to = "SD")
Resids2melt$Environment <- sapply(strsplit(Resids2melt$Sample, "_"),`[`,2)
Resids2melt$Strain <- sapply(strsplit(Resids2melt$Sample, "_"),`[`,1)

ResidsSDs <- left_join(Residsmelt, Resids2melt, by = c("Gene","Environment","Strain"))
Residspivot <- pivot_wider(ResidsSDs, names_from = Environment, values_from = c(Fitness, SD), id_cols = c(Gene, Strain))

ggplot(data = Residspivot, aes(x = Fitness_MCF, y = Fitness_SDS)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = Strain)) +
  geom_errorbar(aes(y = Fitness_SDS, x = Fitness_MCF, ymax = Fitness_SDS + SD_SDS, ymin = Fitness_SDS - SD_SDS, color = Strain)) +
  geom_errorbarh(aes(y = Fitness_SDS, x = Fitness_MCF, xmax = Fitness_MCF + SD_MCF, xmin = Fitness_MCF - SD_MCF, color = Strain)) +
  scale_color_manual(values = c("mediumvioletred","goldenrod3","cyan3","forestgreen"), labels = c("P75016","P75063",'AMS5192',"L26")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Fitness in MCF") +
  ylab("Fitness in SDS")
ggsave("AllMCFSDSComp.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)
AllMCFSDSMod <- lm(Residspivot$Fitness_SDS ~ Residspivot$Fitness_MCF)
summary(AllMCFSDSMod)

#Comparing residuals between MCF and SDS for just genes that are significant in both.
for (i in 1:nrow(Residsmelt)) {
  Residsmelt[i,"SIG"] <- ifelse(Residsmelt$Gene[i] %in% TopResids[[Residsmelt$Environment[i]]][[Residsmelt$Strain[i]]], "SIG", "NS")
}
SIGpivot <- pivot_wider(Residsmelt, names_from = Environment, values_from = c(SIG), id_cols = c(Gene, Strain))
colnames(SIGpivot)[3:9] <- paste0("SIG_",colnames(SIGpivot)[3:9])
Residspivot <- left_join(Residspivot, SIGpivot, by = c("Gene","Strain"))

ggplot(data = Residspivot[Residspivot$SIG_SDS == "SIG" & Residspivot$SIG_MCF == "SIG",], aes(x = Fitness_MCF, y = Fitness_SDS)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = Strain)) +
  geom_errorbar(aes(y = Fitness_SDS, x = Fitness_MCF, ymax = Fitness_SDS + SD_SDS, ymin = Fitness_SDS - SD_SDS, color = Strain)) +
  geom_errorbarh(aes(y = Fitness_SDS, x = Fitness_MCF, xmax = Fitness_MCF + SD_MCF, xmin = Fitness_MCF - SD_MCF, color = Strain)) +
  scale_color_manual(values = c("mediumvioletred","goldenrod3","cyan3","forestgreen"), labels = c("P75016","P75063",'AMS5192',"L26")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Fitness in MCF") +
  ylab("Fitness in SDS") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')
ggsave("SigMCFSDSComp.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)
SDSMCFSigMod <- lm(Residspivot[Residspivot$SIG_SDS == "SIG" & Residspivot$SIG_MCF == "SIG",]$Fitness_SDS ~ Residspivot[Residspivot$SIG_SDS == "SIG" & Residspivot$SIG_MCF == "SIG",]$Fitness_MCF)
summary(SDSMCFSigMod)

#Now doing the same for FLC and MCF
ggplot(data = Residspivot[Residspivot$SIG_FLC == "SIG" & Residspivot$SIG_MCF == "SIG",], aes(x = Fitness_MCF, y = Fitness_FLC)) +
  geom_smooth(method = "lm") +
  geom_point(aes(color = Strain)) +
  geom_errorbar(aes(y = Fitness_FLC, x = Fitness_MCF, ymax = Fitness_FLC + SD_FLC, ymin = Fitness_FLC - SD_FLC, color = Strain)) +
  geom_errorbarh(aes(y = Fitness_FLC, x = Fitness_MCF, xmax = Fitness_MCF + SD_MCF, xmin = Fitness_MCF - SD_MCF, color = Strain)) +
  scale_color_manual(values = c("mediumvioletred","goldenrod3","cyan3","forestgreen"), labels = c("P75016","P75063",'AMS5192',"L26")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Fitness in MCF") +
  ylab("Fitness in FLC") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ylim(-1,1) #This excludes the CDR1 that is only significant in micafungin in P75016
ggsave("SigMCFFLCComp.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)
FLCMCFSigMod <- lm(Residspivot[Residspivot$SIG_FLC == "SIG" & Residspivot$SIG_MCF == "SIG",]$Fitness_FLC ~ Residspivot[Residspivot$SIG_FLC == "SIG" & Residspivot$SIG_MCF == "SIG",]$Fitness_MCF)
summary(FLCMCFSigMod)

# Scatterplots comparing fitness between YPAD and other environments (visualizing the basis of the residuals)
ggplot(CompmultirepsMean, aes(x = BF_YPAD_T1_Mean, y = BF_FLC_T1_Mean)) +
  geom_point(aes(color = gRNA %in% TopResids$FLC$BF)) +
  geom_errorbarh(aes(xmin = BF_YPAD_T1_Mean - BF_YPAD_T1_SD, xmax = BF_YPAD_T1_Mean + BF_YPAD_T1_SD, color = gRNA %in% TopResids$FLC$BF)) +
  geom_errorbar(aes(ymin = BF_FLC_T1_Mean - BF_FLC_T1_SD, ymax = BF_FLC_T1_Mean + BF_FLC_T1_SD, color = gRNA %in% TopResids$FLC$BF)) +
  scale_color_manual(values = c("#839faaff","#ca1d58ff")) +
  #xlim(-2.25,1) +
  #ylim(-2.25,3) +
  geom_smooth(method = "lm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(legend.position = "none") +
  ylab("OE Fitness in FLC") +
  xlab("OE Fitness in Rich Media") +
  ggtitle("AMS5192") +
  geom_abline(slope = 1, intercept = 0)
ggsave("BFYPADFLC.tiff", path = fig.dir, plot = last_plot(), width = 3, height = 3)

##########################
# Upset plots for differences between strains/environments

TopResidsT <- transpose(TopResids)
elements <- lapply(TopResidsT$BF, function(x) {
  is_empty(x)
})
upset(fromList(TopResidsT$BF[which(elements == FALSE)]), sets = c(names(TopResidsT$BF[which(elements == FALSE)])), text.scale = 1.5)

tiff(filename = paste0(fig.dir,"BFENVUPSET.tiff"), width = 1400, height = 1000, units = "px", res = 300)
upset(fromList(TopResidsT$BF[which(elements == FALSE)]), sets = c(names(TopResidsT$BF[which(elements == FALSE)])), text.scale = 1.5)
dev.off()

elements <- lapply(TopResidsT$L26, function(x) {
  is_empty(x)
})
upset(fromList(TopResidsT$L26[which(elements == FALSE)]), sets = c(names(TopResidsT$L26[which(elements == FALSE)])))

tiff(filename = paste0(fig.dir,"L26ENVUPSET.tiff"), width = 1200, height = 1000, units = "px", res = 300)
upset(fromList(TopResidsT$L26[which(elements == FALSE)]), sets = c(names(TopResidsT$L26[which(elements == FALSE)])), text.scale = 1.5)
dev.off()

elements <- lapply(TopResidsT$`63`, function(x) {
  is_empty(x)
})
upset(fromList(TopResidsT$`63`[which(elements == FALSE)]), sets = c(names(TopResidsT$`63`[which(elements == FALSE)])), text.scale = 1.5)
P63Upset <- fromList(TopResidsT$`63`[which(elements == FALSE)])

tiff(filename = paste0(fig.dir,"P63ENVUPSET.tiff"), width = 1200, height = 1000, units = "px", res = 300)
upset(fromList(TopResidsT$`63`[which(elements == FALSE)]), sets = c(names(TopResidsT$`63`[which(elements == FALSE)])), text.scale = 1.5)
dev.off()

elements <- lapply(TopResidsT$`16`, function(x) {
  is_empty(x)
})
upset(fromList(TopResidsT$`16`[which(elements == FALSE)]), sets = c(names(TopResidsT$`16`[which(elements == FALSE)])), text.scale = 1.5)
P16Upset <- fromList(TopResidsT$`16`[which(elements == FALSE)])

tiff(filename = paste0(fig.dir,"P16ENVUPSET.tiff"), width = 1200, height = 1000, units = "px", res = 300)
upset(fromList(TopResidsT$`16`[which(elements == FALSE)]), sets = c(names(TopResidsT$`16`[which(elements == FALSE)])), text.scale = 1.5)
dev.off()

# Looking at upset-plots that are done environment-wise - making this a loop - 
# Unfortunately, the loop isn't working. Have to reset i manually and tab through...
for (i in 1:length(Conditionvec)) {
  set <- TopResids[[i]]
  elements <- lapply(set, function(x) {
    is_empty(x)
  })
  tiff(filename = paste0(fig.dir,names(TopResids)[i],"STRUPSET.tiff"), width = 1200, height = 1000, units = "px", res = 300)
  upset(fromList(TopResids[[i]][which(elements == FALSE)]), sets = c(names(TopResids[[i]][which(elements == FALSE)])), text.scale = 1.5)
  dev.off()
}

###########################
# Connecting back to CNV positions

# Reading in gene positional information
library(ape) #If this library is loaded with those above it will produce a conflict with the source code
A22features <- "/C_albicans_SC5314_version_A22-s07-m01-r213_features.gff"
A22features <- read.gff(paste0(data.dir,A22features))
getAttributeField <- function (x, field, attrsep = ";") { #From the 'ballgown' package, which appears to be deprecated
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
A22features$A22_name <- getAttributeField(A22features$attributes, field = "ID") #This is in a different format than the names in the datafile
A22features$gRNA <- gsub("_","", A22features$A22_name)
Positionkey <- A22features[,c("gRNA", "seqid", "start")]
Positionkey$seqid <- as.character(Positionkey$seqid)
Positionkey$Chrom <- ifelse(Positionkey$seqid == "Ca22chr1A_C_albicans_SC5314", "Chr1", Positionkey$seqid)
Positionkey$Chrom <- ifelse(Positionkey$seqid == "Ca22chr3A_C_albicans_SC5314" & Positionkey$start < 1000000, "Chr3L", Positionkey$Chrom)
Positionkey$Chrom <- ifelse(Positionkey$seqid == "Ca22chr3A_C_albicans_SC5314" & Positionkey$start > 1000000, "Chr3R", Positionkey$Chrom)
Positionkey$Chrom <- ifelse(Positionkey$seqid == "Ca22chr4A_C_albicans_SC5314", "Chr4", Positionkey$Chrom)

colnames(Positionkey)[1] <- "Gene"
ResidsInfo <- left_join(Residsmelt, Positionkey, by = "Gene")
ResidsInfo$Environment <- as.character(ResidsInfo$Environment)
colnames(ResidsSDs)[which(colnames(ResidsSDs) == "Sample.x")] <- "Sample"
ResidsInfo <- left_join(ResidsInfo, ResidsSDs[,c("Gene","Sample","SD")], by = c("Gene","Sample"))

# Creating a variable to be able to plot all and color by strain or environment where it is significant
for (i in 1:nrow(ResidsInfo)) {
  ResidsInfo[i,"SIGStrain"] <- ifelse(ResidsInfo[i,"Gene"] %in% TopResids[[paste(ResidsInfo[i,"Environment"])]][[paste(ResidsInfo[i,"Strain"])]], ResidsInfo[i,"Strain"], "NS")
  ResidsInfo[i,"SIGEnv"] <- ifelse(ResidsInfo[i,"Gene"] %in% TopResids[[paste(ResidsInfo[i,"Environment"])]][[paste(ResidsInfo[i,"Strain"])]], ResidsInfo[i,"Environment"], "NS")
}

ggplot(data = ResidsInfo[ResidsInfo$Environment == "FLC",], aes(x = start/1000000, y = Fitness)) +
  geom_point(aes(color = SIGStrain), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness-SD, ymax = Fitness+SD, color = SIGStrain), alpha = 0.7) +
  geom_point(data = ResidsInfo[ResidsInfo$SIGEnv != "NS" & ResidsInfo$Environment == "FLC",], aes(x = start/1000000, y = Fitness, color = SIGStrain), alpha = 0.7, size = 2) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "white")) +
  xlab("Chromosomal Position (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("mediumorchid","goldenrod3","cyan3","forestgreen","grey"), labels = c("P75016","P75063","AMS5192","L26","NS")) +
  labs(color = "Background")
ggsave("FLCALLPOS.tiff", path = fig.dir, plot = last_plot(), width = 7.5, height = 3)

#Zoomed in version
ggplot(data = ResidsInfo[ResidsInfo$Environment == "FLC",], aes(x = start/1000000, y = Fitness, label = Gene)) +
  geom_point(aes(color = SIGStrain), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness-SD, ymax = Fitness+SD, color = SIGStrain), alpha = 0.7) +
  #geom_text() +
  geom_point(data = ResidsInfo[ResidsInfo$SIGEnv != "NS" & ResidsInfo$Environment == "FLC",], aes(x = start/1000000, y = Fitness, color = SIGStrain), alpha = 0.7, size = 2) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "white")) +
  xlab("Chromosomal Position (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("mediumorchid","goldenrod3","cyan3","forestgreen","grey"), labels = c("P75016","P75063","AMS5192","L26","NS")) +
  labs(color = "Background") +
  ylim(-1,1)
ggsave("FLCALLPOSzoom.tiff", path = fig.dir, plot = last_plot(), width = 7.5, height = 3)

ggplot(data = ResidsInfo[ResidsInfo$Environment == "MCF",], aes(x = start/1000000, y = Fitness, label = Gene)) +
  geom_point(aes(color = SIGStrain), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness-SD, ymax = Fitness+SD, color = SIGStrain), alpha = 0.7) +
  geom_point(data = ResidsInfo[ResidsInfo$SIGEnv != "NS" & ResidsInfo$Environment == "MCF",], aes(x = start/1000000, y = Fitness, color = SIGStrain), alpha = 0.7, size = 2) +
  #geom_text() +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "white")) +
  xlab("Chromosomal Position (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("mediumorchid","goldenrod3","cyan3","forestgreen","grey"), labels = c("P75016","P75063","AMS5192","L26","NS")) +
  labs(color = "Background")
ggsave("MCFALLPOS.tiff", path = fig.dir, plot = last_plot(), width = 7.5, height = 3)

ggplot(data = ResidsInfo[ResidsInfo$Environment == "SDS",], aes(x = start/1000000, y = Fitness, label = Gene)) +
  geom_point(aes(color = SIGStrain), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness-SD, ymax = Fitness+SD, color = SIGStrain), alpha = 0.7) +
  geom_point(data = ResidsInfo[ResidsInfo$SIGEnv != "NS" & ResidsInfo$Environment == "SDS",], aes(x = start/1000000, y = Fitness, color = SIGStrain), alpha = 0.7, size = 2) +
  #geom_text() +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_rect(fill = "white")) +
  xlab("Chromosomal Position (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("mediumorchid","goldenrod3","cyan3","forestgreen","grey"), labels = c("P75016","P75063","AMS5192","L26","NS")) +
  labs(color = "Background") +
  ylim(-1,1)
ggsave("SDSALLPOS.tiff", path = fig.dir, plot = last_plot(), width = 7.5, height = 3)

#Plotting a few genes of interest across environmental conditions
Justforthis <- ResidsInfo
Justforthis$Environment <- factor(Justforthis$Environment, levels = c("MCF","FLC",'SDS'))
ggplot(data = Justforthis[Justforthis$Environment %in% c("FLC","MCF","SDS") & Justforthis$Gene %in% c("C303460CA","C114140CA","C305650WA"),], aes(x = Environment, y = Fitness, group = interaction(Strain,Gene))) +
  geom_point(aes(color = Gene)) +
  geom_linerange(aes(ymin = Fitness - SD, ymax = Fitness + SD, color = Gene)) +
  geom_line(aes(color = Gene)) +
  scale_color_manual(values = c("firebrick","dodgerblue","green3"), labels = c("SET3", "C3_03460C","YCK2")) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Environment") +
  ylab("Fitness") +
  geom_hline(yintercept = 0)
ggsave("GENETRACK.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

# Plotting fitness for one background in CNV regions and colored according to environment

# The region that AMS5778 amplifies
ggplot(data = ResidsInfo[ResidsInfo$Strain == "BF"& ResidsInfo$Chrom == "Chr4"& ResidsInfo$start >= 532000 & ResidsInfo$start <= 703000,], aes(x = start/1000000, y = Fitness)) +
  geom_point(aes(color = SIGEnv), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness - SD, ymax = Fitness + SD, color = SIGEnv), alpha = 0.7) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Chromosome 4 (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("purple4","forestgreen","grey"))
ggsave("AMS5778ENVPos.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

ggplot(data = ResidsInfo[ResidsInfo$Strain == "63" & ResidsInfo$Chrom == "Chr3R",], aes(x = start/1000000, y = Fitness)) + 
  geom_point(aes(color = SIGEnv), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness - SD, ymax = Fitness + SD, color = SIGEnv), alpha = 0.7) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Chromosome 3 (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","grey","cyan3",'violetred',"goldenrod3"))
ggsave("P633RENVPos.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

ggplot(data = ResidsInfo[ResidsInfo$Strain == "63" & ResidsInfo$Chrom == "Chr1",], aes(x = start/1000000, y = Fitness)) + 
  geom_point(aes(color = SIGEnv), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness - SD, ymax = Fitness + SD, color = SIGEnv), alpha = 0.7) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black"))  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Chromosome 1 (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","grey","goldenrod3"))
ggsave("P631ENVPos.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

ggplot(data = ResidsInfo[ResidsInfo$Strain == "16" & ResidsInfo$Chrom == "Chr3R",], aes(x = start/1000000, y = Fitness)) + 
  geom_point(aes(color = SIGEnv), alpha = 0.7, size = 2) +
  geom_linerange(aes(ymin = Fitness - SD, ymax = Fitness + SD, color = SIGEnv), alpha = 0.7) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Chromosome 3 (Mb)") +
  ylab("Fitness") +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","grey",'violetred',"goldenrod3"))
ggsave("P163RENVPos.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

##################### 
# Comparison to CNV growth curve data

ResidsInfo$gRNA <- ResidsInfo$Gene

# Calculating the sums of all significant genes within a region in each environment
SigSums7083 <- ResidsInfo %>%
  filter(Chrom == "Chr1" & start >= 2115000 & start <= 3125000 & Strain == "63" & SIGEnv != "NS") %>%
  group_by(Environment) %>%
  summarise(FitSum = sum(Fitness))
SigSums7083$Strain <- rep("AMS7083", nrow(SigSums7083))

SigSums5778 <- ResidsInfo %>%
  filter(Chrom == "Chr4" & start >= 532000 & start <= 703000 & Strain == "BF" & SIGEnv != "NS") %>%
  group_by(Environment) %>%
  summarise(FitSum = sum(Fitness))
SigSums5778$Strain <- rep("AMS5778", nrow(SigSums5778))

SigSums4397 <- ResidsInfo %>%
  filter(Chrom == "Chr3R" & start >= 1050000 & start <= 1455000 & Strain == "16" & SIGEnv != "NS") %>%
  group_by(Environment) %>%
  summarise(FitSum = sum(Fitness))
SigSums4397$Strain <- rep("AMS4397", nrow(SigSums4397))

SigSums7084 <- ResidsInfo %>%
  filter(Chrom == "Chr3R" & start >= 1115000 & start <= 1451000 & Strain == "63" & SIGEnv != "NS") %>%
  group_by(Environment) %>%
  summarise(FitSum = sum(Fitness))
SigSums7084$Strain <- rep("AMS7084", nrow(SigSums7084))

SigSums <- rbind(SigSums4397, SigSums5778, SigSums7083, SigSums7084)
colnames(SigSums)[1] <- "Condition"

# Read in the growth curve data for the 'clean' CNV strains, with no other CNV or aneuploidy
CleanNorm <- read.table(paste0(output.dir,"/CNVGC.txt"), sep = "\t", header = 1)

# Joining the two data frames
TOGjoin <- full_join(SigSums, CleanNorm, by = join_by(Strain, Condition), relationship = 'many-to-many')

# Adding in zeros for the environments where no genes are significantly different. 
TOGjoin$FitSum <- ifelse(is.na(TOGjoin$FitSum), 0, TOGjoin$FitSum)

# Removing the baseline condition
TOGjoinlim <- TOGjoin[TOGjoin$Condition != "YPAD",]
colnames(TOGjoinlim)[19] <- "Relative_Growth"

#AMS7084
ggplot(data = TOGjoinlim[TOGjoinlim$Strain == "AMS7084",], aes(x = FitSum, y = Relative_Growth)) +
  geom_smooth(method = "lm") +
  geom_point(aes(x = FitSum, y = Relative_Growth, color = Condition), size = 2) +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","cyan3",'violetred',"goldenrod3","grey")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab('Growth of CNV Strain \n(AUC Relative to Euploid)') +
  xlab('Sum of Individual Gene OE Fitness')
ggsave("CNVCRISPRComp7084.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

#Getting the linear model
lmmod <- lm(TOGjoinlim[TOGjoinlim$Strain == "AMS7084",]$Relative_Growth ~ TOGjoinlim[TOGjoinlim$Strain == "AMS7084",]$FitSum)
sum <- summary(lmmod)
sum$coefficients #To get actual p-value

# AMS7083
ggplot(data = TOGjoinlim[TOGjoinlim$Strain == "AMS7083",], aes(x = FitSum, y = Relative_Growth)) +
  #geom_rect(aes(xmin=bottom_2.5_percent_value, xmax=top_2.5_percent_value, ymin=-Inf, ymax=Inf), fill = "grey", alpha = 0.01) +
  geom_smooth(method = "lm") +
  geom_point(aes(x = FitSum, y = Relative_Growth, color = Condition), size = 2) +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","cyan3",'violetred',"goldenrod3","grey")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab('Growth of CNV Strain \n(AUC Relative to Euploid)') +
  xlab('Sum of Individual Gene OE Fitness')
ggsave("CNVCRISPRComp7083.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

lmmod <- lm(TOGjoinlim[TOGjoinlim$Strain == "AMS7083",]$Relative_Growth ~ TOGjoinlim[TOGjoinlim$Strain == "AMS7083",]$FitSum)
sum <- summary(lmmod)
sum$coefficients #To get actual p-value

# AMS4397
ggplot(data = TOGjoinlim[TOGjoinlim$Strain == "AMS4397",], aes(x = FitSum, y = Relative_Growth)) + #Two points are so close they look like one for FLC
  #geom_rect(aes(xmin=bottom_2.5_percent_value, xmax=top_2.5_percent_value, ymin=-Inf, ymax=Inf), fill = "grey", alpha = 0.01) +
  geom_smooth(method = "lm") +
  geom_point(aes(x = FitSum, y = Relative_Growth, color = Condition), size = 2) +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","cyan3",'violetred',"goldenrod3","grey")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab('Growth of CNV Strain \n(AUC Relative to Euploid)') +
  xlab('Sum of Individual Gene OE Fitness')
ggsave("CNVCRISPRComp4397.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

#Getting the linear model
lmmod <- lm(TOGjoinlim[TOGjoinlim$Strain == "AMS4397",]$Relative_Growth ~ TOGjoinlim[TOGjoinlim$Strain == "AMS4397",]$FitSum)
sum <- summary(lmmod)
sum$coefficients #To get actual p-value

ggplot(data = TOGjoinlim[TOGjoinlim$Strain == "AMS5778",], aes(x = FitSum, y = Relative_Growth)) +
  #geom_rect(aes(xmin=bottom_2.5_percent_value, xmax=top_2.5_percent_value, ymin=-Inf, ymax=Inf), fill = "grey", alpha = 0.01) +
  geom_smooth(method = "lm") +
  geom_point(aes(x = FitSum, y = Relative_Growth, color = Condition), size = 2) +
  scale_color_manual(values = c("firebrick","darkorange","purple4","forestgreen","cyan3",'violetred',"goldenrod3","grey")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  ylab('Growth of CNV Strain \n(AUC Relative to Euploid)') +
  xlab('Sum of Individual Gene OE Fitness')
ggsave("CNVCRISPRComp5778.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3)

#Getting the linear model
lmmod <- lm(TOGjoinlim[TOGjoinlim$Strain == "AMS5778",]$Relative_Growth ~ TOGjoinlim[TOGjoinlim$Strain == "AMS5778",]$FitSum)
sum <- summary(lmmod)
sum$coefficients #To get actual p-value

########################
## Predicting fitness effects of 'simulated CNVs'

#Adding on positional information to the CompmultirepsMean dataframe
CompmultirepsMean <- left_join(CompmultirepsMean, ResidsInfo[,c("gRNA","Chrom","start")], by = "gRNA")

#First split up the dataframe by chromosome to make paneling across each easier.
Chr1Genes <- CompmultirepsMean[CompmultirepsMean$Chrom == "Chr1",c(grep("Resid", colnames(CompmultirepsMean)),1,159)]
Chr1Genes <- as.data.frame(Chr1Genes)
# Replacing nonsignificant genes with 0 - This only works because these vectors are in the same order as the column names!
n <- 1
for (i in 1:length(Strainvec)) {
  for (j in 1:length(Conditionvec)) {
    if (length(ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene) == 0) {
      Chr1Genes[,n] <- rep(0, nrow(Chr1Genes))
    } else {
      Chr1Genes[,n] <- ifelse(Chr1Genes$gRNA %in% ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene, Chr1Genes[,n], 0)
    }
    n <- n +1
  }
}


Chr3LGenes <- CompmultirepsMean[CompmultirepsMean$Chrom == "Chr3L",c(grep("Resid", colnames(CompmultirepsMean)),1,159)]
Chr3LGenes <- as.data.frame(Chr3LGenes)
n <- 1
for (i in 1:length(Strainvec)) {
  for (j in 1:length(Conditionvec)) {
    if (length(ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene) == 0) {
      Chr3LGenes[,n] <- rep(0, nrow(Chr3LGenes))
    } else {
      Chr3LGenes[,n] <- ifelse(Chr3LGenes$gRNA %in% ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene, Chr3LGenes[,n], 0)
    }
    n <- n +1
  }
}

Chr3RGenes <- CompmultirepsMean[CompmultirepsMean$Chrom == "Chr3R",c(grep("Resid", colnames(CompmultirepsMean)),1,159)]
Chr3RGenes <- as.data.frame(Chr3RGenes)
n <- 1
for (i in 1:length(Strainvec)) {
  for (j in 1:length(Conditionvec)) {
    if (length(ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene) == 0) {
      Chr3RGenes[,n] <- rep(0, nrow(Chr3RGenes))
    } else {
      Chr3RGenes[,n] <- ifelse(Chr3RGenes$gRNA %in% ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene, Chr3RGenes[,n], 0)
    }
    n <- n +1
  }
}

Chr4Genes <- CompmultirepsMean[CompmultirepsMean$Chrom == "Chr4",c(grep("Resid", colnames(CompmultirepsMean)),1,159)]
Chr4Genes <- as.data.frame(Chr4Genes)
n <- 1
for (i in 1:length(Strainvec)) {
  for (j in 1:length(Conditionvec)) {
    if (length(ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene) == 0) {
      Chr4Genes[,n] <- rep(0, nrow(Chr4Genes))
    } else {
      Chr4Genes[,n] <- ifelse(Chr4Genes$gRNA %in% ResidsInfo[ResidsInfo$SIGStrain == Strainvec[i] & ResidsInfo$SIGEnv == Conditionvec[j],]$Gene, Chr4Genes[,n], 0)
    }
    n <- n +1
  }
}

SLIDINGFIT <- function(x, wind, shift) {
  n <- min(x$start)
  maxn <- max(x$start)
  Rows <- which(x$start >= n & x$start <= n+wind)
  WindSums <- as.data.frame(t(colSums(x[Rows,grep("Resid", colnames(x))])))
  n <- n+shift
  while(n < maxn) {
    Rows <- which(x$start >= n & x$start <= n+wind)
    newrow <- colSums(x[Rows,grep("Resid", colnames(x))])
    WindSums <- rbind(WindSums, newrow)
    n <- n+shift
  }
  return(WindSums)
}

# A function to apply this across different window sizes for all chromosomes
Windows <- c(50000, 100000, 200000, 300000, 400000, 500000) #From 50 to 500KB
WINDOW <- function(x) {
  Chr1Slide <- SLIDINGFIT(Chr1Genes, wind = x, shift = 1000)
  Chr1Slide$ID <- c(1:nrow(Chr1Slide))
  Chr1Slidemelt <- pivot_longer(Chr1Slide, cols = colnames(Chr1Slide)[1:ncol(Chr1Slide)-1])
  Chr1Slidemelt$Strain <- sapply(strsplit(Chr1Slidemelt$name, "_"), "[", 1)
  Chr1Slidemelt$Environment <- sapply(strsplit(Chr1Slidemelt$name, "_"), "[", 2)
  Chr1Slidemelt$Chrom <- "Chr1"
  Chr3RSlide <- SLIDINGFIT(Chr3RGenes, wind = x, shift = 1000) 
  Chr3RSlide$ID <- c(1:nrow(Chr3RSlide))
  Chr3RSlidemelt <- pivot_longer(Chr3RSlide, cols = colnames(Chr3RSlide)[1:ncol(Chr3RSlide)-1])
  Chr3RSlidemelt$Strain <- sapply(strsplit(Chr3RSlidemelt$name, "_"), "[", 1)
  Chr3RSlidemelt$Environment <- sapply(strsplit(Chr3RSlidemelt$name, "_"), "[", 2)
  Chr3RSlidemelt$Chrom <- "Chr3R"
  Chr3LSlide <- SLIDINGFIT(Chr3LGenes, wind = x, shift = 1000) 
  Chr3LSlide$ID <- c(1:nrow(Chr3LSlide))
  Chr3LSlidemelt <- pivot_longer(Chr3LSlide, cols = colnames(Chr3LSlide)[1:ncol(Chr3LSlide)-1])
  Chr3LSlidemelt$Strain <- sapply(strsplit(Chr3LSlidemelt$name, "_"), "[", 1)
  Chr3LSlidemelt$Environment <- sapply(strsplit(Chr3LSlidemelt$name, "_"), "[", 2)
  Chr3LSlidemelt$Chrom <- "Chr3L"
  Chr4Slide <- SLIDINGFIT(Chr4Genes, wind = x, shift = 1000) 
  Chr4Slide$ID <- c(1:nrow(Chr4Slide))
  Chr4Slidemelt <- pivot_longer(Chr4Slide, cols = colnames(Chr4Slide)[1:ncol(Chr4Slide)-1])
  Chr4Slidemelt$Strain <- sapply(strsplit(Chr4Slidemelt$name, "_"), "[", 1)
  Chr4Slidemelt$Environment <- sapply(strsplit(Chr4Slidemelt$name, "_"), "[", 2)
  Chr4Slidemelt$Chrom <- "Chr4"
  TOGSlide <- rbind(Chr1Slidemelt, Chr3RSlidemelt, Chr3LSlidemelt, Chr4Slidemelt)
  return(TOGSlide)
}
Windowsize <- lapply(Windows, WINDOW)

#Now plotting each of these
a <- ggplot(data = Windowsize[[1]][Windowsize[[1]]$Environment %in% c("FLC") & Windowsize[[1]]$Strain == "63",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[1]][Windowsize[[1]]$Environment %in% c("MCF") & Windowsize[[1]]$Strain == "63",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")

b <- ggplot(data = Windowsize[[2]][Windowsize[[2]]$Environment %in% c("FLC") & Windowsize[[2]]$Strain == "63",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[2]][Windowsize[[2]]$Environment %in% c("MCF") & Windowsize[[2]]$Strain == "63",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")

c <- ggplot(data = Windowsize[[3]][Windowsize[[3]]$Environment %in% c("FLC") & Windowsize[[3]]$Strain == "63",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[3]][Windowsize[[3]]$Environment %in% c("MCF") & Windowsize[[3]]$Strain == "63",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")

d <- ggplot(data = Windowsize[[4]][Windowsize[[4]]$Environment %in% c("FLC") & Windowsize[[4]]$Strain == "63",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[4]][Windowsize[[4]]$Environment %in% c("MCF") & Windowsize[[4]]$Strain == "63",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")

e <- ggplot(data = Windowsize[[5]][Windowsize[[5]]$Environment %in% c("FLC") & Windowsize[[5]]$Strain == "63",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[5]][Windowsize[[5]]$Environment %in% c("MCF") & Windowsize[[5]]$Strain == "63",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")

f <- ggplot(data = Windowsize[[6]][Windowsize[[6]]$Environment %in% c("FLC") & Windowsize[[6]]$Strain == "63",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[6]][Windowsize[[6]]$Environment %in% c("MCF") & Windowsize[[6]]$Strain == "63",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  facet_grid(~Chrom, scales = "free_x", space = "free_x") +
  theme_bw() + theme(axis.text = element_text(size = 12, angle = 45, hjust = 1, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")

plot_grid(a, b, c, d, e, f, nrow = 6)
ggsave("SizeComp.tiff", path = fig.dir, plot = last_plot(), width = 7, height = 15)

plot_grid(a, b, d, f, nrow = 4)
ggsave("SizeComplim.tiff", path = fig.dir, plot = last_plot(), width = 7, height = 10)

#A few individual examples
ggplot(data = Windowsize[[1]][Windowsize[[1]]$Environment %in% c("FLC") & Windowsize[[1]]$Strain == "63" & Windowsize[[1]]$Chrom == "Chr3R",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[1]][Windowsize[[1]]$Environment %in% c("MCF") & Windowsize[[1]]$Strain == "63"& Windowsize[[1]]$Chrom == "Chr3R",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")
ggsave("Chr3RSlide50K.tiff", path = fig.dir, width = 4, height = 3, plot = last_plot())

ggplot(data = Windowsize[[6]][Windowsize[[6]]$Environment %in% c("FLC") & Windowsize[[6]]$Strain == "63"& Windowsize[[6]]$Chrom == "Chr3R",], aes(x = ID, y = value)) +
  geom_col(fill = "purple4", alpha = 0.7) +
  geom_col(data = Windowsize[[6]][Windowsize[[6]]$Environment %in% c("MCF") & Windowsize[[6]]$Strain == "63"& Windowsize[[6]]$Chrom == "Chr3R",], aes(x = ID, y = value), fill = "forestgreen", alpha = 0.5) +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  xlab("Window Start Position") +
  ylab("Additive Fitness")
ggsave("Chr3RSlide500K.tiff", path = fig.dir, width = 4, height = 3, plot = last_plot())


#Now, how might these be constrained by the positions of long inverted repeats?
#Uploading the data from Supplementary File 2 from Todd et al 2019
LIRS <- read.xlsx(paste0(data.dir,"elife-45954-supp2-v2.xlsx"))

# I just want things that are both on the same chromosome, where they would recombine with each other to amplify that region
# Also, I just need chromosomes 1, 3, and 4 for this, and only where both start and stop are in the same chromosome
LIRSlim <- LIRS[LIRS$Chromosome.Sequence.1 == LIRS$Chromosome.Sequence.2 & LIRS$Chromosome.Sequence.1 %in% c(1,3,4),]

Chr1LIRS <- LIRSlim[LIRSlim$Chromosome.Sequence.1 == 1 & (LIRSlim$Start.Sequence.1 >= min(Chr1Genes$start) & LIRSlim$Start.Sequence.1 <= max(Chr1Genes$start)) & (LIRSlim$Start.Sequence.2 >= min(Chr1Genes$start) & LIRSlim$Start.Sequence.2 <= max(Chr1Genes$start)),]
CNVsFit <- data.frame(Chr1Genes[1,grep("Resid",colnames(Chr1Genes))]) #Just to set up the right colnames
for(i in 1:nrow(Chr1LIRS)) {
  Rows <- which((Chr1Genes$start >= Chr1LIRS[i,"Start.Sequence.1"] & Chr1Genes$start <= Chr1LIRS[i,"Start.Sequence.2"]) | (Chr1Genes$start <= Chr1LIRS[i,"Start.Sequence.1"] & Chr1Genes$start >= Chr1LIRS[i,"Start.Sequence.2"]))
  Fitsums <- colSums(Chr1Genes[Rows,grep("Resid", colnames(Chr1Genes))])
  CNVsFit <- rbind(CNVsFit, Fitsums)
}
Chr1LIRS <- cbind(Chr1LIRS, CNVsFit[2:nrow(CNVsFit),]) #Getting rid of the first row only used to set up the dataframe.

Chr3LLIRS <- LIRSlim[LIRSlim$Chromosome.Sequence.1 == 3 & (LIRSlim$Start.Sequence.1 >= min(Chr3LGenes$start) & LIRSlim$Start.Sequence.1 <= max(Chr3LGenes$start)) & (LIRSlim$Start.Sequence.2 >= min(Chr3LGenes$start) & LIRSlim$Start.Sequence.2 <= max(Chr3LGenes$start)),]
CNVsFit <- data.frame(Chr3LGenes[1,grep("Resid",colnames(Chr3LGenes))]) #Just to set up the right colnames
for(i in 1:nrow(Chr3LLIRS)) {
  Rows <- which((Chr3LGenes$start >= Chr3LLIRS[i,"Start.Sequence.1"] & Chr3LGenes$start <= Chr3LLIRS[i,"Start.Sequence.2"]) | (Chr3LGenes$start <= Chr3LLIRS[i,"Start.Sequence.1"] & Chr3LGenes$start >= Chr3LLIRS[i,"Start.Sequence.2"]))
  Fitsums <- colSums(Chr3LGenes[Rows,grep("Resid", colnames(Chr3LGenes))])
  CNVsFit <- rbind(CNVsFit, Fitsums)
}
Chr3LLIRS <- cbind(Chr3LLIRS, CNVsFit[2:nrow(CNVsFit),])

Chr3RLIRS <- LIRSlim[LIRSlim$Chromosome.Sequence.1 == 3 & (LIRSlim$Start.Sequence.1 >= min(Chr3RGenes$start) & LIRSlim$Start.Sequence.1 <= max(Chr3RGenes$start)) & (LIRSlim$Start.Sequence.2 >= min(Chr3RGenes$start) & LIRSlim$Start.Sequence.2 <= max(Chr3RGenes$start)),]
CNVsFit <- data.frame(Chr3RGenes[1,grep("Resid",colnames(Chr3RGenes))]) #Just to set up the right colnames
for(i in 1:nrow(Chr3RLIRS)) {
  Rows <- which((Chr3RGenes$start >= Chr3RLIRS[i,"Start.Sequence.1"] & Chr3RGenes$start <= Chr3RLIRS[i,"Start.Sequence.2"]) | (Chr3RGenes$start <= Chr3RLIRS[i,"Start.Sequence.1"] & Chr3RGenes$start >= Chr3RLIRS[i,"Start.Sequence.2"]))
  Fitsums <- colSums(Chr3RGenes[Rows,grep("Resid", colnames(Chr3RGenes))])
  CNVsFit <- rbind(CNVsFit, Fitsums)
}
Chr3RLIRS <- cbind(Chr3RLIRS, CNVsFit[2:nrow(CNVsFit),])

Chr4LIRS <- LIRSlim[LIRSlim$Chromosome.Sequence.1 == 4 & (LIRSlim$Start.Sequence.1 >= min(Chr4Genes$start) & LIRSlim$Start.Sequence.1 <= max(Chr4Genes$start)) & (LIRSlim$Start.Sequence.2 >= min(Chr4Genes$start) & LIRSlim$Start.Sequence.2 <= max(Chr4Genes$start)),]
CNVsFit <- data.frame(Chr4Genes[1,grep("Resid",colnames(Chr4Genes))]) #Just to set up the right colnames
for(i in 1:nrow(Chr4LIRS)) {
  Rows <- which((Chr4Genes$start >= Chr4LIRS[i,"Start.Sequence.1"] & Chr4Genes$start <= Chr4LIRS[i,"Start.Sequence.2"]) | (Chr4Genes$start <= Chr4LIRS[i,"Start.Sequence.1"] & Chr4Genes$start >= Chr4LIRS[i,"Start.Sequence.2"]))
  Fitsums <- colSums(Chr4Genes[Rows,grep("Resid", colnames(Chr4Genes))])
  CNVsFit <- rbind(CNVsFit, Fitsums)
}
Chr4LIRS <- cbind(Chr4LIRS, CNVsFit[2:nrow(CNVsFit),])

ALLLIRS <- rbind(Chr1LIRS, Chr3LLIRS, Chr3RLIRS,Chr4LIRS)

ggplot(data = ALLLIRS, aes(x = `BF_FLC_Resid`, y = `BF_MCF_Resid`)) +
  geom_point(aes(shape = Chromosome.Sequence.1), size = 3, color = "cyan3") +
  geom_point(aes(x = `X63_FLC_Resid`, y = `X63_MCF_Resid`, shape = Chromosome.Sequence.1), size = 3, color = "goldenrod3") +
  geom_point(aes(x = `X16_FLC_Resid`, y = `X16_MCF_Resid`, shape = Chromosome.Sequence.1), size = 3, color = "mediumorchid") +
  geom_point(aes(x = `L26_FLC_Resid`, y = `L26_MCF_Resid`, shape = Chromosome.Sequence.1), size = 3, color = "forestgreen") +
  theme_bw() + theme(axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12, color = "black"), plot.title = element_text(size = 15, hjust = 0.5, color = "black")) +
  labs(shape = "Chr") +
  xlab("Fitness in Fluconazole") +
  ylab("Fitness in Micafungin") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
ggsave("RepeatCNVSim.tiff", path = fig.dir, plot = last_plot(), width = 4, height = 3.5)
