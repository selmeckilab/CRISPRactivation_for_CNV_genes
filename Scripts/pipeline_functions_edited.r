############
# This is the source file containing functions for use with 'AnalysisandFigs.R'
#
# These functions are edited from Martinson et al, 2023. There they note:
# Many of these functions are re-writes of code from Wetmore (et al. 2015)
# pipeline or the Morin (et al. 2019) pipeline.
#
############

library(tidyr)

get_barcodes_with_any_below_threshold = function(counts, T0s, threshold){
  # this returns the list of barcodes where any of the T0 samples
  # have a count below the threshold
  counts %>%
    select(barcode, T0s) %>%
    pivot_longer(-barcode, names_to = "rep", values_to = "count") %>% 
    group_by(barcode) %>%
    summarize(below_threshold = any(count < threshold)) %>%
    filter(below_threshold) %>% 
    select(barcode) %>%
    unlist() %>% unname
}

get_replicates_with_low_median_gene_reads = function(counts, threshold){
  test <- counts %>%
    pivot_longer(-barcode, names_to = "replicate", values_to = "count") %>% 
    group_by(replicate) %>%
    summarize(median_reads_per_barcode = median(count)) %>%
    ungroup() %>%
    filter(median_reads_per_barcode < threshold) %>%
    select(replicate) %>%
    unname() %>% unlist()
}

add_pseudocount = function(counts, pseudocount = 0.1){
  counts %>%
    pivot_longer(cols = c(-barcode),
                 names_to = "replicate", values_to = "count") %>%
    mutate(raw_count = count, 
           count = count + pseudocount) 
}

filter_barcodes = function(counts, T0_thresh = 0, f_min = 0, f_max = 1){
  counts %>%
    filter(T0 > T0_thresh & f > f_min & f < f_max)
}

normalize_using_reference_genes = function(counts, ref){
  counts %>%
    left_join(counts %>%
    filter(barcode %in% ref) %>%
    mutate(overall_mean_ref_count = mean(raw_count, na.rm = TRUE)) %>%
    group_by(replicate) %>%
    summarize(norm_factor = first(overall_mean_ref_count) / mean(raw_count),
              mean_ref_count = mean(raw_count),
              overall_mean_ref_count = first(overall_mean_ref_count))) %>%
    mutate(count_normalized = count * norm_factor)
}

calculate_strain_fitness_from_T0 = function(counts, T0_name = "T0"){
  # this calculates strain fitness and strain fitness variance
  # it uses the depth-normalized counts for fitness, and
  # it uses the raw counts for variance
  counts %>%
    mutate(log2count = log2(count_normalized)) %>%
    group_by(barcode) %>%
    mutate(strain_fitness = log2count - log2count[replicate == T0_name]) %>%
    mutate(strain_fitness_var = ((1 / (raw_count + 1)) + (1 / (raw_count[replicate == T0_name] + 1))) / log(2)^2) %>%
    ungroup() 
}

calculate_weighted_gene_fitness_edit = function(fitness, count_cap = 50){
  # My edited version of this
  cap = 1 / (((1/(count_cap+1))+(1/(count_cap+1)))/(log(2)^2))
  fitness %>% 
    mutate(strain_weight = 1/strain_fitness_var) %>%
    mutate(strain_weight_capped = ifelse(strain_weight > cap, 
                                         cap, strain_weight)) %>%
    group_by(replicate, gRNA) %>%
    mutate(fitness = sum(strain_weight_capped*strain_fitness)/sum(strain_weight_capped),
           fitness_left = sum(strain_weight_capped[DIST <= 300] *strain_fitness[DIST <= 300])/sum(strain_weight_capped[DIST <= 300]),
           fitness_right = sum(strain_weight_capped[DIST > 300] *strain_fitness[DIST > 300])/sum(strain_weight_capped[DIST > 300])) %>%
    ungroup() 
}

# A modification of the UpsetR 'fromList' function that retains original value names
fromList <- function (input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input)))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  row.names(data) <- elements
  return(data)
}

