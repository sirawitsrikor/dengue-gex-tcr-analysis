# Load libraries
library(dplyr)
library(ggplot2)
library(scRepertoire)

# 1) read & split
contigs <- read.csv("path/single_A11_24_10x_filtered_contig_annotation_plusmeta.csv")
individual_dfs <- split(contigs, contigs$sample)  # named list by 'sample'

# 2) build samples from donor_id and keep only those that exist
donor_ids <- unique(contigs$donor_id)
samples <- c(paste0(donor_ids, "_Acute"), paste0(donor_ids, "_Convalescence"))
samples <- samples[samples %in% names(individual_dfs)]

# **critical**: subset & reorder list to match `samples`
individual_dfs <- individual_dfs[samples]

# sanity check
stopifnot(length(individual_dfs) == length(samples))

# 3) combine
combined.TCR <- combineTCR(
  individual_dfs,
  samples     = samples,
  removeNA    = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

acute_names <- sub("_Acute$", "", names(combined.TCR)[grepl("Acute$", names(combined.TCR))])
conv_names  <- sub("_Convalescence$", "", names(combined.TCR)[grepl("Convalescence$", names(combined.TCR))])
donors <- intersect(acute_names, conv_names)

for (d in donors) {
  pair <- paste0(d, c("_Acute", "_Convalescence"))
  p <- clonalCompare(
    combined.TCR,
    top.clones = 20,
    samples    = pair,
    cloneCall  = "aa",    # change if your AA column has a different name
    graph      = "alluvial",
    palette    = "viridis"
  )
  print(p)
}
