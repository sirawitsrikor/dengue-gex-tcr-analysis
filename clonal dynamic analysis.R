# Load libraries
library(dplyr)
library(ggplot2)
library(scRepertoire)

# Load filtered contig file with metadata
contigs <- read.csv("path/to/single_A11_24_10x_filtered_contig_annotation_plusmeta.csv")

# Define donors and their sample names
donor_ids <- c("T006", "I085", "I159", "I147")
sample_names <- unlist(lapply(donor_ids, function(id) c(paste0(id, "_Acute"), paste0(id, "_Convalescence"))))

# Subset contigs for selected donors
contigs_filtered <- contigs %>% filter(sample %in% sample_names)

# Split into sample-wise list for clonal analysis
combined.TCR <- split(contigs_filtered, contigs_filtered$sample)

# Loop through each donor and plot clonal comparison
for (donor in donor_ids) {
  clonalCompare(
    combined.TCR,
    top.clones = 10,
    samples = c(paste0(donor, "_Acute"), paste0(donor, "_Convalescence")),
    cloneCall = "aa",
    graph = "alluvial",
    group.by = "clinical_phase",
    palette = "viridis"
  )
}
