# Load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# Load reference Seurat object
immune_ref <- LoadH5Seurat("path/to/CCA_integrate_only_dengue_gts_ss2_10x.h5Seurat")

# Load 10X raw count data
raw_counts <- Read10X(data.dir = "path/to/matrix_files_new_10x_CD8_rawcount")
cd8_10x <- CreateSeuratObject(counts = raw_counts, project = "CD8_10X")

# Add metadata
metadata <- read.csv("path/to/matrix_files_new_10x_CD8_rawcount_metadata.csv", row.names = "index")
cd8_10x <- AddMetaData(cd8_10x, metadata = metadata)

# Normalize
cd8_10x <- NormalizeData(cd8_10x)

# Split by donor
donor_ids <- unique(cd8_10x$donor_id)
query_list <- lapply(donor_ids, function(id) subset(cd8_10x, subset = donor_id == id))
names(query_list) <- donor_ids

# Transfer annotations for each donor
annotated_list <- lapply(query_list, function(query) {
  anchors <- FindTransferAnchors(
    reference = immune_ref, query = query,
    reference.reduction = "pca", dims = 1:30, k.anchor = 6
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = immune_ref$seurat_clusters,
    dims = 1:30, k.weight = 40
  )
  AddMetaData(query, metadata = predictions)
})
names(annotated_list) <- donor_ids

# Merge all annotated data
annotated_merged <- Reduce(function(x, y) merge(x, y), annotated_list)

# Set severity factor level
annotated_merged$severity <- factor(annotated_merged$severity, levels = c("AD", "DF", "DHF"))

# Plot predicted cell types per donor
p <- ggplot(annotated_merged@meta.data, aes(x = donor_id, fill = predicted.id)) +
  geom_bar(position = "fill") +
  facet_wrap(~severity, scales = "free_x") +
  theme_minimal() +
  labs(y = "Proportion", title = "Predicted Cell Types by Donor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)

# Export metadata
meta_export <- data.frame(
  cell_id = rownames(annotated_merged@meta.data),
  predicted_id = annotated_merged$predicted.id,
  donor_id = annotated_merged$donor_id,
  severity = annotated_merged$severity
)
write.csv(meta_export, file = "10x_asymp_donor_map_gts_ka6_kw40.csv", row.names = FALSE)
