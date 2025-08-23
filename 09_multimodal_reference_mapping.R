# Load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# Load reference Seurat object
immune.combine <- readRDS("path/to/CCA_integrate_only_dengue_gts_ss2_10x.rds")

# Load data and metadata 
cd8t_10x.data <- Read10X(
  data.dir = "path/matrix_files_new_10x_CD8_rawcount"
)

cd8t_10x <- CreateSeuratObject(counts = cd8t_10x.data, project = "cd8t_10x")

metadata <- read.csv(
  "path/matrix_files_new_10x_CD8_rawcount_metadata_annotation.csv",
  row.names = "index"
)
cd8t_10x <- AddMetaData(cd8t_10x, metadata = metadata)

# Prepare query data 
data.query <- NormalizeData(cd8t_10x)

# Split by donor 
donor_ids <- unique(data.query@meta.data$donor_id)
query_list <- lapply(donor_ids, function(donor) {
  subset(data.query, subset = donor_id == donor)
})
names(query_list) <- donor_ids

# Annotate cell type by donor 
annotated_list <- lapply(query_list, function(query_obj) {
  anchors <- FindTransferAnchors(
    reference = immune.combine,
    query = query_obj,
    dims = 1:30,
    reference.reduction = "pca",
    k.anchor = 6
  )
  preds <- TransferData(
    anchorset = anchors,
    refdata = immune.combine$seurat_clusters,
    dims = 1:30,
    k.weight = 40
  )
  AddMetaData(query_obj, metadata = preds)
})
names(annotated_list) <- donor_ids

#  Merge annotated objects 
data.query.merged <- merge(annotated_list[[1]], y = annotated_list[-1])

# Set factor levels 
data.query.merged@meta.data$severity <- factor(
  data.query.merged@meta.data$severity,
  levels = c("AD", "DF", "DHF")
)

# Visualize cell proportion across donor
ggplot(data.query.merged@meta.data, aes(x = donor_id, fill = predicted.id)) +
  geom_bar(position = "fill") +
  facet_wrap(~ severity, scales = "free_x") +
  theme_minimal() +
  labs(y = "Proportion", title = "Predicted Cell Types by Donor") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Extract metadata
metadata_df <- data.query.merged@meta.data
meta_export <- data.frame(
  cell_id = rownames(metadata_df),
  predicted_id = metadata_df$predicted.id,
  donor_id = metadata_df$donor_id,
  severity = metadata_df$severity
)

# Write to CSV
write.csv(meta_export, file = "path/10x_asymp_donor_map_gts_ka6_kw40_new_5_cluster.csv", row.names = FALSE)

