# Load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

#---------------------------#
# Load 10x GTS Dataset
#---------------------------#
gts_10x_counts <- Read10X("path/to/matrix_files_10x_gts_rawcount")
gts_10x_meta <- read.csv("path/to/matrix_files_10x_gts_rawcount_metadata_pair.csv", row.names = "index")
gts_10x <- CreateSeuratObject(counts = gts_10x_counts, meta.data = gts_10x_meta, project = "GTS_10X", min.cells = 3, min.features = 200)
gts_10x_list <- SplitObject(gts_10x, split.by = "Donor_ID")

gts_10x_list <- lapply(gts_10x_list, function(x) {
  NormalizeData(x) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
})

#---------------------------#
# Load SS2 GTS Dataset
#---------------------------#
gts_ss2_counts <- Read10X("path/to/matrix_files_ss2_gts_rawcount")
gts_ss2_meta <- read.csv("path/to/matrix_files_ss2_gts_rawcount_metadata.csv", row.names = 1)
gts_ss2 <- CreateSeuratObject(counts = gts_ss2_counts, meta.data = gts_ss2_meta, project = "GTS_SS2", min.cells = 3, min.features = 200)
gts_ss2_list <- SplitObject(gts_ss2, split.by = "Plate_ID")

gts_ss2_list <- lapply(gts_ss2_list, function(x) {
  NormalizeData(x) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
})

#---------------------------#
# Merge and subset
#---------------------------#
gts_10x_merged <- Reduce(merge, gts_10x_list)
gts_ss2_merged <- Reduce(merge, gts_ss2_list)
gts_10x_paired <- subset(gts_10x_merged, subset = pair == "paired")

integration_list <- list(gts_10x_paired, gts_ss2_merged)

#---------------------------#
# Integration
#---------------------------#
features <- SelectIntegrationFeatures(integration_list)

integration_list <- lapply(integration_list, function(x) {
  ScaleData(x, features = features, verbose = FALSE) %>%
    RunPCA(features = features, npcs = 50, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(
  object.list = integration_list,
  anchor.features = features,
  reduction = "cca",
  k.anchor = 4,
  k.filter = 200
)

immune.combine <- IntegrateData(anchorset = anchors, k.weight = 200)

#---------------------------#
# Post-integration Analysis
#---------------------------#
DefaultAssay(immune.combine) <- "integrated"

immune.combine <- immune.combine %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters(resolution = 0.7)

#---------------------------#
# Marker Analysis
#---------------------------#
immune.combine.marker <- FindAllMarkers(
  immune.combine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)

top2 <- immune.combine.marker %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
print(top2)

#---------------------------#
# Heatmap of Top10 Markers
#---------------------------#
top10 <- immune.combine.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap_plot <- DoHeatmap(immune.combine, features = top10$gene) + NoLegend()
ggsave("figures/6sup_heatmap_small.pdf", heatmap_plot, width = 7, height = 5)
print(heatmap_plot)

#---------------------------#
# Gene Expression Density Plots
#---------------------------#
plot_genes <- c("CX3CR1", "CXCR6", "CREM", "IL7R", "MKI67", "STMN1")

for (gene in plot_genes) {
  p <- plot_density(immune.combine, gene)
  ggsave(paste0("figures/6sup_density_plot_", tolower(gene), ".pdf"), p, width = 4.5, height = 3)
  print(p)
}

#---------------------------#
# UMAP Plot
#---------------------------#
umap_plot <- DimPlot(immune.combine, reduction = "umap")
ggsave("figures/6sup_umap.pdf", umap_plot, width = 4.5, height = 3)
print(umap_plot)
