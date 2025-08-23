# Load libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# Load 10x data and preprocessing
gts_10x.data <- Read10X("path/matrix_files_10x_gts_rawcount")
metadata_10x      <- read.csv("path/matrix_files_10x_gts_rawcount_metadata_pair_annotation.csv", row.names = "index")
gts_10x <- CreateSeuratObject(gts_10x.data, project = "gts_10x", meta.data = metadata_10x, min.cells = 3, min.features = 200)
gts_10x.list <- SplitObject(gts_10x, split.by = "donor_id")
gts_10x.list <- lapply(gts_10x.list, function(x) { x <- NormalizeData(x); FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) })
features <- SelectIntegrationFeatures(gts_10x.list)
gts_10x.list <- lapply(gts_10x.list, function(x) { x <- ScaleData(x, features = features, verbose = FALSE); RunPCA(x, features = features, verbose = FALSE) })
gts_10x_merge <- Reduce(function(x, y) merge(x, y, add.cell.ids = NULL, merge.data = TRUE), gts_10x.list)
gts_10x_merge_subset <- subset(gts_10x_merge, subset = epitope == "gts")

# Load SS2 data and preprocessing
metadata_ss2   <- read.csv("path/matrix_files_ss2_gts_rawcount_metadata_annotation.csv", row.names = "X")
gts_ss2.data <- Read10X("path/matrix_files_ss2_gts_rawcount")
gts_ss2 <- CreateSeuratObject(gts_ss2.data, project = "gts_ss2", meta.data = metadata_ss2, min.cells = 3, min.features = 200)
gts_ss2.list <- SplitObject(gts_ss2, split.by = "plate_id")
gts_ss2.list <- lapply(gts_ss2.list, function(x) { x <- NormalizeData(x); FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) })
gts_ss2_merge <- Reduce(function(x, y) merge(x, y, add.cell.ids = NULL, merge.data = TRUE), gts_ss2.list)

# Integrate data
full.list <- c(gts_10x_merge_subset, gts_ss2_merge)
features  <- SelectIntegrationFeatures(full.list)
full.list <- lapply(full.list, function(x) { x <- ScaleData(x, features = features, verbose = FALSE); RunPCA(x, features = features, verbose = FALSE, npcs = 50) })
immune.anchors  <- FindIntegrationAnchors(object.list = full.list, anchor.features = features, reduction = "cca", k.anchor = 4, k.filter = 200)
immune.combined <- IntegrateData(anchorset = immune.anchors, k.weight = 200)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:50)
immune.combined <- FindClusters(immune.combined, resolution = 0.7)

## Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "donor_id")
p2 <- DimPlot(immune.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
p1+p2

## Find Markers 
markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

## Heatmap of top10 per cluster 
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(immune.combined, features = unique(top10$gene)) + NoLegend()

## Density plots for selected genes
genes <- c("CX3CR1","CXCR6","CREM","IL7R","MKI67")
invisible(lapply(genes, function(g) print(plot_density(immune.combined, g))))
