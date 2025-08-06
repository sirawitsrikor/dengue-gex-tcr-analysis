# Load required libraries
library(CellChat)
library(future)
library(ggplot2)
library(patchwork)

# Load metadata
metadata <- read.csv("path/to/allcelltype_10x_asymp_cd8subset_mapgts_bydonor.csv", row.names = "index")
rownames(metadata) <- gsub("-1$", "", rownames(metadata))
metadata <- subset(metadata, predicted_celltype != "")

# Load CellChat objects
cellchatAD <- readRDS("path/to/cellchat_AD.rds")
cellchatDHF <- readRDS("path/to/cellchat_DHF.rds")

# Ensure 'images' slot exists
for (obj in list(cellchatAD, cellchatDHF)) {
  if (!"images" %in% slotNames(obj)) {
    slot(obj, "images") <- list()
  }
}

# Attach metadata
cellchatAD@meta <- metadata
cellchatDHF@meta <- metadata

# Update CellChat objects
cellchatAD <- updateCellChat(cellchatAD)
cellchatDHF <- updateCellChat(cellchatDHF)

# Store in list
cellchat.list <- list(AD = cellchatAD, DHF = cellchatDHF)

# Set CellChatDB
CellChatDB <- CellChatDB.human
options(future.globals.maxSize = 1000 * 1024^2)  # Increase memory limit

# Process each condition
plan("multisession", workers = 4)
for (name in names(cellchat.list)) {
  obj <- cellchat.list[[name]]

  # Filter valid metadata and match to expression data
  valid <- !is.na(obj@meta$predicted_celltype)
  obj@meta <- obj@meta[valid, , drop = FALSE]
  matched_cells <- intersect(rownames(obj@meta), colnames(obj@data))
  obj@meta <- obj@meta[matched_cells, , drop = FALSE]
  obj@data <- obj@data[, matched_cells]

  # Set identity and options
  obj@idents <- droplevels(factor(obj@meta$predicted_celltype))
  obj@options$group.by <- "predicted_celltype"
  obj@DB <- CellChatDB

  # Run CellChat pipeline
  obj <- subsetData(obj)
  obj <- identifyOverExpressedGenes(obj)
  obj <- identifyOverExpressedInteractions(obj)
  obj <- computeCommunProb(obj)
  obj <- filterCommunication(obj, min.cells = 10)
  obj <- computeCommunProbPathway(obj)
  obj <- aggregateNet(obj)

  # Save back to list
  cellchat.list[[name]] <- obj
}

# Optional: Save cellchat.list
# saveRDS(cellchat.list, "cellchatlist_processed.rds")

# Load processed list (if already saved)
# cellchat.list <- readRDS("path/to/cellchatlist_processed.rds")

# Merge conditions
cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

# Plot communication bubble (CD8+ T cell subsets to Bmem/CD16 Mono)
targets.use <- c("CX3CR1+ CD8", "CXCR6+ CD8", "Int CD8")

gg1 <- netVisual_bubble(
  cellchat, sources.use = 3, targets.use = targets.use,
  comparison = c(1, 2), max.dataset = 1,
  title.name = "Increased signaling in AD: B memory", angle.x = 45, remove.isolate = TRUE
)

gg2 <- netVisual_bubble(
  cellchat, sources.use = 6, targets.use = targets.use,
  comparison = c(1, 2), max.dataset = 1,
  title.name = "Increased signaling in AD: CD16 Mono", angle.x = 45, remove.isolate = TRUE
)

# Save combined plot
combined_plot <- gg2 + gg1
ggsave("figures/6sup_cellchat_plot_cd8subset_bmem_cd16mono_hiad.pdf", combined_plot, width = 7, height = 4)
