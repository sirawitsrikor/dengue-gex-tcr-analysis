# Load required libraries
library(CellChat)
library(future)
library(ggplot2)
library(patchwork)

# Load metadata 
metadata <- read.csv(
  "path/to/cellchat_metadata.csv",
  row.names = "index"
)
rownames(metadata) <- gsub("-1$", "", rownames(metadata))
metadata <- metadata[metadata$predicted_celltype != "", ]

# Load CellChat objects 
cellchatAD  <- readRDS("path/to/cellchat_AD.rds")
cellchatDHF <- readRDS("path/to/cellchat_DHF.rds")

# Ensure 'images' slot exists
if (!"images" %in% slotNames(cellchatAD))  slot(cellchatAD,  "images") <- list()
if (!"images" %in% slotNames(cellchatDHF)) slot(cellchatDHF, "images") <- list()

# Attach metadata
cellchatAD@meta  <- metadata
cellchatDHF@meta <- metadata

# Update objects
cellchatAD  <- updateCellChat(cellchatAD)
cellchatDHF <- updateCellChat(cellchatDHF)

# List + DB + options
cellchat.list <- list(AD = cellchatAD, DHF = cellchatDHF)
CellChatDB    <- CellChatDB.human
options(future.globals.maxSize = 1000 * 1024^2)
future::plan("multisession", workers = 4)

# Rebuild each object with new annotation
for (i in seq_along(cellchat.list)) {
  # Filter rows with non-NA predicted_celltype
  valid_cells <- !is.na(cellchat.list[[i]]@meta$predicted_celltype)
  cellchat.list[[i]]@meta <- cellchat.list[[i]]@meta[valid_cells, , drop = FALSE]

  # Match metadata rows to data columns
  cell_names <- intersect(rownames(cellchat.list[[i]]@meta), colnames(cellchat.list[[i]]@data))
  cellchat.list[[i]]@meta <- cellchat.list[[i]]@meta[cell_names, , drop = FALSE]
  cellchat.list[[i]]@data <- cellchat.list[[i]]@data[, cell_names, drop = FALSE]

  # Set idents and grouping
  celltype_clean <- droplevels(factor(cellchat.list[[i]]@meta$predicted_celltype))
  cellchat.list[[i]]@idents <- celltype_clean
  cellchat.list[[i]]@options$group.by <- "predicted_celltype"

  # Update + pipeline
  cellchat.list[[i]] <- updateCellChat(cellchat.list[[i]])
  cellchat.list[[i]]@DB <- CellChatDB
  cellchat.list[[i]] <- subsetData(cellchat.list[[i]])
  cellchat.list[[i]] <- identifyOverExpressedGenes(cellchat.list[[i]])
  cellchat.list[[i]] <- identifyOverExpressedInteractions(cellchat.list[[i]])
  cellchat.list[[i]] <- computeCommunProb(cellchat.list[[i]])
  cellchat.list[[i]] <- filterCommunication(cellchat.list[[i]], min.cells = 10)
  cellchat.list[[i]] <- computeCommunProbPathway(cellchat.list[[i]])
  cellchat.list[[i]] <- aggregateNet(cellchat.list[[i]])
}

# ---- merge ----
cellchat <- mergeCellChat(cellchat.list, add.names = names(cellchat.list))

# ---- Plots ----
targets.use <- c("CX3CR1+ CD8", "CXCR6+ CD8", "Int CD8")
gg1 <- netVisual_bubble(
  cellchat, sources.use = 3, targets.use = targets.use,
  comparison = c(1, 2), max.dataset = 1,
  title.name = "Increased signaling in AD: B memory",
  angle.x = 45, remove.isolate = TRUE
)
gg2 <- netVisual_bubble(
  cellchat, sources.use = 6, targets.use = targets.use,
  comparison = c(1, 2), max.dataset = 1,
  title.name = "Increased signaling in AD: CD16 Mono",
  angle.x = 45, remove.isolate = TRUE
)

combined_plot <- gg2 + gg1
combined_plot


# Save combined plot
combined_plot <- gg2 + gg1
ggsave("figures/6sup_cellchat_plot_cd8subset_bmem_cd16mono_hiad.pdf", combined_plot, width = 7, height = 4)
