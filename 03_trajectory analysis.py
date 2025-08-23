import numpy as np
import scanpy as sc

adata= sc.read(‘path/dengue_ss2_celltype_annotaion.h5ad')

# Step 1: Run PAGA on Leiden clusters
sc.tl.paga(adata, groups="leiden")
sc.pl.paga(adata, color="leiden", frameon=False)

# Step 2: Initialize graph layout with PAGA positions
sc.tl.draw_graph(adata, init_pos="paga")

# Optional: custom color palette (adjust if needed)
custom_palette = {
    '0': '#d62728',   # CX3CR1+ CD8
    '1': '#ff7f0e',   # Prolif CD8
    '2': '#279e68',   # Naive/CM CD8
    '3': '#1f77b4',   # Int CD8
    '4': '#aa40fc'    # CXCR6+ CD8
}

# Step 3: Visualize draw_graph layout and PAGA connectivity
sc.pl.draw_graph(adata, color="leiden", legend_loc="on data", palette=custom_palette)
sc.pl.paga_compare(
    adata,
    threshold=0.03,
    title="",
    right_margin=0.2,
    size=10,
    edge_width_scale=0.5,
    legend_fontsize=12,
    fontsize=12,
    frameon=False,
    edges=True,
)

# Step 4: Define root cluster and compute pseudotime
adata.uns["iroot"] = np.flatnonzero(adata.obs["leiden"] == "2")[0]  # use Naive/CM CD8 as root
sc.tl.dpt(adata)

# Step 5: Plot pseudotime on graph layout
sc.pl.draw_graph(
    adata,
    color=["leiden", "dpt_pseudotime"],
    legend_loc="on data",
    palette=custom_palette,
)

