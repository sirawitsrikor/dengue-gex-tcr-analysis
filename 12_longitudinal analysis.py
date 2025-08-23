import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import matplotlib.pyplot as plt
import os

from scipy import stats
from itertools import product

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['savefig.format'] = 'pdf'
sc.settings.verbosity = 1

# ---- Load Data SS2, 10x & longitudinal TCR ---- #
# Define paths
base_path = “path”
tcr_data = {
    "adata_tcr2": ir.io.read_tracer(base_path + "samples_100000_results/tracer_assemble/in_asm/"),
    "adata_tcr1": ir.io.read_tracer(base_path + "samples_not_100000_results/tracer_assemble/in_asm/"),
    "bdata_tcr1": ir.io.read_tracer(base_path + "batch1_results/tracer_assemble"),
    "bdata_tcr2": ir.io.read_tracer(base_path + "run_41802_lane_6_results/tracer_assemble")
}

gex_data = {
    "adata_1": sc.read_text(base_path + "samples_not_100000_results/combined/study6369-tic1718-star-fc-genecounts.txt").T,
    "adata_2": sc.read_text(base_path + "samples_100000_results/combined/study6369-tic1718-star-fc-genecounts.txt").T,
    "bdata_1": sc.read_text(base_path + "batch1_results/combined/study6369-tic1718-star-fc-genecounts.txt").T,
    "bdata_2": sc.read_text(base_path + "run_41802_lane_6_results/combined/study6369-tic1718-star-fc-genecounts.txt").T,
}

# Clean cell names
def clean_cell_ids(adata, patterns):
    for pattern in patterns:
        adata.obs.index = adata.obs.index.str.replace(pattern, '', regex=False)

for key in tcr_data:
    clean_cell_ids(tcr_data[key], ['out-'])

# Merge TCR with GEX
ir.pp.merge_with_ir(gex_data["adata_1"], tcr_data["adata_tcr1"])
ir.pp.merge_with_ir(gex_data["adata_2"], tcr_data["adata_tcr2"])
ir.pp.merge_with_ir(gex_data["bdata_1"], tcr_data["bdata_tcr1"])
ir.pp.merge_with_ir(gex_data["bdata_2"], tcr_data["bdata_tcr2"])

# Concatenate
adata_A24 = gex_data["adata_1"].concatenate(gex_data["adata_2"])
adata_A11 = gex_data["bdata_1"].concatenate(gex_data["bdata_2"])
adata_ss2 = adata_A11.concatenate(adata_A24)
clean_cell_ids(adata_ss2, ['-0-0', '-0-1', '-1-0', '-1-1'])

adata_10x = ir.io.read_10x_vdj(
    "path/to/DENV_TCR_merged.csv",
    filtered=True,
)

long_10x = ir.io.read_10x_vdj(
    "path/to/longitudinal_TCR_combine.csv",
    filtered=True,
)
long_10x.obs.index = long_10x.obs.index.str.replace("-1", "", regex=False)

# 7) Concatenate data
adata = adata_10x.concatenate(adata_ss2, long_10x)

# Strip common suffixes that sometimes linger after merges
adata.obs.index = (
    adata.obs.index.str.replace("-0", "", regex=False)
                   .str.replace("-1", "", regex=False)
                   .str.replace("-2", "", regex=False)
)

# Add metadata and harmonize final indices
metadata = pd.read_csv("longitudinal_metadata.csv", index_col=0)
adata.obs = adata.obs.join(metadata)


# Filter to CD8 T cells & basic TCR QC / cleanup
adata = adata[~adata.obs.get("predicted_celltype").isna()].copy()

# TCR quality cleanup
ir.tl.chain_qc(adata)

# Remove flanking 'TGT' at start and 'TTT/TTC' at end in junctions to harmonize across technological platforms
for col in ["IR_VDJ_1_junction", "IR_VJ_1_junction"]:
    if col in adata.obs.columns:
        adata.obs[col] = adata.obs[col].str.replace(r"^TGT(.*)TT[TC]$", r"\1", regex=True)

# Drop problematic chain-pairings
if "chain_pairing" in adata.obs.columns:
    drop_states = {"multichain", "no IR", "orphan VDJ", "orphan VJ"}
    adata = adata[~adata.obs["chain_pairing"].isin(drop_states)].copy()

# Compute IR distances & define clonotypes
ir.pp.ir_dist(adata)
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")

# Longitudinal subset & consistent colors
long_donors = ["I010", "I076", "I085", "I160", "T002", "T010", "YH027", "YH034"]
adata_long = adata[adata.obs["donor_id"].isin(long_donors)].copy()

celltype_colors = {
    "Prolif CD8":   "#e66101",
    "CXCR6+ CD8":   "#5e3c99",
    "CX3CR1+ CD8":  "#b2182b",
    "Int CD8":      "#0571b0",
    "Naive/CM CD8": "#4daf4a",
}

if "predicted_celltype" in adata_long.obs.columns:
    adata_long.obs["predicted_celltype"] = adata_long.obs["predicted_celltype"].astype("category")
    cats = list(adata_long.obs["predicted_celltype"].cat.categories)
    adata_long.uns["predicted_celltype_colors"] = [celltype_colors.get(ct, "#999999") for ct in cats]


# Create plot dunction
def plot_expansion(data, severity):
    """
    Stacks clonal (clone_id_size>1) abundance by day and CD8 subset,
    normalized within-day. Expects columns: 'severity', 'clone_id_size',
    'day', and 'predicted_celltype'.
    """
    subset = data[(data.obs["severity"] == severity) & (data.obs["clone_id_size"] > 1)]
    # If you have a known day order, define here (example):
    day_order = ["D1","F2"]
    ax = ir.pl.group_abundance(
        subset,
        groupby="day",
        target_col="predicted_celltype",
        sort=day_order,  
        normalize=True,
        max_cols=10,
    )
    return ax

# Example usage:
plot_expansion(adata_long, severity="DF")
plot_expansion(adata_long, severity="DHF")


# Donor I085: shared clones across full timepoints & per-clone composition
adata_i085 = adata[adata.obs['donor_id'] == 'I085'].copy()

# TCR processing
ir.pp.ir_dist(adata_i085)
ir.tl.define_clonotypes(adata_i085, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonotype_network(adata_i085, min_cells=2)

# Keep clones present on ≥2 different days
shared_ids = (
    adata_i085.obs.groupby('clone_id')['day'].nunique()
    .loc[lambda s: s > 2].index
)
adata_shared = adata_i085[adata_i085.obs['clone_id'].isin(shared_ids)].copy()

# Counts + % per clone × day × cell type
celltypes = list(celltype_colors) if 'celltype_colors' in locals() else \
    ["Prolif CD8","CXCR6+ CD8","CX3CR1+ CD8","Int CD8","Naive/CM CD8"]

obs = adata_shared.obs.copy()
obs["predicted_celltype"] = pd.Categorical(obs["predicted_celltype"], categories=celltypes)

df_i085 = (
    obs.groupby(["clone_id","day","predicted_celltype"], observed=False).size()
       .unstack("predicted_celltype").reindex(columns=celltypes, fill_value=0)
       .stack().rename("cell_count").reset_index()
       .rename(columns={"clone_id":"clone","predicted_celltype":"cell_type"})
)
g = df_i085.groupby(["clone","day"])["cell_count"]
df_i085["total_cells_per_clone_day"] = g.transform("sum")
df_i085["percentage"] = 100 * df_i085["cell_count"] / df_i085["total_cells_per_clone_day"]
df_i085 = df_i085.sort_values(["clone","day","cell_type"], kind="stable").reset_index(drop=True)

