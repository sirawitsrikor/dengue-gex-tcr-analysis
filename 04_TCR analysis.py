import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scirpy as ir
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy import sparse
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# Plot settings
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams["savefig.format"] = "pdf"
sc.settings.verbosity = 1
np.set_printoptions(linewidth=180)

# ---- Load Data ---- #
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

# Merge metadata
meta_ss2 = pd.read_csv('path/ss2_tcr_metadata.csv',
                  index_col = 0)
adata_tcr.obs = adata_ss2.obs.join(meta_ss2)
adata_tcr.obs

# Chain QC & filtering
ir.tl.chain_qc(adata_tcr)
ax = ir.pl.group_abundance(adata_tcr, groupby="chain_pairing", target_col = 'hla')
adata_tcr = adata_tcr[~adata_tcr.obs["chain_pairing"].isin(["multichain", "no IR", "orphan VDJ", "orphan VJ"])]
adata_tcr = adata_tcr[adata_tcr.obs["hla"].isin(['A*11', 'A*24 & A*11']), :].copy()

# Clonotype definition and clustering
ir.pp.ir_dist(adata_tcr)
ir.tl.define_clonotypes(adata_tcr, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonotype_network(adata_tcr, min_cells=2)

var = ['donor_id','severity','clinical_phase','hla','reactivity']
for group in var:
    ax = ir.pl.clonotype_network(
    adata_tcr, base_size=10, panel_size=(5, 5),color=group
) 

# Clonotype network clustering
ir.pp.ir_dist(adata_tcr, metric="alignment", sequence="aa", cutoff=15,)
ir.tl.define_clonotype_clusters(adata_tcr, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(adata_tcr, min_cells=2, sequence="aa", metric="alignment")

for group in var:
    ir.pl.clonotype_network(
    adata_tcr, color=group, label_fontsize=9, panel_size=(7, 7), base_size=10
)
    print(group)

# Update receptor_subtype and chain pairing 
mask_tra = (
    adata_tcr.obs['IR_VJ_1_j_call'].str.startswith('TRA', na=False) &
    adata_tcr.obs['IR_VJ_1_v_call'].str.startswith('TRA', na=False)
)

adata_tcr.obs.loc[
    (adata_tcr.obs['receptor_subtype'] == 'ambiguous') & mask_tra,
    'receptor_subtype'
] = 'TRA+TRB'

adata_tcr.obs.loc[
    (adata_tcr.obs['chain_pairing'] == 'ambiguous') & mask_tra,
    'chain_pairing'
] = 'extra VJ'

# Create VJ pair of TRA and TRB chain
adata_tcr.obs['TRA_VJ'] = adata_tcr.obs['IR_VJ_1_v_call'].str.cat(
    adata_tcr.obs['IR_VJ_1_j_call'], sep='_'
)

adata_tcr.obs['TRB_VJ'] = adata_tcr.obs['IR_VDJ_1_v_call'].str.cat(
    adata_tcr.obs['IR_VDJ_1_j_call'], sep='_'
)

def plot_vj_heatmap(adata, vj_col="TRA_VJ", phase="Acute", save_path=None, dpi=300):
    """
    Plot log10 heatmap of VJ usage by severity and cell type.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with `.obs` containing:
        - 'clinical_phase'
        - 'severity'
        - 'original_celltype'
        - VJ column (e.g., 'TRA_VJ' or 'TRB_VJ')
    vj_col : str
        Column name in `.obs` for the VJ usage.
    phase : str
        Clinical phase to filter on (e.g., 'Acute', 'Convalescence').
    save_path : str, optional
        File path to save the figure (e.g., "./figures/vj_heatmap.pdf").
        If None, the figure will not be saved.
    dpi : int, default=300
        Resolution (dots per inch) for saving the figure.
    """

    # Step 1: Filter and convert types
    df = adata.obs[
        (adata.obs['clinical_phase'] == phase) &
        (~adata.obs['original_celltype'].isna())
    ].copy()
    df['severity'] = df['severity'].astype(str)
    df['original_celltype'] = df['original_celltype'].astype(str)

    # Step 2: Count VJ per Severity
    vj_counts = df.groupby(['severity', vj_col]).size().reset_index(name='count')

    # Step 3: Pivot to get VJ × Severity matrix
    vj_matrix = vj_counts.pivot(index=vj_col, columns='severity', values='count').fillna(0)

    # Step 4: Define severity order and weighted sorting rule
    severity_weights = {'DHF': 100, 'DF': 10, 'AD': 1}
    for s in severity_weights:
        if s not in vj_matrix.columns:
            vj_matrix[s] = 0
    vj_matrix['sort_score'] = sum(vj_matrix[s] * w for s, w in severity_weights.items())

    # Step 5: Get VJ in sorted order
    vj_sorted = vj_matrix.sort_values('sort_score', ascending=False).index.tolist()

    # Step 6: Grouped raw counts (no normalization)
    df_grouped = df.groupby(['severity', 'original_celltype', vj_col]).size().reset_index(name='count')
    df_grouped['group_label'] = df_grouped['severity'] + ' | ' + df_grouped['original_celltype']

    # Step 7: Pivot to heatmap matrix
    heatmap_data = df_grouped.pivot(index=vj_col, columns='group_label', values='count').fillna(0)

    # Step 8: Reorder rows
    heatmap_data = heatmap_data.reindex(index=vj_sorted).fillna(0)

    # Step 9: Reorder columns to include all severity × cell type combos
    severity_order = ['AD', 'DF', 'DHF']
    cell_type_order = ['CXCR6+ CD8', 'Prolif CD8', 'CX3CR1+ CD8', 'Int CD8', 'Naive/CM CD8']
    all_cols = [f"{s} | {ct}" for s in severity_order for ct in cell_type_order]
    heatmap_data = heatmap_data.reindex(columns=all_cols).fillna(0)

    # Step 10: Log scale transform
    heatmap_log = np.log10(heatmap_data.replace(0, 1e-4))

    # Plot
    plt.figure(figsize=(14, 10))
    sns.heatmap(heatmap_log, cmap='viridis', vmin=-4)
    plt.title(f'Log10 Raw {vj_col} usage by Cell Type ({phase})')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Severity | Cell Type')
    plt.ylabel(f'{vj_col} (ordered by DHF > DF > AD)')
    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, format="pdf", dpi=dpi, bbox_inches="tight")

    plt.show()

#Visualize VJ pair heatmap
plot_vj_heatmap(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], vj_col="TRA_VJ", phase="Acute",save_path=save_path)
plot_vj_heatmap(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], vj_col="TRB_VJ", phase="Acute",save_path=save_path)

def plot_vj_clustermap(
    adata,
    vj_col="TRA_VJ",                 # or "TRB_VJ"
    phase="Acute",                   # or "Convalescence"
    exclude_celltypes=("Naive/CM CD8",),  # set to () or None to keep all
    phase_col="clinical_phase",
    severity_col="severity",
    celltype_col="original_celltype",
    cmap="viridis",
    figsize=(14, 12),
    tiny=1e-4,                       # value for zeros before log10
    save_path=None,                  # file path for saving (e.g., "./figures/heatmap.pdf")
    dpi=300,                         # resolution when saving
    return_data=False                # return (grid, heatmap_df, log_df)
):
    """
    Unsupervised clustermap of VJ usage (TRA/TRB) by Severity | Cell Type.
    """
    # --- Step 1: Filter data ---
    obs = adata.obs
    mask = (obs[phase_col] == phase) & (~obs[celltype_col].isna())
    if exclude_celltypes:
        mask &= (~obs[celltype_col].isin(list(exclude_celltypes)))
    df = obs.loc[mask, [severity_col, celltype_col, vj_col]].copy()

    # Ensure strings
    df[severity_col] = df[severity_col].astype(str)
    df[celltype_col] = df[celltype_col].astype(str)
    df[vj_col] = df[vj_col].astype(str)

    # --- Step 2: Count per (Severity, CellType, VJ) ---
    df_grouped = (
        df.groupby([severity_col, celltype_col, vj_col])
          .size()
          .reset_index(name='count')
    )
    df_grouped['group_label'] = df_grouped[severity_col] + ' | ' + df_grouped[celltype_col]

    # --- Step 3: Pivot to VJ × (Severity | Cell Type) matrix ---
    heatmap_data = (
        df_grouped.pivot(index=vj_col, columns='group_label', values='count')
                  .fillna(0)
    )

    # --- Step 4: Remove rows/cols that are all zero (optional cleanup) ---
    heatmap_data_clean = heatmap_data.loc[
        (heatmap_data != 0).any(axis=1),
        (heatmap_data != 0).any(axis=0)
    ]

    # --- Step 5: Log10 transform (keep zeros as tiny value so background isn't white) ---
    heatmap_log = np.log10(heatmap_data_clean.replace(0, tiny))

    # --- Step 6: Clustered heatmap ---
    g = sns.clustermap(
        heatmap_log,
        cmap=cmap,
        metric='correlation',
        figsize=figsize,
        linewidths=0
    )
    plt.suptitle(
        f'Unsupervised Clustering of {vj_col} Usage by Severity | Cell Type ({phase})',
        fontsize=16, y=1.02
    )

    # Save if path provided
    if save_path:
        g.savefig(save_path, format=save_path.split('.')[-1], dpi=dpi, bbox_inches="tight")

    plt.show()

    if return_data:
        return g, heatmap_data_clean, heatmap_log

plot_vj_clustermap(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], vj_col="TRA_VJ", phase="Acute", save_path=save_path)
plot_vj_clustermap(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], vj_col="TRB_VJ", phase="Acute", save_path=save_path)

def plot_reactivity_by_group(adata, phase="Acute", reactivity_palette=None, save_path=None, dpi=300):
    """
    Plot normalized reactivity proportions by Severity | Cell Type.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with `.obs` containing:
        - 'clinical_phase'
        - 'severity'
        - 'original_celltype'
        - 'reactivity'
    phase : str, default="Acute"
        Clinical phase to filter on.
    reactivity_palette : dict, optional
        Dictionary mapping reactivity categories to colors.
        Defaults to current/previous/cross-reactive colors.
    save_path : str, optional
        File path to save the figure (e.g., "./figures/reactivity_plot.pdf").
        If None, the figure is not saved.
    dpi : int, default=300
        Resolution (dots per inch) for saving the figure.
    """

    # Step 1: Filter and copy
    data = adata.obs[
        (adata.obs['clinical_phase'] == phase) &
        (~adata.obs['original_celltype'].isna())
    ].copy()

    # Define full desired group labels in order
    severity_order = ['AD', 'DF', 'DHF']
    cell_type_order = ['CXCR6+ CD8', 'Prolif CD8', 'CX3CR1+ CD8', 'Int CD8', 'Naive/CM CD8']
    group_order = [f"{sev} | {ct}" for sev in severity_order for ct in cell_type_order]

    # Step 2: Group and count reactivity
    count_df = (
        data.groupby(['severity', 'original_celltype'])['reactivity']
        .value_counts()
        .unstack(fill_value=0)
    )

    # Step 3: Convert to string-based index
    count_df.index = [f"{sev} | {ct}" for sev, ct in count_df.index]

    # Step 4: Reindex to include all desired rows
    count_df = count_df.reindex(group_order, fill_value=0)

    # Step 5: Normalize rows to sum to 1
    normalized_df = count_df.div(count_df.sum(axis=1), axis=0)

    # Default color palette
    if reactivity_palette is None:
        reactivity_palette = {
            'current': '#1f77b4',
            'previous': '#ff7f0e',
            'cross-reactive': '#2ca02c',
        }

    # Step 6: Plot
    ax = normalized_df.plot(
        kind='bar',
        stacked=True,
        figsize=(12, 4),
        color=reactivity_palette,
        width=0.9
    )
    plt.ylabel("Proportion of reactivity")
    plt.xlabel("Severity | Cell type")
    plt.title(f"Normalized reactivity by Severity and Cell Type ({phase})")
    plt.xticks(rotation=90)
    plt.legend(title='Reactivity', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path, format=save_path.split('.')[-1], dpi=dpi, bbox_inches="tight")

    plt.show()

    return normalized_df

plot_reactivity_by_group(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], phase="Acute", save_path= save_path)

plot_vj_heatmap(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], vj_col="TRA_VJ", phase="Convalescence", save_path= save_path)
plot_vj_heatmap(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], vj_col="TRB_VJ", phase="Convalescence", save_path= save_path)
plot_reactivity_by_group(adata_tcr[adata_tcr.obs['receptor_subtype']=='TRA+TRB'], phase="Convalescence", save_path= save_path)

severity = ['AD', 'DF', 'DHF']
for var in severity:
    ir.pl.vdj_usage(
        adata_tcr[adata_tcr.obs['severity'] == var],
        full_combination=False,
        max_segments=None,
        max_ribbons=30,
    )

