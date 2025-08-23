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
adata_ss2.obs = adata_ss2.obs.join(meta_ss2)
adata_ss2.obs

# Load 10x TCR and filter CD8 cells
tcr_10x = ir.io.read_10x_vdj('path/DENV_TCR_merged.csv', filtered=True)
meta_10x = pd.read_csv('/path/asymp_10x_tcr_metadata.csv',
                  index_col = 0)
tcr_10x.obs = tcr_10x.obs.join(meta_10x)
exclude = [
    'NK', 'B naive', 'MAIT', 'CD4 TCM', 'CD4 Naive', 'NK Proliferating', 'Plasmablast',
    'CD16 Mono', 'ILC', 'B intermediate', 'NK_CD56bright', 'CD14 Mono', 'CD4 TEM', 'dnT',
    'Eryth', 'cDC2', 'Platelet', 'Treg', 'HSPC', 'CD4 CTL', 'CD4 Proliferating',
    'B memory', 'pDC', 'gdT', 'ASDC', 'cDC1'
]
cd8_adata = tcr_10x[~tcr_10x.obs['original_celltype'].isin(exclude)].copy()

# Final combined dataset
adata_tcr = cd8_adata.concatenate(adata_ss2)

# Chain QC & filtering
ir.tl.chain_qc(adata_tcr)
adata_tcr = adata_tcr[~adata_tcr.obs["chain_pairing"].isin(["multichain", "no IR", "orphan VDJ", "orphan VJ"])]
adata_tcr = adata_tcr[adata_tcr.obs["hla"].isin(['A*11', 'A*24 & A*11']), :].copy()

# ---- TCR Processing ---- #
# Clean junctions
for col in ['IR_VDJ_1_junction', 'IR_VJ_1_junction']:
    adata_tcr.obs[col] = adata_tcr.obs[col].str.replace(r'^TGT(.*)TT[TC]$', r'\1', regex=True)
# Define by amino acid
for col in ['IR_VDJ_1_junction_aa', 'IR_VJ_1_junction_aa']:
    adata_tcr.obs[col] = adata_tcr.obs[col].str.replace(r'^C(.*)F$', r'\1', regex=True)

# Clonotype definition and clustering
ir.pp.ir_dist(adata_tcr)
ir.tl.define_clonotypes(adata_tcr, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonotype_network(adata_tcr, min_cells=2)

ir.pp.ir_dist(adata_tcr, metric="levenshtein", sequence="aa", cutoff=2)
ir.tl.define_clonotype_clusters(adata_tcr, sequence="aa", metric="levenshtein", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(adata_tcr, sequence="aa", metric="levenshtein", min_cells=2)

# ---- Identify Shared Clonotypes across technological platforms ---- #
shared_ids = (
    adata_tcr.obs.groupby('cc_aa_levenshtein')['batch'].nunique()
    .reset_index().query("batch > 1")['cc_aa_levenshtein']
)
data = adata_tcr[(adata_tcr.obs['technology'] == 'ss2') | (adata_tcr.obs['cc_aa_levenshtein'].isin(shared_ids))]
data.obs['epitope']='gts'

# Recluster for analysis
ir.pp.ir_dist(data, metric="levenshtein", sequence="aa", cutoff=2)
ir.tl.define_clonotype_clusters(data, sequence="aa", metric="levenshtein", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(data, min_cells=5, sequence="aa", metric="levenshtein")

# ---- Visualization ---- #
severity_palette = {'AD': '#469d64', 'DF': '#f17a39', 'DHF': '#d6322d'}
data.obs['severity'] = data.obs['severity'].astype('category').cat.set_categories(['AD', 'DF', 'DHF'], ordered=True)
data.uns['severity_colors'] = [severity_palette[cat] for cat in data.obs['severity'].cat.categories]

var = ['donor_id', 'severity', 'clinical_phase', 'technology']
for col in var:
    ir.pl.clonotype_network(data, color=col, label_fontsize=9, panel_size=(5,5), base_size=15, show_labels=False)

# Group abundance plot
donor_order = ['IH006','IH019','YH089','T006','I147','Y016','I085','I029','I159','Y078']
ir.pl.group_abundance(data, groupby="donor_id", target_col="technology", normalize=False, sort=donor_order)



df_epitope = data.obs[['epitope']]
adata_tcr.obs = adata_tcr.obs.join(df_epitope)

adata_tcr.obs['epitope'] = (
    adata_tcr.obs['epitope']
    .astype("category")
    .cat.add_categories("other")
    .fillna("other")
)

# Generate group abundance plot
ir.pl.group_abundance(
    adata_tcr[adata_tcr.obs['clinical_phase'] == 'Acute'],
    target_col="epitope",
    groupby="severity",
    normalize=True, sort=['AD','DF', 'DHF']
)
plt.show()

##Visualize integrated celltype proportion 
# Step 1: Copy data and define expansion
df = data.obs[(~data.obs['integrate_celltype'].isna())&(df['clinical_phase']=='Acute')]
df['expansion'] = df['clone_id_size'].apply(lambda x: 'non-expanded' if x == 1 else 'expanded')

# Step 2: Group and count
grouped = df.groupby(['severity', 'expansion', 'integrate_celltype']).size().reset_index(name='count')

# Step 3: Normalize each (severity, expansion) group to sum to 1
grouped['fraction'] = grouped.groupby(['severity', 'expansion'])['count'].transform(lambda x: x / x.sum())

# Step 4: Pivot to get structure suitable for plotting later
pivot = grouped.pivot_table(index=['severity', 'expansion'], columns='integrate_celltype', values='fraction', fill_value=0)

# Define consistent colors for 5 seurat_clusters (cell types)
cell_types = pivot.columns.tolist()
color_map = {
    
   'CX3CR1+ CD8': '#d62728',  # red
    'Prolif CD8': '#ff7f0e',  # orange
    'Int CD8': '#1f77b4',  # blue
    'Naive/CM CD8': '#2ca02c',  # green 
    'CXCR6+ CD8': '#9467bd'   # purple
}

severity_order = ['AD', 'DF', 'DHF']
expansion_order = ['non-expanded', 'expanded']

# Define gap between groups
gap = 0.05
bar_width = 0.35
x = np.arange(len(severity_order))  # Base positions for severity groups

# Plot setup
plt.figure(figsize=(8, 4))
for i, exp_type in enumerate(expansion_order):
    offset = (-bar_width / 2 - gap / 2) if exp_type == 'non-expanded' else (bar_width / 2 + gap / 2)
    alpha = 1 if exp_type == 'expanded' else 0.7

    subset = pivot.loc[(slice(None), exp_type), :].reindex(severity_order, level=0)
    bottom = np.zeros(len(severity_order))

    for cell in cell_types:
        values = subset[cell].values
        
        # Convert hex to RGBA only for facecolor
        face_color = mcolors.to_rgba(color_map[cell], alpha=alpha)
        
        plt.bar(
            x + offset,
            values,
            bottom=bottom,
            width=bar_width,
            color=face_color,
            edgecolor='black',     # stays fully opaque
            linewidth=0.5,
            label=cell if i == 0 else None
        )
        bottom += values

# Custom legend handles
cell_type_handles = [Patch(facecolor=color_map[cell], edgecolor='black', label=cell) for cell in cell_types]
expansion_handles = [
    Patch(facecolor='white', edgecolor='black', hatch='.', label='non-expanded'),
    Patch(facecolor='white', edgecolor='black', label='expanded')
]

# Legends
legend1 = plt.legend(handles=cell_type_handles, title='Cell Type (predicted_celltype)', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.gca().add_artist(legend1)  # Keep first legend
plt.legend(handles=expansion_handles, title='Expansion', bbox_to_anchor=(1.05, 0.5), loc='upper left')

plt.xticks(x, severity_order)
plt.ylabel("Proportion")
plt.title("Normalized Stacked Cell Subsets per Severity and Expansion")
plt.tight_layout()
plt.show()
