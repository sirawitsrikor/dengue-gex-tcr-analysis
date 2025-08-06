mport numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scirpy as ir
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from scipy import sparse

# Plot settings
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams["savefig.format"] = "pdf"
sc.settings.verbosity = 1
np.set_printoptions(linewidth=180)

# ---- Load Data ---- #
# Define paths
base_path = "/mnt/icbs_shared_storage_vc_general/Dengue_asym_SMARTseq2/Dengue_asym_SMARTseq2_sirawit_2022/Smart-Seq2/"
tcr_data = {
    "adata_tcr2": ir.io.read_tracer(base_path + "samples_100000_results/tracer_assemble/in_asm/"),
    "adata_tcr1": ir.io.read_tracer(base_path + "samples_not_100000_results/tracer_assemble/in_asm/"),
    "bdata_tcr1": ir.io.read_tracer(base_path + "batch1_results/tracer_assemble"),
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

# Merge external csv-based TCR
bdata_tcr2 = pd.read_csv('fastqs_for_tracer 2_A.csv').set_index('cell_id')
bdata_tcr2.index = bdata_tcr2.index.str.replace('out-', '', regex=False)
gex_data["bdata_2"].obs = gex_data["bdata_2"].obs.join(bdata_tcr2)

# Concatenate
adata_A24 = gex_data["adata_1"].concatenate(gex_data["adata_2"])
adata_A11 = gex_data["bdata_1"].concatenate(gex_data["bdata_2"])
adata = adata_A11.concatenate(adata_A24)
clean_cell_ids(adata, ['-0-0', '-0-1', '-1-0', '-1-1'])

# Merge metadata
meta = pd.read_csv('./mauscript_csv/final_df_A11_A24_update_reactivity_tetramer_CD69_avidity_celltype.csv')
meta = meta.drop_duplicates('Cell_ID').set_index('Cell_ID')
adata.obs = adata.obs.join(meta)

# Load 10x TCR and filter CD8 cells
tcr_10x = ir.io.read_10x_vdj('DENV_TCR_merged.csv', filtered=True)
adata_full = sc.read('/mnt/icbs_shared_storage_vc_general/Dengue_asym_SMARTseq2/Dengue_asym_SMARTseq2_sirawit_2022/asymdenv_annotated_rawcount.h5ad')
exclude = [...]  # (Same list of non-CD8 types)
cd8_adata = adata_full[~adata_full.obs['CellType'].isin(exclude)].copy()
ir.pp.merge_with_ir(cd8_adata, tcr_10x)

# Final combined dataset
adata_tcr = cd8_adata.concatenate(adata)

# ---- Metadata Harmonization ---- #
def merge_obs_column(adata, target, backup):
    if backup in adata.obs:
        if target not in adata.obs:
            adata.obs[target] = adata.obs[backup]
        else:
            adata.obs[target] = adata.obs[target].fillna(adata.obs[backup])
        adata.obs.drop(columns=backup, inplace=True)

merge_obs_column(adata_tcr, 'donor_id', 'donor_ID')
merge_obs_column(adata_tcr, 'severity', 'Severity')
merge_obs_column(adata_tcr, 'age', 'Age')
merge_obs_column(adata_tcr, 'cell_type', 'CellType')
merge_obs_column(adata_tcr, 'serotype', 'RT_PCR_serotype')
adata_tcr.obs.rename(columns={'Clinical_phase': 'clinical_phase'}, inplace=True)

# Format batch & clinical phase
adata_tcr.obs['batch'] = adata_tcr.obs['batch'].astype(str).replace({'0': '10x', '1': 'ss2'}).astype('category')
adata_tcr.obs.loc[adata_tcr.obs['batch'] == '10x', 'clinical_phase'] = 'Acute'

# Chain QC & filtering
ir.tl.chain_qc(adata_tcr)
adata_tcr = adata_tcr[~adata_tcr.obs["chain_pairing"].isin(["multichain", "no IR", "orphan VDJ", "orphan VJ"])]
adata_tcr = adata_tcr[~adata_tcr.obs["donor_id"].isin([...])]

# ---- TCR Processing ---- #
# Clean junctions
for col in ['IR_VDJ_1_junction', 'IR_VJ_1_junction']:
    adata_tcr.obs[col] = adata_tcr.obs[col].str.replace(r'^TGT(.*)TT[TC]$', r'\1', regex=True)

# Clonotype definition and clustering
ir.pp.ir_dist(adata_tcr)
ir.tl.define_clonotypes(adata_tcr, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonotype_network(adata_tcr, min_cells=2)

# Define by amino acid
for col in ['IR_VDJ_1_junction_aa', 'IR_VJ_1_junction_aa']:
    adata_tcr.obs[col] = adata_tcr.obs[col].str.replace(r'^C(.*)F$', r'\1', regex=True)

ir.pp.ir_dist(adata_tcr, metric="levenshtein", sequence="aa", cutoff=2)
ir.tl.define_clonotype_clusters(adata_tcr, sequence="aa", metric="levenshtein", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(adata_tcr, sequence="aa", metric="levenshtein", min_cells=2)

# ---- Identify Shared Clonotypes ---- #
shared_ids = (
    adata_tcr.obs.groupby('cc_aa_levenshtein')['batch'].nunique()
    .reset_index().query("batch > 1")['cc_aa_levenshtein']
)
data = adata_tcr[(adata_tcr.obs['batch'] == 'ss2') | (adata_tcr.obs['cc_aa_levenshtein'].isin(shared_ids))]

# ---- Cell Type Predictions (CCA Integration) ---- #
predict = pd.read_csv('mauscript_csv/cca_integration_only_dengue_10x_ss2_celltype.csv', index_col=0)
predict.index = predict.index.to_series().apply(lambda x: f"{x}-1" if x.startswith("TA") else f"{x}-0")
data.obs = data.obs.join(predict)
data.obs['cell_predict'] = data.obs['seurat_clusters'].fillna('unidentified').astype('category')

# Recluster for analysis
data_1 = data[data.obs['cell_predict'] != 'unidentified'].copy()
ir.pp.ir_dist(data_1, metric="levenshtein", sequence="aa", cutoff=2)
ir.tl.define_clonotype_clusters(data_1, sequence="aa", metric="levenshtein", receptor_arms="all", dual_ir="any")

# ---- Visualization ---- #
severity_palette = {'AD': '#469d64', 'DF': '#f17a39', 'DHF': '#d6322d'}
data_1.obs['severity'] = data_1.obs['severity'].astype('category').cat.set_categories(['AD', 'DF', 'DHF'], ordered=True)
data_1.uns['severity_colors'] = [severity_palette[cat] for cat in data_1.obs['severity'].cat.categories]

for col in ['donor_id', 'severity', 'clinical_phase', 'batch']:
    ir.pl.clonotype_network(data_1, color=col, label_fontsize=9, panel_size=(5,5), base_size=15, show_labels=False)

# Group abundance plot
donor_order = ['IH006','IH019','YH089','T006','I147','Y016','I085','I029','I159','Y078']
ir.pl.group_abundance(data_1, groupby="donor_id", target_col="batch", normalize=False, sort=donor_order)
