import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import re

# Load raw .h5ad file
adata = sc.read('path/A11_A24_GEX_TCR_gene_annotate_unprocessed.h5ad')

#add metadata
meta = pd.read_csv('path/dengue_ss2_metadata.csv', index_col=0)
adata.obs = adata.obs.join(meta)

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.total_counts < 4e6, :]
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 17.5, :]

sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)

# Save raw counts layer
adata.layers["counts"] = adata.X.copy()

# Normalization and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

# Highly variable gene selection within batch
sc.pp.highly_variable_genes(adata, batch_key="plate")

# Exclude TCR/BCR related genes from HVGs
for gene in adata.var.index:
    if re.search(r'^MT\.', gene):  # Mitochondrial genes
        adata.var.at[gene, 'highly_variable'] = False
    if re.search(r'^TR[DGAB][VDJC]', gene):  # TCR/BCR genes
        adata.var.at[gene, 'highly_variable'] = False
    if gene == 'TRAV1-2':  # Include MAIT-specific gene
        adata.var.at[gene, 'highly_variable'] = True

# Visualize HVGs
sc.pl.highly_variable_genes(adata)

# Subset to HVGs
adata = adata[:, adata.var.highly_variable]

#Score cell cycle genes
cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
 sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# Regress out unwanted sources of variation
sc.pp.regress_out(adata, ['n_genes_by_counts', 'S_score', 'G2M_score'])

# Scale the data
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Scanorama integration
sc.external.pp.scanorama_integrate(
    adata, 
    key='batch', 
    basis='X_pca', 
    adjusted_basis='X_scanorama',
    knn=20, 
    sigma=3, 
    approx=True, 
    alpha=0.1, 
    batch_size=5000
)

sc.pp.neighbors(adata, n_neighbors=13, n_pcs=6, use_rep = 'X_scanorama')
sc.tl.umap(adata)
sc.pl.umap(adata)

# Clustering
sc.tl.leiden(adata, resolution=0.32)

# Rename Leiden clusters to celltype labels
cluster_to_celltype = {
    '0': 'CX3CR1+ CD8',
    '1': 'Prof CD8',
    '2': 'Naive:CM CD8',
    '3': 'IntCD8',
    '4': 'CXCR6+ CD8'
}

# Create a new column in adata.obs for cell type annotations
adata.obs['celltype'] = adata.obs['leiden'].map(cluster_to_celltype)

# UMAP visualization
sc.pl.umap(adata, color='celltype')

# Define CD8_Cytotoxic gene module 
CD8_Cytotoxic = [
    'CCL5', 'GZMK', 'GNLY', 'TRGC2', 'FGFBP2', 'C1orf21', 'KLRF1', 'FCGR3A', 'PTGDR', 'KLRC2',
    'EOMES', 'S1PR5', 'CLIC3', 'AOAH', 'CADM1', 'TRGC1', 'DTHD1', 'LILRB1', 'SAMD3', 'ZNF683',
    'KLRD1', 'NCR1', 'FAM49A', 'KLRG1', 'CTSW', 'CD244', 'CMC1', 'APOBEC3H', 'CST7', 'CX3CR1',
    'FCRL6', 'TMCC3', 'PLA2G16', 'TYROBP', 'TPRG1', 'C12orf75', 'PLCG2', 'PLEK', 'RCAN2', 'DKK3',
    'ADRB2', 'FCRL3', 'NKG7', 'PPP2R2B', 'SYNGR1', 'KLRC4', 'HLA-DPB1', 'DAPK2', 'F2R', 'KIR3DL2',
    'B3GAT1', 'CD8B', 'TTC16', 'GALNT3', 'SCD5', 'PDGFD', 'ABCB1', 'MXRA7', 'CTBP2', 'CD8A', 'ZEB2',
    'SYTL2', 'CHN2', 'FGR', 'TGFBR3', 'SETBP1', 'COLGALT2', 'KIR2DL4', 'FKBP1B', 'ADGRG1'
]

# Score the gene module
sc.tl.score_genes(adata, gene_list=CD8_Cytotoxic, score_name='CD8_Cytotoxic_score')

# Visualize the score across clusters 
with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata, 'CD8_Cytotoxic_score', groupby='celltype', stripplot=False, inner='box')

#savesfile
adata.write(‘path/dengue_ss2_celltype_annotaion.h5ad’)
