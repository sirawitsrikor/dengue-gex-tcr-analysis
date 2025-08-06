import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# Optional: customize scanpy/matplotlib styles
sc.set_figure_params(dpi=100, facecolor="white")

# Load raw .h5ad file
adata = sc.read('/content/drive/MyDrive/smart_seq2/A11_A24_GEX_TCR_gene_annotate_unprocessed.h5ad')

# Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.total_counts < 4e6, :]
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 17.5, :]

# Save raw counts layer
adata.layers["counts"] = adata.X.copy()

# Normalization and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

# Highly variable gene selection within batch
sc.pp.highly_variable_genes(adata, batch_key="batches")

# Exclude unwanted genes from HVGs
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
sc.pl.umap(adata, color='leiden')
