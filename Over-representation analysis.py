import pandas as pd
import time
Import scanpy as sc
import gseapy as gp

adata= sc.read(‘path/dengue_ss2_celltype_annotaion.h5ad')

adata.uns['log1p']["base"] = None

#find marker genes
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=30, sharey=False)

# Get DEG results
result = adata_a.uns['rank_genes_groups']
groups = result['names'].dtype.names

# Convert to DataFrame
degs = pd.DataFrame({
    f"{g}_{k}": result[k][g]
    for g in groups
    for k in ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
})

# Define group names and labels
group_keys = {
    '0': 'CX3CR1+ CD8',
    '1': 'Prolif CD8',
    '2': 'Naive/CMCD8',
    '3': 'int CD8',
    '4': 'CXCR6+ CD8'
}

gene_set = 'MSigDB_Hallmark_2020'
all_results = []

for gid, label in group_keys.items():
    padj = degs[f'{gid}_pvals_adj']
    logfc = degs[f'{gid}_logfoldchanges']
    genes = degs[f'{gid}_names']

    up_genes = genes[(padj < 0.05) & (logfc > 0)].dropna().tolist()
    down_genes = genes[(padj < 0.05) & (logfc < 0)].dropna().tolist()

    if not up_genes and not down_genes:
        print(f"Skipping {label}: no significant genes")
        continue

    for gene_list, direction in zip([up_genes, down_genes], ['UP', 'DOWN']):
        if not gene_list:
            continue
        try:
            enr = gp.enrichr(gene_list=gene_list, gene_sets=gene_set, outdir=None)
            df = enr.res2d.copy()
            df['Term'] = df['Term'].str.split(" \(GO").str[0]
            df['group'] = label
            df['UP_DW'] = direction
            all_results.append(df)
            time.sleep(1.5)
        except Exception as e:
            print(f"Enrichr failed for {label} {direction}: {e}")

# Combine results
df = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()

