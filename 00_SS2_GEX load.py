import scanpy as sc
import pandas as pd

# Load Data 
adata_1 = sc.read_text('path/to/batch1_A11_results/batch1_results/combined/study6369-tic1718-star-fc-genecounts.txt').T
adata_2 = sc.read_text('path/to/batch1_A11_results/run_41802_lane_6_results/combined/study6369-tic1718-star-fc-genecounts.txt').T
bdata_1 = sc.read_text('path/to/home/sirawit/Smart-Seq2/samples_not_100000_results/combined/study6369-tic1718-star-fc-genecounts.txt').T
bdata_2 = sc.read_text('path/to/samples_100000_results/combined/study6369-tic1718-star-fc-genecounts.txt').T

# Concatenate
adata = adata_1.concatenate(adata_2, bdata_1, bdata_2)

#Make Gene Names Unique
adata.var_names_make_unique()

#Annotate Genes 
annot = sc.queries.biomart_annotations(
    "hsapiens",
    ["ensembl_gene_id", "external_gene_name"],
).set_index("ensembl_gene_id")

adata.var = adata.var.join(annot, how="left")

# Remove Genes Without Names 
adata = adata[:, adata.var["external_gene_name"].notna()].copy()

#Set Gene Names as Index
adata.var.set_index("external_gene_name", inplace=True)

# Save Final AnnData 
adata.write("path/to/GEX_A11+A24_gene_annotate_unprocess.h5ad")

