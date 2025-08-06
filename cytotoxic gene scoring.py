import scanpy as sc
from matplotlib.pyplot import rc_context

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

# Visualize the score across clusters (use 'celltype' if previously defined)
with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata, 'CD8_Cytotoxic_score', groupby='leiden', stripplot=False, inner='box')
