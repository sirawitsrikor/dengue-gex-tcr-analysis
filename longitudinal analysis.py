import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import matplotlib.pyplot as plt
import os

from scipy import stats
from itertools import product

# --- Settings --- #
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['savefig.format'] = 'pdf'
sc.settings.verbosity = 1
os.makedirs("figures", exist_ok=True)

# --- Load metadata --- #
metadata = pd.read_csv("metadata_for_longitudinal_ss2_10x.csv", index_col=0)

# --- Load datasets --- #
# Load annotated 10X CD8+ T cells
adata_10x = sc.read("/path/to/asymdenv_annotated_rawcount.h5ad")
non_cd8 = [...]  # same list of non-CD8 types as before
adata_cd8 = adata_10x[~adata_10x.obs["CellType"].isin(non_cd8)].copy()
tcr_10x = ir.io.read_10x_vdj("DENV_TCR_merged.csv", filtered=True)
ir.pp.merge_with_ir(adata_cd8, tcr_10x)

# Load Smart-seq2 GEX + TCR
# [Insert the concise version of Smart-seq2 merging from previous scripts if needed]

# Load longitudinal TCR
long_10x = ir.io.read_10x_vdj("longitudinal_TCR_combine.csv", filtered=True)
long_10x.obs.index = long_10x.obs.index.str.replace("-1", "")

# Merge everything
adata = adata_cd8.concatenate(adata, long_10x)
adata.obs.index = adata.obs.index.str.replace("-0", "").str.replace("-1", "").str.replace("-2", "")

# Join metadata
adata.obs = adata.obs.drop(columns=[col for col in adata.obs.columns if col in metadata.columns], errors='ignore')
adata.obs = adata.obs.join(metadata)
adata = adata[adata.obs["severity"].isin(["AD", "DF", "DHF"])]

# Recode days
adata.obs["day"] = adata.obs["day"].replace({"D2": "D1", "D3": "D1"}).astype("category")

# --- Join cell type annotation --- #
df_map = pd.concat([
    pd.read_csv("mauscript_csv/10x_asymp_donor_map_gts_ka6_kw40.csv", index_col=0),
    pd.read_csv("mauscript_csv/10x_long_donor_map_gts_ka6_kw40.csv", index_col=0).rename(index=lambda x: x.replace("-0", "").replace("-1", "").replace("-2", "")),
    pd.read_csv("mauscript_csv/ss2_donor_map_gts_ka6_kw40.csv", index_col=0)
])
df_map = df_map[["predicted_id"]].rename(columns={"predicted_id": "predicted_celltype"})
df_map["predicted_celltype"] = df_map["predicted_celltype"].replace({
    0: "CX3CR1+ CD8", 1: "Int CD8", 2: "CXCR6+ CD8", 3: "Naive/CM CD8", 4: "Prolif CD8", 5: "Prolif CD8"
})
adata.obs = adata.obs.join(df_map)

# --- Clean TCR junctions and filter --- #
for col in ["IR_VDJ_1_junction", "IR_VJ_1_junction"]:
    adata.obs[col] = adata.obs[col].str.replace(r"^TGT(.*)TT[TC]$", r"\1", regex=True)

adata = adata[~adata.obs["chain_pairing"].isin(["multichain", "no IR", "orphan VDJ", "orphan VJ"])]
ir.pp.ir_dist(adata)
ir.tl.define_clonotypes(adata, receptor_arms="all", dual_ir="primary_only")

# --- Define longitudinal donors --- #
long_donors = ["I010", "I076", "I085", "I160", "T002", "T010", "YH027", "YH034"]
adata_long = adata[(adata.obs["Donor_ID"].isin(long_donors)) & (adata.obs["day"] != "F1")].copy()

# --- Set colors for CD8 subsets --- #
celltype_colors = {
    "Prolif CD8": "#e66101",
    "CXCR6+ CD8": "#5e3c99",
    "CX3CR1+ CD8": "#b2182b",
    "Int CD8": "#0571b0",
    "Naive/CM CD8": "#4daf4a"
}
adata_long.obs["predicted_celltype"] = adata_long.obs["predicted_celltype"].astype("category")
adata_long.uns["predicted_celltype_colors"] = [celltype_colors[ct] for ct in adata_long.obs["predicted_celltype"].cat.categories]

# --- Plot longitudinal clonal expansion --- #
def plot_expansion(data, severity, filename):
    ax = ir.pl.group_abundance(
        data[(data.obs.severity == severity) & (data.obs.clone_id_size > 1)],
        groupby="day",
        target_col="predicted_celltype",
        sort=["D1", "F2"],
        normalize=True,
        max_cols=10
    )
    plt.savefig(f"figures/{filename}", format="pdf", dpi=300, bbox_inches="tight")
    plt.close()

plot_expansion(adata_long, severity="DF", filename="sup7_expansion_long_cell_DF.pdf")
plot_expansion(adata_long, severity="DHF", filename="sup7_expansion_long_cell_DHF.pdf")

print("✅ Longitudinal expansion plots saved in ./figures/")
