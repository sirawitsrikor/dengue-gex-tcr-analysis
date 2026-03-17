import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# settings
group_cols = ["donor_id", "severity"]
obs = adata_tcr.obs.copy()

# remove allele suffix, e.g. TRBV7-9*01 -> TRBV7-9
def strip_allele(x):
    return np.nan if pd.isna(x) else str(x).split("*")[0]

# required columns
col_map = {
    "TRA_V": "IR_VJ_1_v_call",
    "TRA_J": "IR_VJ_1_j_call",
    "TRB_V": "IR_VDJ_1_v_call",
    "TRB_J": "IR_VDJ_1_j_call",
    "TRB_D": "IR_VDJ_1_d_call",   # optional
}

# copy and clean columns
for new, old in col_map.items():
    if old in obs.columns:
        obs[new] = obs[old].map(strip_allele)

# keep rows with main TCR info
obs2 = obs.dropna(subset=["TRA_V", "TRA_J", "TRB_V", "TRB_J"] + group_cols).copy()

# helper: group-level frequency table
def usage_table(df, feature, prefix):
    tab = (
        df.groupby(group_cols + [feature])
          .size()
          .unstack(fill_value=0)
    )
    return tab.div(tab.sum(axis=1), axis=0).add_prefix(prefix)

# build feature matrix
parts = [
    usage_table(obs2, "TRA_V", "TRA_V_"),
    usage_table(obs2, "TRA_J", "TRA_J_"),
    usage_table(obs2, "TRB_V", "TRB_V_"),
    usage_table(obs2, "TRB_J", "TRB_J_"),
]

if "TRB_D" in obs2.columns and obs2["TRB_D"].notna().any():
    parts.append(usage_table(obs2, "TRB_D", "TRB_D_"))

X = pd.concat(parts, axis=1).fillna(0)

print("Feature matrix shape:", X.shape)


pca = PCA(n_components=2, random_state=0)
pca_df = pd.DataFrame(
    pca.fit_transform(StandardScaler().fit_transform(X)),
    index=X.index,
    columns=["PC1", "PC2"]
).reset_index()

plt.figure(figsize=(6, 5))
for label, sub in pca_df.groupby(group_cols[-1]):
    plt.scatter(sub["PC1"], sub["PC2"], label=label, alpha=0.8)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.legend(title=group_cols[-1], bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()


from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

Z = linkage(pca_df[["PC1", "PC2"]], method="ward")
labels = [f"{d} | {g}" for d, g in zip(pca_df["donor_id"], pca_df["severity"])]

plt.figure(figsize=(7, 4))
dendrogram(Z, labels=labels, leaf_rotation=90)
plt.title("Hierarchical clustering (PC1 + PC2)")
plt.ylabel("Distance")
plt.tight_layout()
plt.show()
