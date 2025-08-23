import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
import os

plt.rcParams["savefig.format"] = "pdf"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---- Load data ---- #
df = pd.read_csv('total_raw_data_CBA_dengue.csv')
df_cell.to_csv('10x_predicted_celltype_expanded.csv', index_col = 0)

# clean cytokine data
df = df.apply(lambda col: col.apply(lambda x: float(str(x).replace('<','').replace('↓','')) / 2
                                    if any(c in str(x) for c in ['<', '↓']) else x))
df = df.apply(lambda col: col.astype(str).str.replace('<','').str.replace('↓',''))

# Calculate CD8 cell type frequency per donor 
counts = df_cell.groupby(["donor_ID", "cell_predicted"]).size().unstack(fill_value=0)
total_counts = counts.sum(axis=1)
percentages = counts.div(total_counts, axis=0) * 100
percentages.reset_index(inplace=True)

# Merge cytokine and CD8 frequency data 
merged = percentages.merge(df, on="donor_ID")
for col in merged.columns:
    converted = pd.to_numeric(merged[col], errors="coerce")
    if converted.notna().sum() / len(merged) > 0.9:
        merged[col] = converted

# Compute Spearman correlation and p-values 
def spearman_corr_pval(df):
    cols = df.columns
    rho = pd.DataFrame(index=cols, columns=cols, dtype=float)
    pval = pd.DataFrame(index=cols, columns=cols, dtype=float)
    for i in cols:
        for j in cols:
            r, p = stats.spearmanr(df[i], df[j])
            rho.loc[i, j] = r
            pval.loc[i, j] = p
    return rho, pval

df_numeric = merged.select_dtypes(include=["float64", "int64"])
corr, pval = spearman_corr_pval(df_numeric)
corr_filtered = corr.where(pval < 0.05, 0)

# CD8-specific correlation subset 
cd8_cols = [
    "CX3CR1+ CD8", "CXCR6+ CD8", "Int CD8", 
    "Naive/CM CD8", "Prolif CD8-1", "Prolif CD8-2"
]
non_cd8 = ~corr_filtered.index.isin(cd8_cols)
cytokine_assoc = corr_filtered.loc[non_cd8, cd8_cols]
associated_cytokines = cytokine_assoc[(cytokine_assoc != 0).any(axis=1)].index.tolist()
subset_vars = cd8_cols + associated_cytokines
subset_corr = corr_filtered.loc[subset_vars, subset_vars]

# Hierarchical clustering 
linkage = hierarchy.linkage(subset_corr, method="ward")
ordered_idx = hierarchy.dendrogram(linkage, no_plot=True)["leaves"]
corr_ordered = subset_corr.iloc[ordered_idx, ordered_idx]

# Plot heatmap 
plt.figure(figsize=(10, 8))
sns.heatmap(
    corr_ordered,
    cmap="coolwarm",
    mask=(corr_ordered == 0),
    linewidths=0.5,
    linecolor="white",
    xticklabels=True,
    yticklabels=True,
    cbar=True,
    square=True
)
plt.xticks(rotation=90, fontsize=10)
plt.yticks(rotation=0, fontsize=10)
plt.title("CD8-Associated Cytokine Correlation (p < 0.05)", fontsize=14)
plt.tight_layout()

# ---- set function ---- 
def plot_var(var):
    r, p = stats.spearmanr(merged[target], merged[var])

    fig, ax = plt.subplots(figsize=(5, 3))
    sns.scatterplot(data=merged, x=target, y=var, hue="severity",
                    palette=palette, edgecolor="black", s=50, ax=ax)
    sns.regplot(data=merged, x=target, y=var, scatter=False, color="black", ax=ax)

    ax.set(xlabel=target, ylabel=var, title=f"{var} vs {target}")
    ax.legend(title="severity", bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    ax.text(0.98, 0.02, f"Spearman r={r:.2f}\np={p:.4f}",
            transform=ax.transAxes, ha="right", va="bottom")

    fig.tight_layout()
    plt.show()


# ---- Plot cytokine correlation---- #
for v in vars_to_plot:
    plot_var(v)

