#define a function

def plot_vj_heatmap_unsupervised(
    adata,
    vj_col="TRB_VJ",
    phase="Acute",
    severity_order=("AD", "DF", "DHF"),
    cluster_rows=True,
    cluster_cols=False,
    metric="correlation",
    method="average",
    row_scale=None,
    pseudocount=1e-4,
    vmin=-4,
    cmap="RdBu_r",
    save_path=None,
    dpi=300,

    # --- layout controls ---
    show_row_labels=False,
    show_row_dendrogram=True,
    fig_width=3.2,
    row_height=0.12,

    # --- NEW border controls ---
    cell_border=False,
    border_width=0.15,
    border_color="0.7",
):

    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    df = adata.obs.loc[
        (adata.obs["clinical_phase"] == phase),
        ["severity", vj_col]
    ].copy()

    df["severity"] = df["severity"].astype(str)
    df = df.dropna(subset=[vj_col])
    df = df[df[vj_col].astype(str).str.len() > 0]

    vj_counts = df.groupby(["severity", vj_col]).size().reset_index(name="count")

    # ---- RAW COUNTS ----
    mat = (
        vj_counts.pivot(index=vj_col, columns="severity", values="count")
        .fillna(0)
        .astype(float)
    )

    for s in severity_order:
        if s not in mat.columns:
            mat[s] = 0.0

    mat = mat.loc[:, list(severity_order)]
    mat = mat.loc[mat.sum(axis=1) > 0]

    # ---- TRANSFORM ----
    mat_plot = np.log10(mat.replace(0, pseudocount))

    if row_scale == "zscore":
        mu = mat_plot.mean(axis=1)
        sd = mat_plot.std(axis=1).replace(0, np.nan)
        mat_plot = mat_plot.sub(mu, axis=0).div(sd, axis=0).fillna(0)

    if metric == "correlation":
        mat_plot = mat_plot.loc[mat_plot.var(axis=1) > 0]

    mat_plot = mat_plot.replace([np.inf, -np.inf], np.nan).dropna(axis=0)

    if mat_plot.shape[0] < 2:
        raise ValueError(
            f"Too few VJ rows to cluster after filtering (n={mat_plot.shape[0]})."
        )

    # ---- border logic ----
    lw = border_width if cell_border else 0

    # ---- size ----
    fig_height = max(2.2, row_height * mat_plot.shape[0])

    g = sns.clustermap(
        mat_plot,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        metric=metric,
        method=method,
        cmap=cmap,
        vmin=vmin if row_scale != "zscore" else None,
        center=0 if row_scale == "zscore" else None,
        linewidths=lw,
        linecolor=border_color,
        figsize=(fig_width, fig_height),
        yticklabels=show_row_labels,
        dendrogram_ratio=(0.12, 0.02),
        cbar_pos=(0.02, 0.85, 0.04, 0.12),
    )

    # ---- labels ----
    g.fig.suptitle(f"Unsupervised {vj_col} usage ({phase})", y=1.01, fontsize=10)
    g.ax_heatmap.set_xlabel("Severity")
    g.ax_heatmap.set_ylabel("" if not show_row_labels else vj_col)

    if not show_row_labels:
        g.ax_heatmap.set_yticks([])
        g.ax_heatmap.tick_params(left=False)

    if not show_row_dendrogram:
        g.ax_row_dendrogram.set_visible(False)

    g.ax_heatmap.tick_params(axis="x", labelrotation=0)

    plt.tight_layout()

    if save_path:
        g.fig.savefig(save_path, dpi=dpi, bbox_inches="tight")

    plt.show()
    return g, mat

#plot

for phase in ["Acute", "Convalescence"]:
    for vj in ["TRA_VJ", "TRB_VJ"]:

        g, raw_mat = plot_vj_heatmap_unsupervised(
            adata_tcr,
            vj_col=vj,
            phase=phase,
            severity_order=("DF","DHF"),
            metric="euclidean",
            cluster_rows=True,
            cluster_cols=False,
            show_row_labels=False,
            fig_width=12,
            row_height=0.1,
            cell_border=False
        )

        g

#Monte Carlo Fisher-Freeman-Halton test

import numpy as np
from scipy.stats import chi2_contingency

def fisher_rxC_mc(table_df, B=20000, seed=1):
    """Monte Carlo Fisher-Freeman-Halton test for RxC table"""

    rng = np.random.default_rng(seed)
    table = table_df.values

    # observed statistic
    chi2_obs = chi2_contingency(table, correction=False)[0]

    # generate row/column labels
    row_labels = np.repeat(np.arange(table.shape[0]), table.sum(axis=1))
    col_labels = np.repeat(np.arange(table.shape[1]), table.sum(axis=0))

    count = 0
    for _ in range(B):
        rng.shuffle(col_labels)

        sim = np.zeros_like(table)
        for r, c in zip(row_labels, col_labels):
            sim[r, c] += 1

        if chi2_contingency(sim, correction=False)[0] >= chi2_obs:
            count += 1

    p = (count + 1) / (B + 1)

    return {"chi2_observed": chi2_obs, "p_value": p, "B": B}

import pandas as pd

files = {
    "Acute_TRA_VJ": "./heatmap_TRA_VJ_raw_counts_Acute.csv",
    "Acute_TRB_VJ": "./heatmap_TRB_VJ_raw_counts_Acute.csv",
    "Convalescence_TRA_VJ": "./heatmap_TRA_VJ_raw_counts_Con.csv",
    "Convalescence_TRB_VJ": "./heatmap_TRB_VJ_raw_counts_Con.csv"
}

results = []

for name, path in files.items():
    df = pd.read_csv(path)

    # keep only severity columns that exist
    cols = [c for c in ["AD", "DF", "DHF"] if c in df.columns]

    # skip if fewer than 2 severity groups
    if len(cols) < 2:
        print(f"Skipping {name}: not enough severity columns")
        continue

    # first column assumed to be cluster / VJ name
    table_df = df.set_index(df.columns[0])[cols]

    # optional: remove rows with all zeros
    table_df = table_df.loc[table_df.sum(axis=1) > 0]

    # skip if table is too small
    if table_df.shape[0] < 2:
        print(f"Skipping {name}: not enough rows")
        continue

    res = fisher_rxC_mc(table_df, B=20000, seed=1)

    results.append({
        "dataset": name,
        "groups_tested": ",".join(cols),
        "chi2_observed": res["chi2_observed"],
        "p_value": res["p_value"],
        "B": res["B"]
    })

# combine results
results_df = pd.DataFrame(results)

print(results_df)
