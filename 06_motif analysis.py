import scanpy as sc
import scirpy as ir
import matplotlib.pyplot as plt
import matplotlib
import os

# ---- Plot settings ---- #
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['savefig.format'] = 'pdf'

# ---- Function: Plot CDR3 motifs ---- #
def plot_majority_cdr3_motifs(adata, cc_levenshtein, output_dir="cdr3_motifs"):
    """
    Plots CDR3 sequence motifs for the most common length per chain
    in a given clonotype cluster (`cc_aa_levenshtein`).

    Parameters:
    - adata: AnnData with IR data (Scirpy schema)
    - cc_levenshtein: str, clonotype cluster ID
    - output_dir: str, path to save plots
    """

    os.makedirs(output_dir, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(12, 2))

    for chain, ax in zip(["VJ_1", "VDJ_1"], axs):
        # Get chain-level CDR3 sequences
        chain_seqs = ir.get.airr(adata, "junction_aa", chain)

        # Filter for specific cluster
        filtered_seqs = chain_seqs[adata.obs["cc_aa_levenshtein"] == cc_levenshtein]
        if filtered_seqs.empty:
            print(f" No {chain} sequences in cluster {cc_levenshtein}")
            continue

        # Identify dominant CDR3 length
        expected_length = filtered_seqs.str.len().value_counts().idxmax()
        mask = (
            (adata.obs["cc_aa_levenshtein"] == cc_levenshtein) &
            (chain_seqs.str.len() == expected_length)
        )

        if mask.sum() == 0:
            print(f"No {chain} sequences of length {expected_length} in cluster {cc_levenshtein}")
            continue

        # Plot motif
        ir.pl.logoplot_cdr3_motif(
            adata[mask],
            chains=chain,
            to_type="information",
            ax=ax
        )

    plt.tight_layout()
    plt.show()

# Load annotated dataset
adata = sc.read("path/gts_cd8_10x_ss2.h5ad")
ir.io.upgrade_schema(adata)  # Upgrade IR format if needed

share = (
    adata.obs.groupby('cc_aa_levenshtein')['donor_id']
    .nunique()          # number of unique donors for each cc_aa_levenshtein
    .gt(2)              # keep only those with > 2 donors
)

share_list = share[share].index.tolist()
print(share_list)


# Plot motifs for selected clusters (shared at least 3 donors)
share_list = ['36', '37', '48', '72', '108', '233']
for cluster_id in clusters:
    plot_majority_cdr3_motifs(adata, cc_levenshtein=share_list)
print(share_list)
