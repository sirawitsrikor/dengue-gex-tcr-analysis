import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

adata_a = adata[adata.obs["Clinical_phase"] == "Acute"].copy()

adata_a.obs["PE-A"] = pd.to_numeric(adata_a.obs["PE-A"], errors="coerce")
adata_a.obs["APC-A"] = pd.to_numeric(adata_a.obs["APC-A"], errors="coerce")

adata_a.obs["CTFI_log"] = np.log10(
    np.sqrt(adata_a.obs["PE-A"]**2 + adata_a.obs["APC-A"]**2) + 1
)

df_donor = (
    adata_a.obs[["donor_ID", "leiden", "CTFI_log"]]
    .dropna()
    .groupby(["donor_ID", "leiden"], as_index=False, observed=True)
    .agg(CTFI_median=("CTFI_log", "median"),
         n_cells=("CTFI_log", "size"))
    .query("n_cells >= 4")
)

# plot
sns.boxplot(data=df_donor, x="leiden", y="CTFI_median", showfliers=False)
sns.stripplot(data=df_donor, x="leiden", y="CTFI_median",
              color="black", size=6, jitter=0.15)

plt.xlabel("Cluster")
plt.ylabel("Donor-level median composite tetramer fluorescence (log)")
plt.tight_layout()
plt.show()

# mixed model
df = adata_a.obs[["donor_ID", "leiden", "CTFI_log"]].dropna().astype(str)
df["CTFI_log"] = adata_a.obs.loc[df.index, "CTFI_log"]

model = smf.mixedlm("CTFI_log ~ C(leiden)", data=df, groups=df["donor_ID"]).fit(reml=False)
print(model.summary())
