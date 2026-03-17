library(dplyr)
library(pheatmap)
library(tibble)

files <- list(
Acute_TRA_VJ = "./heatmap_TRA_VJ_raw_counts_Acute.csv",
Acute_TRB_VJ = "./heatmap_TRB_VJ_raw_counts_Acute.csv",
Convalescence_TRA_VJ = "./heatmap_TRA_VJ_raw_counts_Con.csv",
Convalescence_TRB_VJ = "./heatmap_TRB_VJ_raw_counts_Con.csv"
)

for (name in names(files)) {

cat("\n==============================\n")
cat("Processing:", name, "\n")

df <- read.csv(files[[name]], row.names = 1, check.names = FALSE)

severity_cols <- intersect(c("AD","DF","DHF"), colnames(df))

zone_table <- df %>%
group_by(zone) %>%
summarise(across(all_of(severity_cols), sum), .groups = "drop") %>%
column_to_rownames("zone")

set.seed(1)
fisher_res <- fisher.test(zone_table, simulate.p.value = TRUE, B = 100000)
chi_res <- chisq.test(zone_table)

std_res <- chi_res$stdres

n <- sum(zone_table)
k <- min(nrow(zone_table)-1, ncol(zone_table)-1)
cramer_v <- sqrt(as.numeric(chi_res$statistic) / (n * k))

max_abs <- max(abs(std_res))

pheatmap(
std_res,
cluster_rows = FALSE,
cluster_cols = FALSE,
color = colorRampPalette(c("blue","white","red"))(100),
breaks = seq(-max_abs, max_abs, length.out = 101),
main = paste0(name,
"\nFisher p=", signif(fisher_res$p.value,3),
", Cramer's V=", round(cramer_v,3))
)

print(fisher_res)
cat("Cramer's V:", round(cramer_v,3), "\n")

}
