#after running run_ldsc.sh in the cluster, run this script for LDSC analyses and correlation heatmap
library(data.table)

ldsc_results_path <- "/gpfs/gibbs/project/polimanti/ejm96/immune/results/pairwise_ldsc/already_munged_sanity_11787120/rg_matrix.tsv"

ldsc_results <- fread(ldsc_results_path)

head(ldsc_results)

rg <- ldsc_results

### Correlation Heatmap
library(ggplot2)


# Convert values such as "nan" to NA, then calculate BH-FDR
rg[, `:=`(
  rg = suppressWarnings(as.numeric(rg)),
  p  = suppressWarnings(as.numeric(p))
)]

rg[, q_fdr := NA_real_]
valid <- which(!is.na(rg$p) & is.finite(rg$p))
rg$q_fdr[valid] <- p.adjust(rg$p[valid], method = "BH")

# Preserve the trait order from the pairwise output
traits <- unique(c(rg$p1, rg$p2))
rg[, `:=`(
  p1_order = match(p1, traits),
  p2_order = match(p2, traits)
)]

# In the pairwise output p1 precedes p2; displaying p1 as columns
# and p2 as rows gives the lower triangle.
lower <- rg[p1_order < p2_order]

lower[, label := fifelse(
  is.na(rg),
  "NA",
  sprintf("%.2f%s", rg, ifelse(!is.na(q_fdr) & q_fdr < 0.05, "*", ""))
)]

ggplot(
  lower,
  aes(
    x = factor(p1, levels = traits),
    y = factor(p2, levels = rev(traits)),
    fill = rg
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey85",
    name = expression(r[g])
  ) +
  coord_fixed() +
  labs(
    x = NULL,
    y = NULL,
    title = "FinnGen genetic correlations",
    subtitle = "* BH-FDR < 0.05"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


