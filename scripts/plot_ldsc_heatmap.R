#!/usr/bin/env Rscript

# Plot the lower triangle of the pairwise LDSC genetic-correlation matrix.
# A single asterisk marks pairs with Benjamini-Hochberg FDR < 0.05.
# In RStudio, set results_dir below and click Source.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

results_dir <- "/gpfs/gibbs/project/polimanti/ejm96/immune/results/pairwise_ldsc/already_munged_sanity_11787120"
output_prefix <- file.path(results_dir, "rg_heatmap_fdr05")

if (length(args) > 2L) stop("Use at most two command-line arguments: results_dir and output_prefix.")
results_file <- file.path(results_dir, "ldsc_all_pairs_labeled_fdr.tsv")
manifest_file <- file.path(results_dir, "input_manifest.tsv")

if (!file.exists(results_file)) stop("Pairwise results table not found: ", results_file)
if (!file.exists(manifest_file)) stop("Input manifest not found: ", manifest_file)

results <- fread(results_file, check.names = FALSE)
manifest <- fread(manifest_file, check.names = FALSE)
normalize <- function(x) gsub("[^a-z0-9]", "", tolower(x))
find_col <- function(nms, wanted) {
  hit <- which(normalize(nms) == normalize(wanted))
  if (length(hit)) nms[[hit[[1L]]]] else NULL
}

required_manifest <- c("phenotype", "source")
missing_manifest <- required_manifest[!vapply(required_manifest, function(x) !is.null(find_col(names(manifest), x)), logical(1))]
if (length(missing_manifest)) {
  stop("Manifest is missing column(s): ", paste(missing_manifest, collapse = ", "))
}
required_results <- c("phenotype_1", "source_1", "phenotype_2", "source_2", "rg", "p_fdr")
missing_results <- required_results[!vapply(required_results, function(x) !is.null(find_col(names(results), x)), logical(1))]
if (length(missing_results)) {
  stop("Labeled results table is missing column(s): ", paste(missing_results, collapse = ", "))
}

manifest[, phenotype := trimws(as.character(get(find_col(names(manifest), "phenotype"))))]
manifest[, source := trimws(as.character(get(find_col(names(manifest), "source"))))]
if ("ldsc_name" %in% names(manifest)) {
  manifest[, ldsc_name := as.character(ldsc_name)]
}
manifest[, plot_label := phenotype]
if (anyDuplicated(manifest$plot_label)) {
  manifest[, plot_label := paste0(phenotype, " [", source, "]")]
}
manifest[, trait_order := .I]

get_result_col <- function(name) {
  column <- find_col(names(results), name)
  if (is.null(column)) stop("Results table is missing column: ", name) #this returns the name of missing columns
  results[[column]]
}

  

plot_data <- data.table(
  phenotype_1 = as.character(get_result_col("phenotype_1")),
  source_1 = as.character(get_result_col("source_1")),
  phenotype_2 = as.character(get_result_col("phenotype_2")),
  source_2 = as.character(get_result_col("source_2")),
  rg = suppressWarnings(as.numeric(get_result_col("rg"))),
  p_fdr = suppressWarnings(as.numeric(get_result_col("p_fdr")))
)
plot_data <- merge(plot_data,
  manifest[, .(phenotype, source, plot_label, trait_order)],
  by.x = c("phenotype_1", "source_1"), by.y = c("phenotype", "source"),
  all.x = TRUE, sort = FALSE
)
setnames(plot_data, c("plot_label", "trait_order"), c("label_1", "order_1"))
plot_data <- merge(plot_data,
  manifest[, .(phenotype, source, plot_label, trait_order)],
  by.x = c("phenotype_2", "source_2"), by.y = c("phenotype", "source"),
  all.x = TRUE, sort = FALSE
)
setnames(plot_data, c("plot_label", "trait_order"), c("label_2", "order_2"))
if (anyNA(plot_data$label_1) || anyNA(plot_data$label_2)) {
  stop("Some result phenotypes could not be matched to input_manifest.tsv.")
}

# Reorient every pair to row >= column so only the lower triangle is plotted.
plot_data[, row_order := pmax(order_1, order_2)]
plot_data[, col_order := pmin(order_1, order_2)]
label_by_order <- manifest$plot_label
plot_data[, row_label := label_by_order[row_order]]
plot_data[, col_label := label_by_order[col_order]]
plot_data[, significance := ifelse(!is.na(p_fdr) & p_fdr < 0.05, "*", "")]

expected_pairs <- choose(nrow(manifest), 2L)
if (nrow(plot_data) != expected_pairs) {
  stop("Expected ", expected_pairs, " pairwise results, but found ", nrow(plot_data),
       ". This plot requires the full pairwise LDSC run.")
}

plot_data[, col_label := factor(col_label, levels = manifest$plot_label)]
plot_data[, row_label := factor(row_label, levels = rev(manifest$plot_label))]

p <- ggplot(plot_data, aes(x = col_label, y = row_label, fill = rg)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = significance), size = 4, na.rm = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "Genetic\ncorrelation (rg)"
  ) +
  coord_fixed() +
  labs(
    title = "Pairwise LDSC genetic correlations",
    subtitle = "Lower triangle shown; * indicates Benjamini–Hochberg FDR < 0.05",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(paste0(output_prefix, ".png"), p, width = 11, height = 10, dpi = 300)
ggsave(paste0(output_prefix, ".pdf"), p, width = 11, height = 10)
message("Saved heatmap: ", output_prefix, ".png and .pdf")
