#!/usr/bin/env Rscript

# Combine completed loneliness-versus-trait results and apply BH-FDR once.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: aggregate_metal_primary_ldsc.R <munged_dir> <results_dir>")
}
munged_dir <- args[[1L]]
results_dir <- args[[2L]]
manifest_path <- file.path(munged_dir, "ldsc_inputs_manifest.tsv")
if (!file.exists(manifest_path)) stop("Missing prepared-input manifest: ", manifest_path)
traits <- fread(manifest_path, sep = "\t", showProgress = FALSE)
if (!all(c("phenotype", "source", "origin", "path") %in% names(traits))) {
  stop("Prepared-input manifest is missing phenotype/source/origin/path columns.")
}
anchor_i <- which(traits$phenotype == "Loneliness isolation PanUKB")
if (length(anchor_i) != 1L) stop("Expected exactly one Loneliness isolation PanUKB row in the manifest.")
others <- setdiff(seq_len(nrow(traits)), anchor_i)
safe_label <- function(x) {
  value <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  value <- gsub("^_+|_+$", "", value)
  ifelse(is.na(value) | !nzchar(value), "trait", value)
}

result_paths <- vapply(others, function(i) {
  file.path(results_dir, sprintf("loneliness_%03d_%s.result.tsv", i, safe_label(traits$phenotype[[i]])))
}, character(1))
missing <- result_paths[!file.exists(result_paths)]
if (length(missing)) {
  stop(length(missing), " of ", length(result_paths), " pair results are missing. ",
       "Check the Slurm array task logs and increase its array range if needed. First missing: ",
       paste(basename(head(missing, 5L)), collapse = ", "))
}

results <- rbindlist(lapply(result_paths, fread), use.names = TRUE, fill = TRUE)
if (nrow(results) != length(others)) stop("Unexpected number of result rows.")
results[, p_fdr := p.adjust(as.numeric(p), method = "BH")]
setorder(results, p_fdr, p)
output_path <- file.path(results_dir, "loneliness_isolation_rg_all_traits_fdr.tsv")
fwrite(results, output_path, sep = "\t", na = "")
message("Wrote ", nrow(results), " FDR-corrected correlations to: ", output_path)
