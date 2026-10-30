#!/usr/bin/env Rscript
# 4_diag_munge.R — inspect columns of traits whose munging failed
# Usage: Rscript diag_munge_fails.R [clean_dir] [n_rows]
# Default clean_dir: /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean

suppressPackageStartupMessages(library(data.table))

clean_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(clean_dir)) clean_dir <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean"
n_rows <- as.integer(commandArgs(trailingOnly = TRUE)[2])
if (is.na(n_rows)) n_rows <- 20000

files <- c("030_Tyrosine")

for (f in files) {
  p <- file.path(clean_dir, paste0(f, "_clean.tsv.gz"))
  cat("\n=====", f, "=====\n")
  if (!file.exists(p)) { cat("  FILE MISSING\n"); next }
  d <- fread(p, nrows = n_rows, showProgress = FALSE, na.strings = c("NA", ""))
  cat("columns:", paste(names(d), collapse = ", "), "\n")
  na_frac <- vapply(d, function(x) round(mean(is.na(x)), 3), numeric(1))
  print(data.table(column = names(d), na_frac = na_frac)[na_frac > 0 | TRUE])
  # signed stat sanity
  if ("Z" %in% names(d)) {
    z <- d$Z[is.finite(d$Z)]
    if (length(z)) cat(sprintf("Z: n=%d median=%.4f mean=%.4f min=%.2f max=%.2f\n",
                               length(z), median(z), mean(z), min(z), max(z)))
    else cat("Z: ZERO finite values\n")
  }
  if ("P" %in% names(d)) cat(sprintf("P range: [%s, %s]\n",
                                     format(min(d$P, na.rm = TRUE)),
                                     format(max(d$P, na.rm = TRUE))))
  if ("INFO" %in% names(d)) {
    ii <- d$INFO[is.finite(d$INFO)]
    cat(sprintf("INFO finite: %d/%d\n", length(ii), nrow(d)))
  }
}
