#!/usr/bin/env Rscript
# qc_munge_logs.R -- sanity-check table from munge_sumstats.py logs
#
# Parses the Metadata block of each munged/<trait>.log into one row per trait:
#   SNPs read, SNPs remaining after filtering, median Z, mean chi^2, lambda GC,
#   max chi^2, genome-wide significant SNPs.
# Use it to spot low-polygenic traits before running LDSR.
#
# Usage (on the cluster):
#   Rscript qc_munge_logs.R                       # MUNGED_DIR env or "munged"
#   Rscript qc_munge_logs.R /path/to/munged       # or pass dir as arg 1
#
# Env vars:
#   MUNGED_DIR   dir with *.sumstats.gz + *.log   (default: munged/ under cwd)

suppressPackageStartupMessages(library(data.table))

munged_dir <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/munged"

logs <- sort(Sys.glob(file.path(munged_dir, "*.log")))
if (length(logs) == 0) stop("No .log files found in ", munged_dir)

# Extract the number captured by the pattern's single capture group
grab_num <- function(txt, pat) {
  m <- regmatches(txt, regexpr(pat, txt))
  if (length(m) == 0) return(NA_real_)
  as.numeric(sub(pat, "\\1", m))
}

rows <- lapply(logs, function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  data.table(
    trait       = sub("\\.log$", "", basename(f)),
    snps_read   = grab_num(txt, "Read ([0-9.]+) SNPs from --sumstats"),
    snps_remain = grab_num(txt, "duplicated rs numbers \\(([0-9.]+) SNPs remain"),
    median_z    = grab_num(txt, "Median value of SIGNED_SUMSTATS was (-?[0-9.]+)"),
    mean_chi2   = grab_num(txt, "Mean chi\\^2 = (-?[0-9.]+)"),
    lambda_gc   = grab_num(txt, "Lambda GC = (-?[0-9.]+)"),
    max_chi2    = grab_num(txt, "Max chi\\^2 = (-?[0-9.]+)"),
    gw_sig_snps = grab_num(txt, "([0-9.]+) Genome-wide significant SNPs"),
    error       = grepl("ERROR|ValueError", txt)
  )
})

out <- rbindlist(rows, fill = TRUE)
setkey(out, trait)
cat("\n== Munge QC summary:", nrow(out), "traits ==\n\n")
print(out, row.names = FALSE)

cat("\nTraits with errors in log:", sum(out$error, na.rm = TRUE), "\n")
cat("Traits with mean chi^2 < 1 (low polygenic signal; OK for LDSC, worth knowing):",
    sum(out$mean_chi2 < 1, na.rm = TRUE), "\n")

fwrite(out, file.path(munged_dir, "munge_qc_summary.tsv"), sep = "\t")
cat("\nWritten:", file.path(munged_dir, "munge_qc_summary.tsv"), "\n")
