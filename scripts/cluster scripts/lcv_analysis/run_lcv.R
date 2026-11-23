#!/usr/bin/env Rscript

# Run six pre-aligned, position-ordered LCV analyses and save each result.
# The summary CSV reports only the specified input direction for each pair.
# Expected input RDS: a named list with T2D_to_TG, TG_to_Loneliness,
# Depression_to_T2D, T2D_to_Loneliness, Loneliness_to_Depression, and
# Depression_to_TG. Each table must contain SNP, Z_1, Z_2, L2, N_1, and N_2.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: run_lcv.R <aligned_ordered_inputs.rds> <output_dir> <LCV_R_dir>")
}

input_rds <- args[[1]]
output_dir <- args[[2]]
lcv_r_dir <- args[[3]]

if (!file.exists(input_rds)) stop("Input RDS not found: ", input_rds)
if (!file.exists(file.path(lcv_r_dir, "RunLCV.R"))) {
  stop("RunLCV.R not found in: ", lcv_r_dir)
}
if (!file.exists(file.path(lcv_r_dir, "MomentFunctions.R"))) {
  stop("MomentFunctions.R not found in: ", lcv_r_dir)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- readRDS(input_rds)

required_pairs <- c(
  "T2D_to_TG", "TG_to_Loneliness", "Depression_to_T2D",
  "T2D_to_Loneliness", "Loneliness_to_Depression", "Depression_to_TG"
)
missing_pairs <- setdiff(required_pairs, names(inputs))
if (length(missing_pairs)) {
  stop("Input RDS is missing pair(s): ", paste(missing_pairs, collapse = ", "))
}

# The supplied pair tables must already have aligned alleles and be sorted in
# genomic order. LCV uses a block jackknife, so preserving that order matters.
run_one_lcv <- function(dat, pair_name, lcv_r_dir) {
  dat <- as.data.table(dat)
  # Accept common names for the per-SNP LD score, then standardize to L2.
  ld_candidates <- names(dat)[tolower(names(dat)) %chin% c(
    "l2", "ell", "ldscore", "ld_score", "l2.x", "ell.x",
    "ldscore.x", "ld_score.x"
  )]
  if (!"L2" %in% names(dat) && length(ld_candidates)) {
    setnames(dat, ld_candidates[1], "L2")
  }

  required_cols <- c("SNP", "Z_1", "Z_2", "L2", "N_1", "N_2")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols)) {
    stop(pair_name, " is missing columns: ", paste(missing_cols, collapse = ", "),
         ". Available columns: ", paste(names(dat), collapse = ", "),
         ". LCV requires a per-SNP LD score; if it was dropped during ",
         "alignment, merge it back by SNP before saving the input RDS.")
  }
  if (!nrow(dat)) stop(pair_name, " has no variants.")
  if (anyDuplicated(dat$SNP)) stop(pair_name, " has duplicated SNP IDs.")
  if (any(!is.finite(dat$Z_1)) || any(!is.finite(dat$Z_2)) ||
      any(!is.finite(dat$L2)) || any(dat$L2 < 0)) {
    stop(pair_name, " contains missing/non-finite Z scores or LD scores.")
  }

  n1 <- median(as.numeric(dat$N_1), na.rm = TRUE)
  n2 <- median(as.numeric(dat$N_2), na.rm = TRUE)
  if (!is.finite(n1) || !is.finite(n2) || n1 <= 0 || n2 <= 0) {
    stop(pair_name, " has invalid sample-size values.")
  }

  message(pair_name, ": ", format(nrow(dat), big.mark = ","),
          " variants; LCV n.1=", round(n1), ", n.2=", round(n2))

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(lcv_r_dir)
  source("RunLCV.R", local = .GlobalEnv)

  RunLCV(
    ell = dat$L2,
    z.1 = dat$Z_1,
    z.2 = dat$Z_2,
    no.blocks = 100,
    crosstrait.intercept = 1,
    ldsc.intercept = 1,
    n.1 = n1,
    n.2 = n2
  )
}

summarise_lcv <- function(result, pair_name, n_variants) {
  pair_traits <- list(
    T2D_to_TG = c("T2D", "TG"),
    TG_to_Loneliness = c("TG", "Loneliness"),
    Depression_to_T2D = c("Depression", "T2D"),
    T2D_to_Loneliness = c("T2D", "Loneliness"),
    Loneliness_to_Depression = c("Loneliness", "Depression"),
    Depression_to_TG = c("Depression", "TG")
  )[[pair_name]]

  if (is.null(pair_traits)) stop("No trait-direction labels defined for ", pair_name)

  gcp <- as.numeric(result$gcp.pm)
  gcp_se <- as.numeric(result$gcp.pse)
  common <- data.table(
    pair_id = pair_name,
    n_variants = n_variants,
    GCP = gcp,
    GCP_SE = gcp_se,
    GCP_CI_low_approx = gcp - 1.96 * gcp_se,
    GCP_CI_high_approx = gcp + 1.96 * gcp_se,
    GCP_Z = as.numeric(result$zscore),
    GCP_P = as.numeric(result$pval.gcpzero.2tailed),
    genetic_correlation = as.numeric(result$rho.est),
    genetic_correlation_SE = as.numeric(result$rho.err),
    h2_Z_trait1 = as.numeric(result$h2.zscore[1]),
    h2_Z_trait2 = as.numeric(result$h2.zscore[2]),
    p_fully_causal_trait1_to_trait2 = as.numeric(result$pval.fullycausal[1]),
    p_fully_causal_trait2_to_trait1 = as.numeric(result$pval.fullycausal[2])
  )

  summary_row <- copy(common)
  summary_row[, `:=`(
    trait1 = pair_traits[1],
    trait2 = pair_traits[2],
    comparison = paste(pair_traits[1], "→", pair_traits[2])
  )]
  summary_row
}

results <- list()
summary_rows <- list()

for (pair_name in required_pairs) {
  pair_rds <- file.path(output_dir, paste0(pair_name, "_lcv.rds"))
  dat <- as.data.table(inputs[[pair_name]])

  if (file.exists(pair_rds)) {
    message("Reusing saved result: ", pair_rds)
    result <- readRDS(pair_rds)
  } else {
    result <- run_one_lcv(dat, pair_name, lcv_r_dir)
    saveRDS(result, pair_rds)
    message("Saved pair result: ", pair_rds)
  }

  results[[pair_name]] <- result
  summary_rows[[pair_name]] <- summarise_lcv(result, pair_name, nrow(dat))
}

summary_table <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
pair_tests <- unique(summary_table[, .(pair_id, GCP_P)])
pair_tests[, GCP_P_FDR := p.adjust(GCP_P, method = "BH")]
summary_table <- merge(summary_table, pair_tests, by = c("pair_id", "GCP_P"),
                       all.x = TRUE, sort = FALSE)
setcolorder(summary_table, c("pair_id", "comparison", "trait1", "trait2",
                             setdiff(names(summary_table),
                                     c("pair_id", "comparison", "trait1", "trait2"))))

fwrite(summary_table, file.path(output_dir, "lcv_results_summary.csv"))
saveRDS(results, file.path(output_dir, "lcv_results_all.rds"))

print(summary_table)
message("All LCV results saved in: ", output_dir)
