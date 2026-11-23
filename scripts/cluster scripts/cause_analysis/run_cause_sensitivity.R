#!/usr/bin/env Rscript

# Run the 12 directional CAUSE analyses from 10_mr_cause.Rmd using 1e-5 as
# both the local LD-clumping threshold and CAUSE's exposure p-value threshold.
# Each completed pair is cached so a rerun can resume without refitting it.

suppressPackageStartupMessages({
  library(data.table)
  library(cause)
  library(ieugwasr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_cause_sensitivity.R <output_dir>")
}

out_dir <- args[[1]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(out_dir, "fit_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

P_THRESHOLD <- 1e-5
CLUMP_KB <- 10000
CLUMP_R2 <- 0.01
PARAM_SAMPLE_N <- 1000000L
SEED <- 100L

PLINK_BIN <- "/vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink"
LD_BFILE <- "/gpfs/gibbs/pi/polimanti/LDSC_1000G/g1000_eur"

GWAS_PATHS <- c(
  Loneliness = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/005_Loneliness_isolation_PanUKB_clean.tsv.gz",
  Triglycerides = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw/Triglycerides_METAL_combined.tsv.gz",
  T2D = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/006_Type_2_Diabetes_clean.tsv.gz",
  Depression = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/012_Depression_PGC_clean.tsv.gz"
)

pair_specs <- list(
  lone_to_tg = c("Loneliness", "Triglycerides"),
  tg_to_lone = c("Triglycerides", "Loneliness"),
  t2d_to_tg = c("T2D", "Triglycerides"),
  tg_to_t2d = c("Triglycerides", "T2D"),
  t2d_to_lone = c("T2D", "Loneliness"),
  lone_to_t2d = c("Loneliness", "T2D"),
  t2d_to_dep = c("T2D", "Depression"),
  dep_to_t2d = c("Depression", "T2D"),
  dep_to_tg = c("Depression", "Triglycerides"),
  tg_to_dep = c("Triglycerides", "Depression"),
  lone_to_dep = c("Loneliness", "Depression"),
  dep_to_lone = c("Depression", "Loneliness")
)

required_sumstats <- c("SNP", "A1", "A2", "BETA", "SE")
missing_paths <- GWAS_PATHS[!file.exists(GWAS_PATHS)]
if (length(missing_paths)) {
  stop("GWAS input file(s) not found: ", paste(missing_paths, collapse = ", "))
}
if (!file.exists(PLINK_BIN)) stop("PLINK executable not found: ", PLINK_BIN)
if (!all(file.exists(paste0(LD_BFILE, c(".bed", ".bim", ".fam"))))) {
  stop("PLINK reference panel is missing one or more .bed/.bim/.fam files: ", LD_BFILE)
}

loo_version <- as.character(packageVersion("loo"))
if (!identical(loo_version, "2.9.0")) {
  stop("This CAUSE setup requires loo 2.9.0; loaded version is ", loo_version,
       ". Load the R library where loo 2.9.0 and CAUSE are installed.")
}

message("CAUSE version: ", as.character(packageVersion("cause")))
message("loo version: ", loo_version)
message("Clumping p threshold: ", P_THRESHOLD)
message("CAUSE pval_thresh: ", P_THRESHOLD)

read_gwas <- function(trait, path) {
  headers <- names(fread(path, nrows = 0, showProgress = FALSE))
  missing <- setdiff(required_sumstats, headers)
  if (length(missing)) {
    stop(trait, " is missing required columns: ", paste(missing, collapse = ", "))
  }

  fread(path, select = required_sumstats, showProgress = TRUE)
}

gwas <- lapply(names(GWAS_PATHS), function(trait) {
  read_gwas(trait, GWAS_PATHS[[trait]])
})
names(gwas) <- names(GWAS_PATHS)

make_pair_data <- function(trait1, trait2) {
  x <- cause::gwas_merge(
    gwas[[trait1]], gwas[[trait2]],
    snp_name_cols = c("SNP", "SNP"),
    beta_hat_cols = c("BETA", "BETA"),
    se_cols = c("SE", "SE"),
    A1_cols = c("A1", "A1"),
    A2_cols = c("A2", "A2")
  )
  x <- as.data.table(x)

  required <- c("snp", "beta_hat_1", "seb1", "beta_hat_2", "seb2")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("CAUSE gwas_merge output is missing: ", paste(missing, collapse = ", "))
  }

  x[, `:=`(
    snp = as.character(snp),
    beta_hat_1 = as.numeric(beta_hat_1),
    seb1 = as.numeric(seb1),
    beta_hat_2 = as.numeric(beta_hat_2),
    seb2 = as.numeric(seb2)
  )]

  x <- x[
    !is.na(snp) & nzchar(snp) &
      is.finite(beta_hat_1) & is.finite(seb1) & seb1 > 0 &
      is.finite(beta_hat_2) & is.finite(seb2) & seb2 > 0
  ]

  duplicate_snps <- x[, .N, by = snp][N > 1L, snp]
  if (length(duplicate_snps)) x <- x[!snp %chin% duplicate_snps]
  if (nrow(x) < 5L) stop("Fewer than five usable shared SNPs for ", trait1, " -> ", trait2)

  x[, p1 := 2 * pnorm(-abs(beta_hat_1 / seb1))]
  x[!is.finite(p1), p1 := 1]
  x
}

extract_parameter_summary <- function(qmat, model_name) {
  out <- list()
  for (parameter in c("gamma", "eta", "q")) {
    values <- if (is.matrix(qmat) && parameter %in% colnames(qmat)) {
      as.numeric(qmat[, parameter])
    } else {
      rep(NA_real_, 3L)
    }
    if (length(values) < 3L) values <- rep(NA_real_, 3L)
    out[[paste(model_name, parameter, "median", sep = "_")]] <- values[1]
    out[[paste(model_name, parameter, "ci_low", sep = "_")]] <- values[2]
    out[[paste(model_name, parameter, "ci_high", sep = "_")]] <- values[3]
  }
  out
}

summarise_fit <- function(fit, analysis, trait1, trait2, metadata) {
  s <- summary(fit, ci_size = 0.95, digits = 8)
  q <- s[["quants"]]
  if (!is.list(q) || length(q) < 2L) {
    stop("Unexpected CAUSE quantile summary for ", analysis)
  }

  row <- c(
    list(
      analysis = analysis,
      trait1 = trait1,
      trait2 = trait2,
      p_threshold = P_THRESHOLD,
      n_shared_usable = metadata$n_shared_usable,
      n_clumped = metadata$n_clumped,
      n_cause_variants = if (!is.null(fit$data)) nrow(fit$data) else NA_integer_,
      model_comparison_z = as.numeric(s[["z"]]),
      model_comparison_p = as.numeric(s[["p"]]),
      causal_model_favored = is.finite(s[["z"]]) && is.finite(s[["p"]]) &&
        s[["z"]] < 0 && s[["p"]] < 0.05,
      run_status = "completed"
    ),
    extract_parameter_summary(q[[1]], "sharing"),
    extract_parameter_summary(q[[2]], "causal")
  )

  elpd <- as.data.table(fit[["elpd"]])
  if (nrow(elpd)) {
    for (i in seq_len(nrow(elpd))) {
      model1 <- if ("model1" %in% names(elpd)) as.character(elpd$model1[i]) else "model"
      model2 <- if ("model2" %in% names(elpd)) as.character(elpd$model2[i]) else paste0("row", i)
      cmp <- gsub("[^A-Za-z0-9]+", "_", paste(model1, "vs", model2, sep = "_"))

      for (field in setdiff(names(elpd), c("model1", "model2"))) {
        value <- elpd[[field]][i]
        if (is.list(value)) value <- value[[1]]
        row[[paste("elpd", cmp, field, sep = "_")]] <- value
      }
    }
  }

  as.data.table(row)
}

cause_fits_p1e5 <- list()
cause_params_p1e5 <- list()
cause_top_vars_p1e5 <- list()
cause_metadata_p1e5 <- list()
summary_rows <- list()
status_rows <- list()

for (analysis in names(pair_specs)) {
  traits <- pair_specs[[analysis]]
  trait1 <- traits[1]
  trait2 <- traits[2]
  message("\n--- ", analysis, ": ", trait1, " -> ", trait2, " ---")

  fit_path <- file.path(cache_dir, paste0("cause_", analysis, "_p1e-5.rds"))
  params_path <- file.path(cache_dir, paste0("params_", analysis, "_p1e-5.rds"))
  vars_path <- file.path(cache_dir, paste0("top_vars_", analysis, "_p1e-5.rds"))
  meta_path <- file.path(cache_dir, paste0("metadata_", analysis, "_p1e-5.rds"))

  result <- tryCatch({
    if (file.exists(fit_path)) {
      message("Reusing cached CAUSE fit: ", fit_path)
      fit <- readRDS(fit_path)
      metadata <- if (file.exists(meta_path)) readRDS(meta_path) else list(
        n_shared_usable = NA_integer_, n_clumped = NA_integer_
      )
      params <- if (file.exists(params_path)) readRDS(params_path) else NULL
      top_vars <- if (file.exists(vars_path)) readRDS(vars_path) else character()
    } else {
      X <- make_pair_data(trait1, trait2)
      message("Usable shared SNPs: ", format(nrow(X), big.mark = ","))

      if (file.exists(params_path)) {
        message("Reusing cached parameter estimates.")
        params <- readRDS(params_path)
      } else {
        set.seed(SEED)
        n_sample <- min(PARAM_SAMPLE_N, nrow(X))
        varlist <- sample(X$snp, size = n_sample, replace = FALSE)
        message("Estimating CAUSE nuisance parameters using ",
                format(n_sample, big.mark = ","), " variants.")
        params <- cause::est_cause_params(X, varlist)
        saveRDS(params, params_path)
      }

      if (file.exists(vars_path)) {
        message("Reusing cached clumped variant list.")
        top_vars <- readRDS(vars_path)
      } else {
        ld_input <- data.frame(
          rsid = X$snp,
          pval = X$p1,
          id = analysis
        )
        ld_input <- ld_input[
          !is.na(ld_input$rsid) & nzchar(ld_input$rsid) &
            is.finite(ld_input$pval) & ld_input$pval <= P_THRESHOLD,
        ]
        ld_input <- ld_input[!duplicated(ld_input$rsid), , drop = FALSE]
        if (!nrow(ld_input)) stop("No variants pass clumping p threshold for ", analysis)

        clumped <- ieugwasr::ld_clump(
          dat = ld_input,
          clump_kb = CLUMP_KB,
          clump_r2 = CLUMP_R2,
          clump_p = P_THRESHOLD,
          bfile = LD_BFILE,
          plink_bin = PLINK_BIN
        )
        top_vars <- unique(as.character(clumped$rsid))
        saveRDS(top_vars, vars_path)
      }

      if (!length(top_vars)) stop("No clumped variants remain for ", analysis)
      metadata <- list(
        analysis = analysis,
        trait1 = trait1,
        trait2 = trait2,
        p_threshold = P_THRESHOLD,
        clump_kb = CLUMP_KB,
        clump_r2 = CLUMP_R2,
        n_shared_usable = nrow(X),
        n_clumped = length(top_vars)
      )

      set.seed(SEED)
      fit <- cause::cause(
        X = X,
        variants = top_vars,
        param_ests = params,
        pval_thresh = P_THRESHOLD
      )
      saveRDS(fit, fit_path)
      saveRDS(metadata, meta_path)
      rm(X)
      gc(verbose = FALSE)
    }

    list(ok = TRUE, fit = fit, params = params, top_vars = top_vars, metadata = metadata)
  }, error = function(e) {
    message("FAILED ", analysis, ": ", conditionMessage(e))
    list(ok = FALSE, error = conditionMessage(e))
  })

  if (isTRUE(result$ok)) {
    cause_fits_p1e5[[analysis]] <- result$fit
    cause_params_p1e5[[analysis]] <- result$params
    cause_top_vars_p1e5[[analysis]] <- result$top_vars
    cause_metadata_p1e5[[analysis]] <- result$metadata
    summary_rows[[analysis]] <- summarise_fit(
      result$fit, analysis, trait1, trait2, result$metadata
    )
    status_rows[[analysis]] <- data.table(
      analysis = analysis, trait1 = trait1, trait2 = trait2,
      run_status = "completed", error = NA_character_
    )
    message("Completed ", analysis)
  } else {
    summary_rows[[analysis]] <- data.table(
      analysis = analysis, trait1 = trait1, trait2 = trait2,
      p_threshold = P_THRESHOLD, run_status = "failed",
      error = result$error
    )
    status_rows[[analysis]] <- data.table(
      analysis = analysis, trait1 = trait1, trait2 = trait2,
      run_status = "failed", error = result$error
    )
  }

  # Refresh outputs after each pair so completed work is retained on interruption.
  cause_summary_p1e5 <- rbindlist(summary_rows, fill = TRUE, use.names = TRUE)
  cause_status_p1e5 <- rbindlist(status_rows, fill = TRUE, use.names = TRUE)
  fwrite(cause_summary_p1e5, file.path(out_dir, "cause_sensitivity_p1e-5_all_outputs.csv"))
  fwrite(cause_status_p1e5, file.path(out_dir, "cause_sensitivity_p1e-5_status.csv"))
  save(
    cause_fits_p1e5, cause_params_p1e5, cause_top_vars_p1e5,
    cause_metadata_p1e5, cause_summary_p1e5, cause_status_p1e5,
    file = file.path(out_dir, "cause_sensitivity_p1e-5_all_results.RData")
  )
}

message("\nCAUSE 1e-5 sensitivity run finished. Outputs: ", out_dir)
