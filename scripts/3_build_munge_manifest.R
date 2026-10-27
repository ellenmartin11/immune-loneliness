#!/usr/bin/env Rscript

# Build the sample-size manifest for the LDSC munging step.
#
# Sample sizes are taken from the CLEANED files themselves (per-SNP N, N_CASE, N_CONTROL medians), which are the ground truth for what munging will see. 
# Curated paper-verified values in n_overrides.tsv take precedenceand are cross-checked against the file; any disagreement >25% is flagged, not silently resolved.
#
# Usage: Rscript 3_build_munge_manifest.R 
#
# Output: <out_dir>/munge_manifest.tsv with one row per phenotype:
#   phenotype, clean_path, trait_type, n_policy,
#   n_median, n_case_median, n_control_median, neff_median (from file),
#   N_total, N_case, N_control, Neff (values to use for munging),
#   munge_note, flag
#
# Neff convention (LDSC / GenomicSEM for case-control):
#   Neff = 4 * N_case * N_control / (N_case + N_control)
# (equals total N when cases and controls are balanced).

suppressPackageStartupMessages(library(data.table))

#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) < 3L) {
 # stop("Usage: build_munge_manifest.R <clean_manifest.tsv> <n_overrides.tsv> <out_dir>")
#}
clean_manifest_path <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/clean_manifest.tsv"
overrides_path      <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/n_overrides.tsv"
out_dir             <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/munging"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


norm <- function(x) gsub("[^a-z0-9]", "", tolower(trimws(x)))

neff_cc <- function(ncase, ncontrol) {
  4 * ncase * ncontrol / (ncase + ncontrol)
}

manifest <- fread(clean_manifest_path)
stopifnot(all(c("phenotype", "clean_path") %in% names(manifest)))
overrides <- fread(overrides_path, na.strings = c("NA", ""))
overrides[, pheno_key := norm(phenotype)]

message("Building N manifest for ", nrow(manifest), " phenotypes")

rows <- vector("list", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  trait <- manifest$phenotype[[i]]
  path  <- manifest$clean_path[[i]]
  flag  <- character(0)
  note  <- character(0)
  
  if (!file.exists(path)) {
    rows[[i]] <- data.table(
      phenotype = trait, clean_path = path, trait_type = NA_character_,
      n_policy = NA_character_, n_median = NA_real_, n_case_median = NA_real_,
      n_control_median = NA_real_, neff_median = NA_real_,
      N_total = NA_real_, N_case = NA_real_, N_control = NA_real_,
      Neff = NA_real_, munge_note = "clean file missing", flag = "FILE_NOT_FOUND"
    )
    next
  }
  
  # ---- Read sample-size columns from the cleaned file (chunked) ------------
  header <- names(fread(path, nrows = 0))
  size_cols <- intersect(c("N", "N_CASE", "N_CONTROL"), header)
  if (!length(size_cols)) {
    rows[[i]] <- data.table(
      phenotype = trait, clean_path = path, trait_type = NA_character_,
      n_policy = NA_character_, n_median = NA_real_, n_case_median = NA_real_,
      n_control_median = NA_real_, neff_median = NA_real_,
      N_total = NA_real_, N_case = NA_real_, N_control = NA_real_,
      Neff = NA_real_, munge_note = "no sample-size columns in clean file",
      flag = "NO_N_COLUMNS"
    )
    next
  }
  
  chunk_rows <- 5e6L
  acc <- list(); n_seen <- 0L; skip <- 0L
  repeat {
    chunk <- fread(path, select = size_cols, skip = skip, nrows = chunk_rows,
                   showProgress = FALSE)
    if (nrow(chunk) == 0L) break
    acc[[length(acc) + 1L]] <- chunk
    n_seen <- n_seen + nrow(chunk)
    skip <- skip + nrow(chunk)
    if (nrow(chunk) < chunk_rows) break
  }
  sz <- rbindlist(acc)
  rm(acc)
  
  med <- function(x) if (all(is.na(x))) NA_real_ else as.numeric(median(x, na.rm = TRUE))
  n_med    <- if ("N" %in% names(sz)) med(sz$N) else NA_real_
  nc_med   <- if ("N_CASE" %in% names(sz)) med(sz$N_CASE) else NA_real_
  nt_med   <- if ("N_CONTROL" %in% names(sz)) med(sz$N_CONTROL) else NA_real_
  neff_med <- if (!is.na(nc_med) && !is.na(nt_med)) med(neff_cc(sz$N_CASE, sz$N_CONTROL)) else NA_real_
  
  # Constant (fabricated) case/control columns: a single distinct value with
  # no variation is a sheet-filled placeholder, not per-SNP data.
  is_constant <- function(x) {
    !is.null(x) && length(unique(na.omit(x))) <= 1L
  }
  nc_constant <- "N_CASE" %in% names(sz) && is_constant(sz$N_CASE)
  nt_constant <- "N_CONTROL" %in% names(sz) && is_constant(sz$N_CONTROL)
  rm(sz)
  
  # ---- Merge curated overrides ---------------------------------------------
  key <- norm(trait)
  ov <- overrides[pheno_key == key]
  if (nrow(ov) > 1L) {
    flag <- c(flag, "DUPLICATE_OVERRIDE")
    ov <- ov[1L]
  }
  
  has_override <- nrow(ov) == 1L
  if (has_override) {
    trait_type <- ov$trait_type
    n_policy   <- ov$n_policy
    N_total    <- ov$N_total
    N_case     <- ov$N_case
    N_control  <- ov$N_control
    Neff       <- ov$Neff
    note       <- ov$note
  } else {
    # Infer trait type from the file: case/control columns present and
    # non-constant -> case-control; otherwise quantitative.
    if (!is.na(nc_med) && !is.na(nt_med) && !nc_constant && !nt_constant) {
      trait_type <- "case_control"
      n_policy   <- "ncase_ncontrol"
      N_case     <- nc_med; N_control <- nt_med
      N_total    <- nc_med + nt_med
      Neff       <- neff_med
    } else if (!is.na(n_med)) {
      trait_type <- "quantitative"
      n_policy   <- "per_snp_n"
      N_total    <- n_med; N_case <- NA_real_; N_control <- NA_real_
      Neff       <- NA_real_
    } else {
      trait_type <- NA_character_; n_policy <- NA_character_
      N_total <- N_case <- N_control <- Neff <- NA_real_
      flag <- c(flag, "NO_USABLE_N")
    }
    note <- "derived from cleaned file (no curated override)"
  }
  
  # ---- Cross-check file vs override / policy -------------------------------
  if (has_override) {
    if (n_policy == "per_snp_n") {
      if (is.na(n_med)) {
        flag <- c(flag, "OVERRIDE_WANTS_PER_SNP_N_BUT_FILE_HAS_NONE")
      } else if (!is.na(N_total) && abs(n_med - N_total) / N_total > 0.25) {
        flag <- c(flag, "FILE_N_VS_OVERRIDE_MISMATCH")
      }
      if ((nc_constant || nt_constant) && !is.na(nc_med)) {
        flag <- c(flag, "FABRICATED_CONSTANT_NCASE_NCONTROL")
        note <- paste(note, "| Drop constant N_CASE/N_CONTROL at munging; use per-SNP N.")
      }
    } else if (n_policy == "ncase_ncontrol") {
      if (!is.na(nc_med) && !is.na(N_case) && abs(nc_med - N_case) / N_case > 0.25) {
        flag <- c(flag, "FILE_NCASE_VS_OVERRIDE_MISMATCH")
      }
      if (!is.na(nt_med) && !is.na(N_control) && abs(nt_med - N_control) / N_control > 0.25) {
        flag <- c(flag, "FILE_NCONTROL_VS_OVERRIDE_MISMATCH")
      }
      if (nc_constant || nt_constant) {
        flag <- c(flag, "CONSTANT_NCASE_NCONTROL_IN_FILE")
        note <- paste(note, "| File N_CASE/N_CONTROL are constant; override values used.")
      }
      if (!is.na(n_med) && !is.na(Neff) && abs(n_med - Neff) / Neff > 0.25) {
        note <- paste(note, "| File N column is not LDSC Neff; do not pass it as --N-total.")
      }
    }
  } else {
    if (nc_constant || nt_constant) {
      flag <- c(flag, "CONSTANT_NCASE_NCONTROL_IN_FILE")
      note <- paste(note, "| Constant N_CASE/N_CONTROL look sheet-filled; verify against the source paper and add an override row.")
    }
  }
  
  rows[[i]] <- data.table(
    phenotype = trait, clean_path = path, trait_type = trait_type,
    n_policy = n_policy,
    n_median = n_med, n_case_median = nc_med, n_control_median = nt_med,
    neff_median = neff_med,
    N_total = N_total, N_case = N_case, N_control = N_control, Neff = Neff,
    munge_note = note, flag = if (length(flag)) paste(flag, collapse = ";") else ""
  )
  message(sprintf("[%d/%d] %s: policy=%s Neff=%s flag=%s",
                  i, nrow(manifest), trait, n_policy,
                  formatC(Neff, format = "f", digits = 0),
                  if (length(flag)) paste(flag, collapse = ";") else ""))
}

out <- rbindlist(rows, fill = TRUE)
out_path <- file.path(out_dir, "munge_manifest.tsv")
fwrite(out, out_path, sep = "\t", na = "NA")

n_flagged <- sum(out$flag != "")
message("Wrote ", out_path, " (", nrow(out), " traits, ", n_flagged, " flagged)")
if (n_flagged) {
  message("Flagged traits: ", paste(out$phenotype[out$flag != ""], collapse = "; "))
}
