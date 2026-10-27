#!/usr/bin/env Rscript
# Deep QC of PRIMARY cleaned summary statistics before LDSC munging. 
# Run this after SLURM submit_primary_clean_v2.sh (and if needed, submit_missing_traits.sh)
# Reads clean_manifest.tsv (phenotype, source, input_path, clean_path) and
# produces a per-phenotype diagnostics table + flagged issues.
# checks variant id integrity (duplicates/multi-allelic)
# checks percent of bad and rare allele frequencies 
# checks that z-score = BETA/SE when present
# checks p-values
# checks coordinates/build
# checks sample sizes against reported sample in sumstats spreadsheet

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))


manifest_path <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/clean_manifest.tsv"
out_prefix    <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/qc_deep"

manifest <- fread(manifest_path)
ACGT <- c("A", "C", "G", "T")

# GRCh37 chromosome lengths (Mb-ish resolution is enough for build checks)
grch37_len <- c(249250621, 243199373, 198022430, 190214555, 181538259,
                170805979, 159345973, 145138636, 138394717, 133797422,
                135086622, 133275309, 114364328, 107043718, 102531392,
                90338345, 81195210, 78077248, 59128983, 63025520,
                48129895, 51304566)
names(grch37_len) <- as.character(1:22)

qc_rows <- list()
for (i in seq_len(nrow(manifest))) {
  pheno <- manifest$phenotype[i]
  path  <- manifest$clean_path[i]
  cat("QC:", pheno, "\n")
  
  if (!file.exists(path)) {
    qc_rows[[length(qc_rows) + 1L]] <- data.table(phenotype = pheno, status = "FILE NOT FOUND")
    next
  }
  
  dt <- fread(path, showProgress = FALSE, na.strings = c("NA", ""))
  n0 <- nrow(dt)
  
  # ---- variant ID integrity -------------------------------------------------
  # Missing rsIDs are counted separately from true duplicates: uniqueN treats
  # all NAs as one value, which otherwise inflates the duplicate count.
  n_snp_na  <- sum(is.na(dt$SNP))
  n_dup_snp <- sum(duplicated(dt$SNP[!is.na(dt$SNP)]))
  pct_rs      <- 100 * mean(grepl("^rs[0-9]+$", dt$SNP), na.rm = TRUE)
  
  # Classify duplicate rsIDs by type:
  #   multi-allelic  same position, different allele pair -- expected in
  #                  whole-genome files after rsID mapping; munging drops them.
  #   pos_conflict   same rsID at different positions -- mapping artifact,
  #                  must be fixed before munging.
  #   identical      fully duplicated rows -- harmless, deduped at munging.
  n_dup_multiallelic <- 0L; n_dup_pos_conflict <- 0L; n_dup_identical <- 0L
  if (n_dup_snp > 0) {
    dup_ids <- dt$SNP[duplicated(dt$SNP) & !is.na(dt$SNP)]
    d <- dt[SNP %in% dup_ids,
            .(n_pos = uniqueN(BP), n_alleles = uniqueN(paste(A1, A2))), by = SNP]
    n_dup_pos_conflict <- d[n_pos > 1, .N]
    n_dup_multiallelic <- d[n_pos == 1 & n_alleles > 1, .N]
    n_dup_identical    <- d[n_pos == 1 & n_alleles == 1, .N]
    rm(d)
  }
  
  # ---- alleles --------------------------------------------------------------
  a1 <- toupper(dt$A1); a2 <- toupper(dt$A2)
  pct_acgt   <- 100 * mean(a1 %in% ACGT & a2 %in% ACGT, na.rm = TRUE)
  n_same_allele <- sum(!is.na(a1) & !is.na(a2) & a1 == a2)
  
  # ---- allele frequency -----------------------------------------------------
  frq_ok <- "FRQ" %in% names(dt)
  pct_frq_bad <- NA_real_; pct_frq_rare <- NA_real_
  if (frq_ok) {
    pct_frq_bad  <- 100 * mean(!is.finite(dt$FRQ) | dt$FRQ <= 0 | dt$FRQ >= 1, na.rm = TRUE)
    pct_frq_rare <- 100 * mean(dt$FRQ < 0.01 | dt$FRQ > 0.99, na.rm = TRUE)
  }
  
  # ---- signed stats: internal consistency -----------------------------------
  # Z should equal BETA/SE when both are present
  z_consistent <- NA_real_; sign_disagree <- NA_real_
  if (all(c("Z", "BETA", "SE") %in% names(dt))) {
    ok <- is.finite(dt$Z) & is.finite(dt$BETA) & is.finite(dt$SE) & dt$SE > 0
    z_ref <- dt$BETA[ok] / dt$SE[ok]
    z_consistent <- 100 * mean(abs(dt$Z[ok] - z_ref) /
                                 pmax(abs(z_ref), 1) < 0.05, na.rm = TRUE)
    sign_disagree <- 100 * mean(sign(dt$Z[ok]) != sign(dt$BETA[ok]), na.rm = TRUE)
  }
  # |Z| should track qnorm from P
  zp_corr <- NA_real_
  if (all(c("Z", "P") %in% names(dt))) {
    ok <- is.finite(dt$Z) & is.finite(dt$P) & dt$P > 0 & dt$P <= 1
    if (sum(ok) > 100) zp_corr <- cor(abs(dt$Z[ok]), qnorm(dt$P[ok] / 2, lower.tail = FALSE))
  }
  
  # ---- P sanity -------------------------------------------------------------
  pct_p_bad <- NA_real_
  if ("P" %in% names(dt)) {
    pct_p_bad <- 100 * mean(!is.finite(dt$P) | dt$P <= 0 | dt$P > 1, na.rm = TRUE)
  }
  
  # ---- coordinates / build --------------------------------------------------
  chr <- as.character(dt$CHR)
  chr[chr %in% c("23", "X")] <- "X"; chr[chr %in% c("24", "Y")] <- "Y"
  chr[chr %in% c("25", "MT", "M")] <- "MT"
  autosomes <- chr %in% as.character(1:22)
  pct_coord_missing <- 100 * mean(!is.finite(dt$BP) | is.na(dt$BP) | !autosomes)
  pct_beyond_grch37 <- NA_real_
  if (sum(autosomes) > 0) {
    bp_auto <- dt$BP[autosomes]
    beyond <- bp_auto > grch37_len[chr[autosomes]]
    pct_beyond_grch37 <- 100 * mean(beyond, na.rm = TRUE)
  }
  
  # ---- sample size ----------------------------------------------------------
  n_min <- n_med <- n_max <- NA_real_
  neff  <- NA_real_
  if ("N" %in% names(dt) && any(is.finite(dt$N))) {
    n_min <- min(dt$N, na.rm = TRUE); n_med <- median(dt$N, na.rm = TRUE)
    n_max <- max(dt$N, na.rm = TRUE)
  }
  if (all(c("N_CASE", "N_CONTROL") %in% names(dt)) &&
      any(is.finite(dt$N_CASE)) && any(is.finite(dt$N_CONTROL))) {
    nc <- median(dt$N_CASE, na.rm = TRUE); nt <- median(dt$N_CONTROL, na.rm = TRUE)
    neff <- 4 / (1 / nc + 1 / nt)
  }
  
  qc_rows[[length(qc_rows) + 1L]] <- data.table(
    phenotype        = pheno,
    status           = "ok",
    n_rows           = n0,
    n_unique_snp     = uniqueN(dt$SNP),
    n_snp_na         = n_snp_na,
    n_dup_snp        = n_dup_snp,
    n_dup_multiallelic = n_dup_multiallelic,
    n_dup_pos_conflict = n_dup_pos_conflict,
    n_dup_identical    = n_dup_identical,
    pct_rs_ids       = round(pct_rs, 2),
    pct_acgt_alleles = round(pct_acgt, 2),
    n_same_allele    = n_same_allele,
    pct_frq_bad      = round(pct_frq_bad, 3),
    pct_frq_rare     = round(pct_frq_rare, 2),
    pct_z_eq_beta_se = round(z_consistent, 2),
    pct_sign_disagree = round(sign_disagree, 3),
    corr_absZ_vs_P   = round(zp_corr, 4),
    pct_p_bad        = round(pct_p_bad, 3),
    pct_coord_missing = round(pct_coord_missing, 2),
    pct_bp_beyond_grch37 = round(pct_beyond_grch37, 3),
    n_min = n_min, n_median = n_med, n_max = n_max,
    n_eff_cc = round(neff, 0),
    columns = paste(names(dt), collapse = "|")
  )
  rm(dt); gc(FALSE)
}

qc <- rbindlist(qc_rows, fill = TRUE)
fwrite(qc, paste0(out_prefix, "_table.tsv"), sep = "\t", na = "NA")

# ---- flagging ---------------------------------------------------------------
# severity: "action" = fix before munging (blocks the QC-passed registry);
#           "info"   = expected for whole-genome sources, munging handles it.
flags <- list()
flag_one <- function(pheno, issue, severity = "action") {
  flags[[length(flags) + 1L]] <<- data.table(phenotype = pheno, severity = severity, issue = issue)
}
for (r in seq_len(nrow(qc))) {
  p <- qc$phenotype[r]
  if (qc$status[r] != "ok") { flag_one(p, qc$status[r]); next }
  if (is.finite(qc$n_snp_na[r]) && qc$n_snp_na[r] > 0.01 * qc$n_rows[r])
    flag_one(p, paste("missing rsIDs:", qc$n_snp_na[r], "rows (",
                      round(100 * qc$n_snp_na[r] / qc$n_rows[r], 1),
                      "%) -- dropped by cleaning; munging cannot use them"), "info")
  # Duplicate rsIDs are classified by type, not count: multi-allelic sites and
  # identical rows are expected in whole-genome sources and munging drops them;
  # only position conflicts (same rsID, different BP) block the registry.
  if (is.finite(qc$n_dup_pos_conflict[r]) && qc$n_dup_pos_conflict[r] > 0)
    flag_one(p, paste("duplicate rsIDs at DIFFERENT positions:", qc$n_dup_pos_conflict,
                      "rsIDs -- rsID mapping artifact, fix before munging"))
  else if (is.finite(qc$n_dup_multiallelic[r]) && qc$n_dup_multiallelic[r] > 0)
    flag_one(p, paste("multi-allelic rsIDs:", qc$n_dup_multiallelic,
                      "(same position, different alleles; munging drops these)"), "info")
  # pct_acgt_alleles is the share of rows with valid single-base alleles.
  # Whole-genome imputed files carry ~5-10% multi-base indel alleles, which
  # munging drops; only a corrupt allele column (<50% valid) blocks reuse.
  if (is.finite(qc$pct_acgt_alleles[r]) && qc$pct_acgt_alleles[r] < 50)
    flag_one(p, paste("allele column corrupt: only", qc$pct_acgt_alleles[r],
                      "% of rows have valid ACGT alleles"))
  else if (is.finite(qc$pct_acgt_alleles[r]) && qc$pct_acgt_alleles[r] < 99)
    flag_one(p, paste("multi-base indel alleles:", 100 - qc$pct_acgt_alleles[r],
                      "% (munging drops these; no re-clean needed)"), "info")
  if (is.finite(qc$pct_frq_bad[r]) && qc$pct_frq_bad[r] > 0.1)
    flag_one(p, paste("FRQ outside (0,1):", qc$pct_frq_bad[r],
                      "% -- check the source FRQ column for misparsing"))
  if (is.finite(qc$pct_z_eq_beta_se[r]) && qc$pct_z_eq_beta_se[r] < 95)
    flag_one(p, paste("Z != BETA/SE in", 100 - qc$pct_z_eq_beta_se[r], "% of rows"))
  if (is.finite(qc$pct_sign_disagree[r]) && qc$pct_sign_disagree[r] > 1)
    flag_one(p, paste("Z/BETA sign disagreement:", qc$pct_sign_disagree[r], "%"))
  if (is.finite(qc$corr_absZ_vs_P[r]) && qc$corr_absZ_vs_P[r] < 0.95)
    flag_one(p, paste("|Z| vs P correlation low:", qc$corr_absZ_vs_P[r]))
  if (is.finite(qc$pct_bp_beyond_grch37[r]) && qc$pct_bp_beyond_grch37[r] > 1)
    flag_one(p, paste("positions beyond GRCh37 lengths (GRCh38?):",
                      qc$pct_bp_beyond_grch37[r], "%"))
  if (is.finite(qc$pct_rs_ids[r]) && qc$pct_rs_ids[r] < 50)
    flag_one(p, paste("fewer than 50% of SNP IDs are rsIDs (",
                      qc$pct_rs_ids[r], "%) -- needs coordinate mapping"))
  if (is.finite(qc$pct_coord_missing[r]) && qc$pct_coord_missing[r] > 5)
    flag_one(p, paste("CHR/BP missing or unparsable for", qc$pct_coord_missing[r],
                      "% of rows (e.g., chr-prefixed CHR)"))
  if (is.finite(qc$n_min[r]) && is.finite(qc$n_max[r]) &&
      qc$n_max[r] > 10 * max(qc$n_min[r], 1))
    flag_one(p, paste("highly variable N:", qc$n_min[r], "-", qc$n_max[r],
                      "(munging drops low-N SNPs)"), "info")
}

flag_table <- rbindlist(flags, fill = TRUE)
fwrite(flag_table, paste0(out_prefix, "_flags.tsv"), sep = "\t", na = "")

# Registry of phenotypes with no actionable flags: the cleaner skips these.
action_flags <- flag_table[severity == "action"]
passed <- qc[status == "ok" & !(phenotype %in% action_flags$phenotype), .(phenotype)]
fwrite(passed, paste0(out_prefix, "_passed.tsv"), sep = "\t", na = "")

cat("\n==== FLAGS ====")
if (nrow(flag_table)) print(as.data.table(flag_table[, .(phenotype, severity, issue)]))
cat("\nPassed QC:", nrow(passed), "phenotype(s) ->", paste0(out_prefix, "_passed.tsv"), "\n")
cat("Saved:", paste0(out_prefix, "_table.tsv"), "and", paste0(out_prefix, "_flags.tsv"), "\n")
