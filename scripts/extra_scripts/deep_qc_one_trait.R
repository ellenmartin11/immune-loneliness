# Interactive deep QC for a single cleaned sumstats file (Obesity).
# Same checks as qc_sumstats_deep.R, but for one file, printing to console.
#
# Use: set the path below, select all, run in your R session.
# (Also works headless: Rscript qc_obesity_one.R <file.tsv.gz>)

path <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/022_Obesity_clean.tsv.gz"

# ---- nothing below needs editing --------------------------------------------
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  path <- commandArgs(trailingOnly = TRUE)[[1L]]
}
suppressPackageStartupMessages(library(data.table))
if (!file.exists(path)) stop("File not found: ", path)

ACGT <- c("A", "C", "G", "T")
grch37_len <- c(249250621, 243199373, 198022430, 190214555, 181538259,
                170805979, 159345973, 145138636, 138394717, 133797422,
                135086622, 133275309, 114364328, 107043718, 102531392,
                90338345, 81195210, 78077248, 59128983, 63025520,
                48129895, 51304566)
names(grch37_len) <- as.character(1:22)

cat("Reading", path, "...\n")
dt <- fread(path, showProgress = FALSE, na.strings = c("NA", ""))
n0 <- nrow(dt)

# ---- variant ID integrity ----------------------------------------------------
n_snp_na  <- sum(is.na(dt$SNP))
n_dup_snp <- sum(duplicated(dt$SNP[!is.na(dt$SNP)]))
pct_rs    <- 100 * mean(grepl("^rs[0-9]+$", dt$SNP), na.rm = TRUE)
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

# ---- alleles -------------------------------------------------------------------
a1 <- toupper(dt$A1); a2 <- toupper(dt$A2)
pct_acgt <- 100 * mean(a1 %in% ACGT & a2 %in% ACGT, na.rm = TRUE)
n_same_allele <- sum(!is.na(a1) & !is.na(a2) & a1 == a2)

# ---- allele frequency (all-NA FRQ reported as NA, not 100% bad) ------------------
pct_frq_bad <- NA_real_; pct_frq_rare <- NA_real_
if ("FRQ" %in% names(dt)) {
  v <- dt$FRQ[!is.na(dt$FRQ)]
  if (length(v)) {
    pct_frq_bad  <- 100 * mean(!is.finite(v) | v <= 0 | v >= 1)
    pct_frq_rare <- 100 * mean(v < 0.01 | v > 0.99)
  }
}

# ---- signed stats ----------------------------------------------------------------
z_consistent <- NA_real_; sign_disagree <- NA_real_
if (all(c("Z", "BETA", "SE") %in% names(dt))) {
  ok <- is.finite(dt$Z) & is.finite(dt$BETA) & is.finite(dt$SE) & dt$SE > 0
  if (any(ok)) {
    z_ref <- dt$BETA[ok] / dt$SE[ok]
    z_consistent <- 100 * mean(abs(dt$Z[ok] - z_ref) / pmax(abs(z_ref), 1) < 0.05)
    sign_disagree <- 100 * mean(sign(dt$Z[ok]) != sign(dt$BETA[ok]))
  } else {
    z_consistent <- NaN; sign_disagree <- NaN
  }
}
zp_corr <- NA_real_
if (all(c("Z", "P") %in% names(dt))) {
  ok <- is.finite(dt$Z) & is.finite(dt$P) & dt$P > 0 & dt$P <= 1
  if (sum(ok) > 100) zp_corr <- cor(abs(dt$Z[ok]), qnorm(dt$P[ok] / 2, lower.tail = FALSE))
}
pct_p_bad <- NA_real_
if ("P" %in% names(dt)) {
  pct_p_bad <- 100 * mean(!is.finite(dt$P) | dt$P <= 0 | dt$P > 1, na.rm = TRUE)
}

# ---- coordinates / build -----------------------------------------------------------
chr <- as.character(dt$CHR)
chr[chr %in% c("23", "X")] <- "X"; chr[chr %in% c("24", "Y")] <- "Y"
chr[chr %in% c("25", "MT", "M")] <- "MT"
autosomes <- chr %in% as.character(1:22)
pct_coord_missing <- 100 * mean(!is.finite(dt$BP) | is.na(dt$BP) | !autosomes)
pct_beyond <- NA_real_
if (sum(autosomes) > 0) {
  beyond <- dt$BP[autosomes] > grch37_len[chr[autosomes]]
  pct_beyond <- 100 * mean(beyond, na.rm = TRUE)
}

# ---- sample size ---------------------------------------------------------------------
n_min <- n_med <- n_max <- neff <- NA_real_
if ("N" %in% names(dt) && any(is.finite(dt$N))) {
  n_min <- min(dt$N, na.rm = TRUE); n_med <- median(dt$N, na.rm = TRUE)
  n_max <- max(dt$N, na.rm = TRUE)
}
if (all(c("N_CASE", "N_CONTROL") %in% names(dt)) &&
    any(is.finite(dt$N_CASE)) && any(is.finite(dt$N_CONTROL))) {
  neff <- 4 / (1 / median(dt$N_CASE, na.rm = TRUE) + 1 / median(dt$N_CONTROL, na.rm = TRUE))
}

# ---- report ----------------------------------------------------------------------------
cat("\n==== Deep QC:", basename(path), "====\n")
cat(sprintf("rows: %s | unique SNP: %s | columns: %s\n",
            format(n0, big.mark = ","), format(uniqueN(dt$SNP), big.mark = ","),
            paste(names(dt), collapse = "|")))

fmt <- function(x, digits = 2) if (is.na(x)) "NA" else formatC(x, format = "f", digits = digits)
cat(sprintf("  rsIDs:                %s%%  (NA: %s, dups: %s = %s multi-allelic / %s pos-conflict / %s identical)\n",
            fmt(pct_rs), format(n_snp_na, big.mark = ","), format(n_dup_snp, big.mark = ","),
            n_dup_multiallelic, n_dup_pos_conflict, n_dup_identical))
cat(sprintf("  ACGT alleles:         %s%%  (same-allele rows: %s)\n", fmt(pct_acgt), n_same_allele))
cat(sprintf("  FRQ bad / rare:       %s%% / %s%%\n", fmt(pct_frq_bad, 3), fmt(pct_frq_rare)))
cat(sprintf("  Z == BETA/SE:         %s%%  (sign disagreement: %s%%)\n",
            fmt(z_consistent), fmt(sign_disagree, 3)))
cat(sprintf("  corr(|Z|, P):         %s\n", fmt(zp_corr, 4)))
cat(sprintf("  P bad:                %s%%\n", fmt(pct_p_bad, 3)))
cat(sprintf("  coord missing:        %s%%  (beyond GRCh37: %s%%)\n",
            fmt(pct_coord_missing), fmt(pct_beyond, 3)))
cat(sprintf("  N min/median/max:     %s / %s / %s  (Neff from Ncase/Ncon: %s)\n",
            fmt(n_min, 0), fmt(n_med, 0), fmt(n_max, 0), fmt(neff, 0)))

# ---- flags (same thresholds as qc_sumstats_deep.R) ---------------------------------------
flags <- list()
add <- function(sev, issue) flags[[length(flags) + 1L]] <<- data.table(severity = sev, issue = issue)
if (is.nan(z_consistent) || is.nan(sign_disagree))
  add("action", "Z/BETA/SE present but ZERO rows with finite values - BETA or SE column likely empty/mis-mapped")
if (n_dup_pos_conflict > 0)
  add("action", paste("duplicate rsIDs at DIFFERENT positions:", n_dup_pos_conflict))
if (pct_acgt < 50) { add("action", paste("allele column corrupt: only", pct_acgt, "% valid ACGT"))
} else if (pct_acgt < 99) { add("info", paste("multi-base indel alleles:", round(100 - pct_acgt, 2), "% (munging drops)")) }
if (is.finite(pct_frq_bad) && pct_frq_bad > 0.1) add("action", paste("FRQ outside (0,1):", pct_frq_bad, "%"))
if (is.finite(z_consistent) && z_consistent < 95) add("action", paste("Z != BETA/SE in", 100 - z_consistent, "% of rows"))
if (is.finite(sign_disagree) && sign_disagree > 1) add("action", paste("Z/BETA sign disagreement:", sign_disagree, "%"))
if (is.finite(zp_corr) && zp_corr < 0.95) add("action", paste("|Z| vs P correlation low:", zp_corr))
if (is.finite(pct_beyond) && pct_beyond > 1) add("action", paste("positions beyond GRCh37 lengths (GRCh38?):", pct_beyond, "%"))
if (pct_rs < 50) add("action", paste("fewer than 50% rsIDs:", pct_rs, "% -- needs coordinate mapping"))
if (pct_coord_missing > 5) add("action", paste("CHR/BP missing/unparsable:", pct_coord_missing, "% of rows"))
if (n_snp_na > 0.01 * n0) add("info", paste("missing rsIDs:", n_snp_na, "rows (dropped by cleaning)"))
if (n_dup_multiallelic > 0) add("info", paste("multi-allelic rsIDs:", n_dup_multiallelic, "(munging drops)"))
if (is.finite(n_min) && is.finite(n_max) && n_max > 10 * max(n_min, 1))
  add("info", paste("highly variable N:", n_min, "-", n_max, "(munging drops low-N SNPs)"))

cat("\n==== FLAGS ====\n")
if (length(flags)) {
  ft <- rbindlist(flags, fill = TRUE)
  print(ft)
  cat("\nVerdict:", if (any(ft$severity == "action")) "ACTION NEEDED - do not munge yet"
      else "PASS (info-level flags only)", "\n")
} else cat("None - fully clean.\n\nVerdict: PASS\n")
rm(dt); gc(FALSE)

