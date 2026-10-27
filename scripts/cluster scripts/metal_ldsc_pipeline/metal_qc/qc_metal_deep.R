#!/usr/bin/env Rscript

# Deep QC for METAL IVW-combined summary statistics (e.g., PanUKB + FinnGen).
#
# Differences from qc_sumstats_deep.R:
#   - Reads the METAL manifest (phenotype, members, N, build, path, variants)
#     instead of clean_manifest.tsv; N comes from the manifest (METAL output
#     has no per-SNP N).
#   - Build-aware coordinate checks (GRCh37 or GRCh38 lengths).
#   - Alias-matches METAL-style headers (MarkerName, Allele1/2, Effect,
#     StdErr, P-value, Freq1, ...) as well as standardized ones.
#   - Z is computed as BETA/SE when absent, then checked against P.
#
# Usage: Rscript qc_metal_deep.R <metal_manifest.tsv> <out_prefix>
# Output: <out_prefix>_table.tsv and <out_prefix>_flags.tsv

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: qc_metal_deep.R <metal_manifest.tsv> <out_prefix>")
manifest_path <- args[[1L]]
out_prefix    <- args[[2L]]

manifest <- fread(manifest_path)
stopifnot(all(c("phenotype", "path", "build", "N") %in% names(manifest)))

ACGT <- c("A", "C", "G", "T")

grch37_len <- c(249250621, 243199373, 198022430, 190214555, 181538259,
                170805979, 159345973, 145138636, 138394717, 133797422,
                135086622, 133275309, 114364328, 107043718, 102531392,
                90338345, 81195210, 78077248, 59128983, 63025520,
                48129895, 51304566)
grch38_len <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                170805979, 159345973, 145138636, 138394717, 133797422,
                135086622, 133275309, 114364328, 107043718, 101991189,
                90338345, 83257441, 80373285, 59419756, 64444167,
                46709983, 50818468)
names(grch37_len) <- names(grch38_len) <- as.character(1:22)

norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Alias maps: standardized names first, METAL output names second.
alias <- list(
  SNP  = c("SNP", "rsid", "rs", "markername", "marker", "variant_id"),
  CHR  = c("CHR", "chromosome", "chrom", "chr"),
  BP   = c("BP", "pos", "position", "bp_hg38", "bp_hg37", "base_pair_location"),
  A1   = c("A1", "allele1", "effect_allele", "tested_allele", "alt", "inc_allele"),
  A2   = c("A2", "allele2", "other_allele", "ref"),
  BETA = c("BETA", "effect", "log_odds", "beta_1", "beta_meta"),
  SE   = c("SE", "stderr", "std_err", "se_1", "se_meta", "stderr_meta"),
  P    = c("P", "pvalue", "p_value", "p", "p_bhat", "pvalue_meta"),
  Z    = c("Z", "zscore", "z_score"),
  FRQ  = c("FRQ", "freq1", "freq", "eaf", "af_alt", "freq_meta",
           "freq_cohort1", "af", "effect_af"),
  INFO = c("INFO", "info_score", "impinfo"),
  DIRECTION = c("Direction", "direction"),
  HET_ISQ   = c("HetISq", "hetisq", "isq"),
  HET_P     = c("HetPVal", "hetpval", "het_p")
)

find_col <- function(headers, candidates) {
  hit <- match(norm(candidates), norm(headers))
  if (all(is.na(hit))) NA_character_ else headers[hit[!is.na(hit)][[1L]]]
}

qc_rows <- list()
for (i in seq_len(nrow(manifest))) {
  pheno <- manifest$phenotype[i]
  path  <- manifest$path[i]
  build <- toupper(manifest$build[i])
  N_man <- suppressWarnings(as.numeric(gsub(",", "", as.character(manifest$N[i]))))
  cat("QC:", pheno, "\n")

  if (!file.exists(path)) {
    qc_rows[[length(qc_rows) + 1L]] <- data.table(phenotype = pheno, status = "FILE NOT FOUND")
    next
  }

  dt <- fread(path, showProgress = FALSE, na.strings = c("NA", ""))
  n0 <- nrow(dt)

  colmap <- vapply(alias, function(cands) find_col(names(dt), cands), character(1))
  getc <- function(key) if (is.na(colmap[[key]])) NULL else dt[[colmap[[key]]]]

  # ---- variant ID: rsIDs vs coordinates -------------------------------------
  snp <- getc("SNP")
  a1 <- toupper(getc("A1")); a2 <- toupper(getc("A2"))

  # If CHR/BP columns are absent but IDs are chr:pos (or chr:pos_ref_alt),
  # parse coordinates out of the ID so build checks still run.
  if (is.null(getc("CHR")) && is.null(getc("BP")) && !is.null(snp)) {
    parsed <- regmatches(as.character(snp),
                         regexec("^(?:chr)?([0-9XYM]+)[:_]([0-9]+)(?:[:_][A-Za-z]+[:_][A-Za-z]+)?$",
                                 snp, perl = TRUE))
    ok_parse <- lengths(parsed) > 0
    if (any(ok_parse)) {
      chr_parsed <- vapply(parsed[ok_parse], function(m) m[[2]], character(1))
      bp_parsed <- as.numeric(vapply(parsed[ok_parse], function(m) m[[3]], character(1)))
      dt[, CHR := NA_character_]; dt[ok_parse, CHR := chr_parsed]
      dt[, BP := NA_real_];       dt[ok_parse, BP := bp_parsed]
      colmap[["CHR"]] <- "CHR"; colmap[["BP"]] <- "BP"
    }
  }

  pct_rs   <- if (is.null(snp)) NA_real_ else 100 * mean(grepl("^rs[0-9]+$", snp), na.rm = TRUE)
  pct_pseudo <- if (is.null(snp)) NA_real_ else
    100 * mean(grepl("^(?:chr)?[0-9XYM]+[:_][0-9]+[:_][A-Za-z]+[:_][A-Za-z]+$", snp), na.rm = TRUE)
  n_snp_na <- if (is.null(snp)) n0 else sum(is.na(snp))
  n_dup_snp <- if (is.null(snp)) 0L else sum(duplicated(snp[!is.na(snp)]))

  # duplicate classification (multi-allelic vs position conflict vs identical)
  # NOTE: build a plain data.table from vectors -- calling getc() inside a
  # dt[i, j] subset would return full-length columns, and the ID column is not
  # necessarily named "SNP" (METAL files use MarkerName).
  n_dup_multiallelic <- 0L; n_dup_pos_conflict <- 0L; n_dup_identical <- 0L
  dup_classifiable <- FALSE
  if (n_dup_snp > 0 && !is.null(snp) && !is.null(getc("BP")) && !is.null(a1) && !is.null(a2)) {
    dup_classifiable <- TRUE
    dup_ids <- snp[duplicated(snp) & !is.na(snp)]
    bp_dup <- suppressWarnings(as.numeric(getc("BP")))
    idx <- snp %in% dup_ids
    d <- data.table(id = snp[idx], bp = bp_dup[idx], al = paste(a1[idx], a2[idx]))
    d <- d[, .(n_pos = uniqueN(bp), n_al = uniqueN(al)), by = id]
    n_dup_pos_conflict <- d[n_pos > 1, .N]
    n_dup_multiallelic <- d[n_pos == 1 & n_al > 1, .N]
    n_dup_identical    <- d[n_pos == 1 & n_al == 1, .N]
    rm(d, bp_dup, idx)
  }

  # ---- alleles ---------------------------------------------------------------
  pct_acgt <- if (is.null(a1) || is.null(a2)) NA_real_ else
    100 * mean(a1 %in% ACGT & a2 %in% ACGT, na.rm = TRUE)
  n_same_allele <- if (is.null(a1) || is.null(a2)) NA_integer_ else
    sum(!is.na(a1) & !is.na(a2) & a1 == a2)

  # ---- frequencies -----------------------------------------------------------
  frq <- getc("FRQ")
  frq <- if (!is.null(frq)) suppressWarnings(as.numeric(gsub(",", "", as.character(frq)))) else NULL
  pct_frq_bad  <- if (is.null(frq)) NA_real_ else 100 * mean(!is.finite(frq) | frq <= 0 | frq >= 1, na.rm = TRUE)
  pct_frq_rare <- if (is.null(frq)) NA_real_ else 100 * mean(frq < 0.01 | frq > 0.99, na.rm = TRUE)

  # ---- effect sizes and Z ------------------------------------------------------
  # as_num_col keeps NULL (column absent) distinct from numeric(0).
  as_num_col <- function(x) if (is.null(x)) NULL else
    suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
  beta <- as_num_col(getc("BETA"))
  se   <- as_num_col(getc("SE"))
  z    <- as_num_col(getc("Z"))
  if (is.null(z) && !is.null(beta) && !is.null(se)) z <- beta / se
  p <- as_num_col(getc("P"))

  z_consistent <- NA_real_; sign_disagree <- NA_real_
  if (!is.null(beta) && !is.null(se) && !is.null(z)) {
    ok <- is.finite(z) & is.finite(beta) & is.finite(se) & se > 0
    if (sum(ok) > 0) {
      z_ref <- beta[ok] / se[ok]
      z_consistent <- 100 * mean(abs(z[ok] - z_ref) / pmax(abs(z_ref), 1) < 0.05, na.rm = TRUE)
      sign_disagree <- 100 * mean(sign(z[ok]) != sign(beta[ok]), na.rm = TRUE)
    }
  }
  zp_corr <- NA_real_
  if (!is.null(z) && !is.null(p)) {
    ok <- is.finite(z) & is.finite(p) & p > 0 & p <= 1
    if (sum(ok) > 100) zp_corr <- cor(abs(z[ok]), qnorm(p[ok] / 2, lower.tail = FALSE))
  }
  pct_p_bad <- if (is.null(p)) NA_real_ else 100 * mean(!is.finite(p) | p <= 0 | p > 1, na.rm = TRUE)

  # ---- coordinates, build-aware ------------------------------------------------
  chr <- toupper(as.character(getc("CHR")))
  chr[chr %in% c("23", "X")] <- "X"; chr[chr %in% c("24", "Y")] <- "Y"
  chr[chr %in% c("25", "MT", "M")] <- "MT"
  bp <- suppressWarnings(as.numeric(getc("BP")))
  # build string may be "GRCh38"/"GRCh37"/"hg38"/... -- match on the number
  ref_len <- if (grepl("38", build)) grch38_len else grch37_len
  autosomes <- chr %in% names(ref_len)
  pct_coord_missing <- 100 * mean(!is.finite(bp) | is.na(bp) | !autosomes, na.rm = TRUE)
  pct_beyond <- NA_real_
  if (sum(autosomes) > 0) {
    beyond <- bp[autosomes] > ref_len[chr[autosomes]]
    pct_beyond <- 100 * mean(beyond, na.rm = TRUE)
  }

  # ---- meta-analysis specifics (optional columns) -------------------------------
  direction <- getc("DIRECTION")
  pct_dir_unclear <- if (is.null(direction)) NA_real_ else
    100 * mean(grepl("?", direction, fixed = TRUE), na.rm = TRUE)
  het_isq <- as_num_col(getc("HET_ISQ"))
  med_het_isq <- if (is.null(het_isq)) NA_real_ else as.numeric(median(het_isq, na.rm = TRUE))

  qc_rows[[length(qc_rows) + 1L]] <- data.table(
    phenotype = pheno, status = "ok", build = build,
    members = manifest$members[i], N_manifest = N_man,
    n_rows = n0, n_unique_snp = uniqueN(snp),
    n_snp_na = n_snp_na, n_dup_snp = n_dup_snp,
    n_dup_multiallelic = n_dup_multiallelic,
    n_dup_pos_conflict = n_dup_pos_conflict,
    n_dup_identical = n_dup_identical,
    pct_rs_ids = round(pct_rs, 2), pct_pseudo_ids = round(pct_pseudo, 2),
    pct_acgt_alleles = round(pct_acgt, 2), n_same_allele = n_same_allele,
    pct_frq_bad = round(pct_frq_bad, 3), pct_frq_rare = round(pct_frq_rare, 2),
    pct_z_eq_beta_se = round(z_consistent, 2),
    pct_sign_disagree = round(sign_disagree, 3),
    corr_absZ_vs_P = round(zp_corr, 4), pct_p_bad = round(pct_p_bad, 3),
    pct_coord_missing = round(pct_coord_missing, 2),
    pct_bp_beyond_build = round(pct_beyond, 3),
    pct_direction_unclear = round(pct_dir_unclear, 2),
    median_het_isq = med_het_isq,
    columns = paste(names(dt), collapse = "|")
  )
  rm(dt); gc(FALSE)
}

qc <- rbindlist(qc_rows, fill = TRUE)
fwrite(qc, paste0(out_prefix, "_table.tsv"), sep = "\t", na = "NA")

# ---- flagging (same severity scheme as the primary deep QC) --------------------
flags <- list()
flag_one <- function(pheno, issue, severity = "action") {
  flags[[length(flags) + 1L]] <<- data.table(phenotype = pheno, severity = severity, issue = issue)
}
for (r in seq_len(nrow(qc))) {
  p <- qc$phenotype[r]
  if (qc$status[r] != "ok") { flag_one(p, qc$status[r]); next }
  if (is.finite(qc$n_dup_pos_conflict[r]) && qc$n_dup_pos_conflict[r] > 0)
    flag_one(p, paste("duplicate IDs at DIFFERENT positions:", qc$n_dup_pos_conflict[r],
                      "-- fix before munging"))
  else if (is.finite(qc$n_dup_multiallelic[r]) && qc$n_dup_multiallelic[r] > 0)
    flag_one(p, paste("multi-allelic IDs:", qc$n_dup_multiallelic[r],
                      "(same position, different alleles; munging drops these)"), "info")
  else if (is.finite(qc$n_dup_snp[r]) && qc$n_dup_snp[r] > 0 &&
             is.finite(qc$n_dup_multiallelic[r]) && is.finite(qc$n_dup_pos_conflict[r]) &&
             qc$n_dup_multiallelic[r] + qc$n_dup_pos_conflict[r] + qc$n_dup_identical[r] == 0)
    flag_one(p, paste("duplicate IDs present (", qc$n_dup_snp[r],
                      ") but could not classify -- missing BP or allele columns"))
  if (is.finite(qc$n_same_allele[r]) && qc$n_same_allele[r] > 0)
    flag_one(p, paste("same-allele rows (A1 == A2):", qc$n_same_allele[r]), "info")
  if (is.finite(qc$pct_acgt_alleles[r]) && qc$pct_acgt_alleles[r] < 50)
    flag_one(p, paste("allele column corrupt: only", qc$pct_acgt_alleles[r], "% valid ACGT"))
  else if (is.finite(qc$pct_acgt_alleles[r]) && qc$pct_acgt_alleles[r] < 99)
    flag_one(p, paste("multi-base indel alleles:", 100 - qc$pct_acgt_alleles[r],
                      "% (munging drops these)"), "info")
  if (is.finite(qc$pct_frq_bad[r]) && qc$pct_frq_bad[r] > 0.1)
    flag_one(p, paste("FRQ outside (0,1):", qc$pct_frq_bad[r], "% -- check the Freq column mapping"))
  if (is.finite(qc$pct_z_eq_beta_se[r]) && qc$pct_z_eq_beta_se[r] < 95)
    flag_one(p, paste("Z != BETA/SE in", 100 - qc$pct_z_eq_beta_se[r], "% of rows"))
  if (is.finite(qc$pct_sign_disagree[r]) && qc$pct_sign_disagree[r] > 1)
    flag_one(p, paste("Z/BETA sign disagreement:", qc$pct_sign_disagree[r], "%"))
  if (is.finite(qc$corr_absZ_vs_P[r]) && qc$corr_absZ_vs_P[r] < 0.95)
    flag_one(p, paste("|Z| vs P correlation low:", qc$corr_absZ_vs_P[r]))
  if (is.finite(qc$pct_bp_beyond_build[r]) && qc$pct_bp_beyond_build[r] > 1)
    flag_one(p, paste("positions beyond", qc$build[r], "lengths:", qc$pct_bp_beyond_build[r],
                      "% -- build mismatch between members?"))
  if (is.finite(qc$pct_coord_missing[r]) && qc$pct_coord_missing[r] > 5)
    flag_one(p, paste("CHR/BP missing or unparsable for", qc$pct_coord_missing[r], "% of rows"))
  if (is.finite(qc$pct_rs_ids[r]) && qc$pct_rs_ids[r] < 50)
    flag_one(p, paste("fewer than 50% rsIDs (", qc$pct_rs_ids[r],
                      "%) -- coordinate mapping needed before munging"))
  if (is.finite(qc$pct_pseudo_ids[r]) && qc$pct_pseudo_ids[r] > 50)
    flag_one(p, paste("SNP IDs are chr:pos_allele pseudo-IDs (", qc$pct_pseudo_ids[r],
                      "%) -- rsID mapping needed before munging"))
  if (is.finite(qc$pct_direction_unclear[r]) && qc$pct_direction_unclear[r] > 50)
    flag_one(p, paste("Direction unclear in", qc$pct_direction_unclear[r],
                      "% of rows -- check member cohort alignment"), "info")
}

if (length(flags) == 0L) {
  flag_table <- data.table(phenotype = character(), severity = character(), issue = character())
} else {
  flag_table <- rbindlist(flags, fill = TRUE)
}
fwrite(flag_table, paste0(out_prefix, "_flags.tsv"), sep = "\t", na = "")

cat("\n==== FLAGS ====")
if (nrow(flag_table)) print(flag_table[, .(phenotype, severity, issue)])
cat("\nSaved:", paste0(out_prefix, "_table.tsv"), "and", paste0(out_prefix, "_flags.tsv"), "\n")
