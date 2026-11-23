#!/usr/bin/env Rscript

# Run bidirectional TwoSampleMR analyses for T2D, Depression, Triglycerides,
# and each of the AoU, PanUKB, and METAL Loneliness datasets.
# Thresholds follow the current notebook: relaxed 5e-6 and stringent 5e-8.

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: run_bidirectional_mr.R <output directory> <scratch cache directory>")
}
out_dir <- args[[1L]]
cache_root <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(cache_root, "pair_cache_v2")
instrument_cache_dir <- file.path(cache_root, "instrument_cache")
normalized_cache_dir <- file.path(cache_root, "normalized_sumstats")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(instrument_cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(normalized_cache_dir, recursive = TRUE, showWarnings = FALSE)

THRESHOLDS <- c(relaxed_5e6 = 5e-6, stringent_5e8 = 5e-8)
CLUMP_KB <- 10000
CLUMP_R2 <- 0.001
SEED <- 100L

PLINK_BIN <- "/vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink"
LD_BFILE <- "/gpfs/gibbs/pi/polimanti/LDSC_1000G/g1000_eur"

# METAL Loneliness and AoU coordinates are GRCh38. PanUKB positions are
# GRCh37, but its rsIDs allow build-independent matching. Other traits are
# matched by rsID; their coordinates are not used for coordinate matching.
TRAIT_PATHS <- c(
  "AoU Loneliness" = "/gpfs/gibbs/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/UCLA_loneliness_scale_eur.tsv",
  "PanUKB Loneliness" = "/gpfs/gibbs/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/005_Loneliness_isolation_PanUKB_clean.tsv.gz",
  "METAL Loneliness" = "/gpfs/gibbs/scratch/polimanti/ejm96/metal_ldsc_pipeline/loneliness_metal/Loneliness_PanUKB_AoU_METAL_combined.tsv.gz",
  "Triglycerides" = "/gpfs/gibbs/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw/Triglycerides_METAL_combined.tsv.gz",
  "Depression" = "/gpfs/gibbs/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/012_Depression_PGC_clean.tsv.gz",
  "Type 2 Diabetes" = "/gpfs/gibbs/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/006_Type_2_Diabetes_clean.tsv.gz"
)

TRAIT_COORD_BUILD <- c(
  "AoU Loneliness" = "GRCh38",
  "PanUKB Loneliness" = "GRCh37",
  "METAL Loneliness" = "GRCh38",
  "Triglycerides" = NA_character_,
  "Depression" = NA_character_,
  "Type 2 Diabetes" = NA_character_
)

CORE_TRAITS <- c("Triglycerides", "Depression", "Type 2 Diabetes")
LONELINESS_TRAITS <- c("AoU Loneliness", "PanUKB Loneliness", "METAL Loneliness")

missing_paths <- TRAIT_PATHS[!file.exists(TRAIT_PATHS)]
if (length(missing_paths)) {
  stop("Input file(s) not found: ", paste(missing_paths, collapse = "; "))
}
if (!file.exists(PLINK_BIN)) stop("PLINK executable not found: ", PLINK_BIN)
if (!all(file.exists(paste0(LD_BFILE, c(".bed", ".bim", ".fam"))))) {
  stop("PLINK reference panel is missing .bed/.bim/.fam files: ", LD_BFILE)
}

needed_bioc <- c("SNPlocs.Hsapiens.dbSNP155.GRCh38", "BSgenome", "GenomicRanges",
                 "IRanges", "S4Vectors", "Biostrings", "GenomeInfoDb")
missing_bioc <- needed_bioc[
  !vapply(needed_bioc, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc)) {
  stop("Install required GRCh38 mapping package(s) in the active R library: ",
       paste(missing_bioc, collapse = ", "))
}
suppressPackageStartupMessages(
  library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
)

message("TwoSampleMR version: ", as.character(packageVersion("TwoSampleMR")))
message("Thresholds: relaxed 5e-6; stringent 5e-8")
message("Coordinate matching is restricted to GRCh38.")

norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
find_col <- function(nms, choices) {
  i <- match(norm(choices), norm(nms), nomatch = 0L)
  i <- i[i > 0L]
  if (length(i)) nms[[i[[1L]]]] else NULL
}
as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
safe_file <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

parse_coordinate_id <- function(ids) {
  matches <- regmatches(
    as.character(ids),
    regexec("^(?:chr)?([^:]+):([0-9]+)(?::.*)?$", as.character(ids),
            perl = TRUE, ignore.case = TRUE)
  )
  ok <- lengths(matches) >= 3L
  chr <- rep(NA_character_, length(ids))
  bp <- rep(NA_real_, length(ids))
  if (any(ok)) {
    chr[ok] <- vapply(matches[ok], `[[`, character(1), 2L)
    bp[ok] <- suppressWarnings(as.numeric(vapply(matches[ok], `[[`, character(1), 3L)))
  }
  list(chr = chr, bp = bp, parsed = ok)
}

normalize_sumstats <- function(x, label, coordinate_build = NA_character_) {
  req <- c("SNP", "A1", "A2", "BETA", "SE", "P")
  missing <- setdiff(req, names(x))
  if (length(missing)) {
    stop(label, " missing required fields: ", paste(missing, collapse = ", "))
  }
  if (!"CHR" %in% names(x)) x[, CHR := NA_character_]
  if (!"BP" %in% names(x)) x[, BP := NA_real_]
  if (!"N" %in% names(x)) x[, N := NA_real_]
  if (!"EAF" %in% names(x)) x[, EAF := NA_real_]

  x[, `:=`(
    SNP = trimws(as.character(SNP)),
    CHR = as.character(CHR), BP = as_num(BP),
    A1 = toupper(trimws(as.character(A1))),
    A2 = toupper(trimws(as.character(A2))),
    BETA = as_num(BETA), SE = as_num(SE), P = as_num(P),
    N = as_num(N), EAF = as_num(EAF)
  )]

  parsed <- parse_coordinate_id(x$SNP)
  needs_coordinates <- (is.na(x$CHR) | !nzchar(x$CHR) | is.na(x$BP)) & parsed$parsed
  if (any(needs_coordinates)) {
    x$CHR[needs_coordinates] <- parsed$chr[needs_coordinates]
    x$BP[needs_coordinates] <- parsed$bp[needs_coordinates]
  }

  invalid <- is.na(x$SNP) | !nzchar(x$SNP) |
    !x$A1 %chin% c("A", "C", "G", "T") |
    !x$A2 %chin% c("A", "C", "G", "T") |
    !is.finite(x$BETA) | !is.finite(x$SE) | x$SE <= 0 |
    !is.finite(x$P) | x$P <= 0 | x$P > 1
  x <- x[!invalid]

  duplicate_ids <- x[, .N, by = SNP][N > 1L, SNP]
  if (length(duplicate_ids)) {
    message(label, ": excluding all rows for ", length(duplicate_ids),
            " duplicated SNP identifiers.")
    x <- x[!SNP %chin% duplicate_ids]
  }
  if (!nrow(x)) stop("No usable variants remain for ", label)
  attr(x, "coordinate_build") <- coordinate_build
  x
}

read_aou <- function(path, label) {
  header <- names(fread(path, nrows = 0L, showProgress = FALSE))
  needed <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "N", "TEST", "BETA", "SE", "LOG10P")
  missing <- setdiff(needed, header)
  if (length(missing)) stop(label, " missing AoU fields: ", paste(missing, collapse = ", "))
  x <- fread(path, select = needed, showProgress = TRUE)
  x <- x[TEST == "ADD", .(
    SNP = as.character(ID), CHR = as.character(CHROM), BP = as_num(GENPOS),
    A1 = toupper(as.character(ALLELE1)), A2 = toupper(as.character(ALLELE0)),
    BETA = as_num(BETA), SE = as_num(SE),
    P = pmax(10^(-as_num(LOG10P)), .Machine$double.xmin),
    N = as_num(N), EAF = NA_real_
  )]
  normalize_sumstats(x, label, coordinate_build = "GRCh38")
}

read_standard <- function(path, label, coordinate_build = NA_character_) {
  header <- names(fread(path, nrows = 0L, showProgress = FALSE, check.names = FALSE))
  pick <- function(choices) find_col(header, choices)
  cols <- unique(na.omit(c(
    pick(c("SNP", "rsid", "rsids", "MarkerName", "variant_id", "ID")),
    pick(c("CHR", "CHROM", "chromosome", "chr", "#chrom")),
    pick(c("BP", "POS", "GENPOS", "base_pair_location", "position")),
    pick(c("A1", "effect_allele", "Allele1", "ALT")),
    pick(c("A2", "other_allele", "Allele2", "REF")),
    pick(c("BETA", "beta", "Effect")), pick(c("SE", "standard_error", "sebeta", "StdErr")),
    pick(c("P", "p_value", "pval", "P.value", "Pvalue")),
    pick(c("LOG10P", "NEG_LOG10_P", "mlogp")),
    pick(c("N", "samplesize", "sample_size")),
    pick(c("EAF", "AF", "effect_allele_frequency", "FRQ", "af_alt"))
  )))
  if (!length(cols)) stop("Could not identify columns in ", label)
  x <- fread(path, select = cols, showProgress = TRUE, check.names = FALSE)
  get <- function(choices, default = NA) {
    nm <- find_col(names(x), choices)
    if (is.null(nm)) rep(default, nrow(x)) else x[[nm]]
  }
  p <- get(c("P", "p_value", "pval", "P.value", "Pvalue"))
  p_col <- find_col(names(x), c("P", "p_value", "pval", "P.value", "Pvalue"))
  if (is.null(p_col)) {
    lp_col <- find_col(names(x), c("LOG10P", "NEG_LOG10_P", "mlogp"))
    if (is.null(lp_col)) stop(label, " needs a P or -log10(P) column.")
    lp <- get(c("LOG10P", "NEG_LOG10_P", "mlogp"))
    p <- pmax(10^(-as_num(lp)), .Machine$double.xmin)
  }
  out <- data.table(
    SNP = get(c("SNP", "rsid", "rsids", "MarkerName", "variant_id", "ID")),
    CHR = get(c("CHR", "CHROM", "chromosome", "chr", "#chrom")),
    BP = get(c("BP", "POS", "GENPOS", "base_pair_location", "position")),
    A1 = get(c("A1", "effect_allele", "Allele1", "ALT")),
    A2 = get(c("A2", "other_allele", "Allele2", "REF")),
    BETA = get(c("BETA", "beta", "Effect")),
    SE = get(c("SE", "standard_error", "sebeta", "StdErr")),
    P = p,
    N = get(c("N", "samplesize", "sample_size")),
    EAF = get(c("EAF", "AF", "effect_allele_frequency", "FRQ", "af_alt"))
  )
  normalize_sumstats(out, label, coordinate_build)
}

trait_cache_path <- function(trait) {
  info <- file.info(TRAIT_PATHS[[trait]])
  fingerprint <- paste(info$size, as.numeric(info$mtime), sep = "_")
  file.path(normalized_cache_dir,
            paste0(safe_file(trait), "__", fingerprint, ".rds"))
}

load_trait <- function(trait) readRDS(trait_cache_path(trait))

# Normalize each large source file once, one at a time, and keep the canonical
# copy in scratch. The analysis then loads only the exposure/outcome pair being
# processed, avoiding holding all six full GWAS files in RAM at once.
qc_rows <- vector("list", length(TRAIT_PATHS))
names(qc_rows) <- names(TRAIT_PATHS)
for (nm in names(TRAIT_PATHS)) {
  cache_path <- trait_cache_path(nm)
  if (file.exists(cache_path)) {
    message("Checking normalized cache for ", nm, ": ", cache_path)
    x <- tryCatch(load_trait(nm), error = function(e) NULL)
    if (is.null(x)) {
      warning("Removing incomplete/unreadable normalized cache for ", nm,
              "; it will be rebuilt: ", cache_path)
      unlink(cache_path)
    }
  }
  if (!file.exists(cache_path)) {
    message("Reading and normalizing ", nm)
    x <- if (identical(nm, "AoU Loneliness")) {
      read_aou(TRAIT_PATHS[[nm]], nm)
    } else {
      read_standard(TRAIT_PATHS[[nm]], nm, TRAIT_COORD_BUILD[[nm]])
    }
    temp_cache <- tempfile(
      pattern = paste0(basename(cache_path), ".tmp_"),
      tmpdir = normalized_cache_dir
    )
    message("Writing compressed normalized cache: ", cache_path)
    tryCatch(
      saveRDS(x, temp_cache, compress = "gzip"),
      error = function(e) {
        unlink(temp_cache)
        stop("Could not write normalized cache for ", nm, " at ", cache_path,
             ". Check scratch quota/free space and directory permissions. Original error: ",
             conditionMessage(e))
      }
    )
    if (!file.rename(temp_cache, cache_path)) {
      unlink(temp_cache)
      stop("Could not move completed normalized cache into place: ", cache_path,
           ". Check scratch quota/free space and directory permissions.")
    }
  }
  if (!exists("x", inherits = FALSE)) x <- load_trait(nm)
  qc_rows[[nm]] <- data.table(
    trait = nm, path = TRAIT_PATHS[[nm]], coordinate_build = TRAIT_COORD_BUILD[[nm]],
    n_rows = nrow(x), n_rsid = sum(grepl("^rs[0-9]+$", x$SNP, ignore.case = TRUE)),
    n_coordinate_id = sum(!grepl("^rs[0-9]+$", x$SNP, ignore.case = TRUE)),
    n_rows_with_chr_bp = sum(!is.na(x$CHR) & is.finite(x$BP))
  )
  rm(x)
  gc(verbose = FALSE)
}
qc_rows <- rbindlist(qc_rows, fill = TRUE)
fwrite(qc_rows, file.path(out_dir, "input_build_id_qc.csv"))
print(qc_rows)

dbsnp <- SNPlocs.Hsapiens.dbSNP155.GRCh38
available_seqnames <- GenomeInfoDb::seqlevels(dbsnp)
iupac <- Biostrings::IUPAC_CODE_MAP

decode_iupac <- function(code) {
  if (is.na(code) || !code %in% names(iupac)) return(character())
  strsplit(unname(iupac[code]), "", fixed = TRUE)[[1L]]
}
complement_alleles <- function(a) chartr("ACGT", "TGCA", toupper(a))

map_coordinates_to_rsid <- function(x, label) {
  x <- data.table::copy(data.table::as.data.table(x))
  req <- c("SNP", "CHR", "BP", "A1", "A2")
  missing <- setdiff(req, names(x))
  if (length(missing)) stop(label, " cannot map coordinate IDs; missing: ", paste(missing, collapse = ", "))
  chr <- sub("^chr", "", as.character(x$CHR), ignore.case = TRUE)
  bp <- as.integer(x$BP)
  seqname <- ifelse(chr %in% available_seqnames, chr, paste0("chr", chr))
  valid <- !is.na(bp) & bp > 0 & seqname %in% available_seqnames
  rows <- which(valid)
  if (!length(rows)) return(x[0])

  query <- GenomicRanges::GRanges(
    seqnames = seqname[rows],
    ranges = IRanges::IRanges(start = bp[rows], end = bp[rows])
  )
  hits <- BSgenome::snpsByOverlaps(dbsnp, query)
  if (!length(hits)) return(x[0])
  overlaps <- GenomicRanges::findOverlaps(query, hits, ignore.strand = TRUE)
  if (!length(overlaps)) return(x[0])

  meta <- S4Vectors::mcols(hits)
  ids <- as.character(meta$RefSNP_id)
  ids <- ifelse(grepl("^rs", ids, ignore.case = TRUE), ids, paste0("rs", ids))
  allele_codes <- as.character(meta$alleles_as_ambig)
  qh <- S4Vectors::queryHits(overlaps)
  sh <- S4Vectors::subjectHits(overlaps)
  matched <- vector("list", length(rows))
  for (k in seq_along(qh)) {
    i <- qh[[k]]
    j <- sh[[k]]
    ref_alleles <- decode_iupac(allele_codes[[j]])
    gwas_alleles <- toupper(c(x$A1[rows[[i]]], x$A2[rows[[i]]]))
    direct <- length(ref_alleles) == 2L && identical(sort(ref_alleles), sort(gwas_alleles))
    flipped <- length(ref_alleles) == 2L &&
      identical(sort(ref_alleles), sort(complement_alleles(gwas_alleles)))
    if (direct || flipped) matched[[i]] <- c(matched[[i]], ids[[j]])
  }
  keep <- integer()
  mapped_ids <- character()
  for (i in seq_along(rows)) {
    these_ids <- unique(matched[[i]])
    if (length(these_ids) == 1L) {
      keep <- c(keep, rows[[i]])
      mapped_ids <- c(mapped_ids, these_ids)
    }
  }
  out <- x[keep]
  if (nrow(out)) out[, SNP := mapped_ids]
  message(label, ": mapped ", nrow(out), " of ", length(rows),
          " coordinate/allele candidates to unique rsIDs.")
  out
}

add_grch38_coordinates_from_rsid <- function(x) {
  x <- data.table::copy(data.table::as.data.table(x))
  ids <- unique(as.character(x$SNP[grepl("^rs[0-9]+$", x$SNP, ignore.case = TRUE)]))
  if (!length(ids)) return(x[0])
  hits <- BSgenome::snpsById(dbsnp, ids, ifnotfound = "drop")
  if (!length(hits)) return(x[0])
  meta <- S4Vectors::mcols(hits)
  hit_ids <- as.character(meta$RefSNP_id)
  hit_ids <- ifelse(grepl("^rs", hit_ids, ignore.case = TRUE), hit_ids, paste0("rs", hit_ids))
  coords <- data.table(
    SNP = hit_ids,
    map_chr = as.character(GenomicRanges::seqnames(hits)),
    map_bp = as.numeric(GenomicRanges::start(hits))
  )
  coords <- unique(coords, by = "SNP")
  x[, .row_id := .I]
  x <- merge(x, coords, by = "SNP", all.x = TRUE, sort = FALSE)
  x[, CHR := map_chr]
  x[, BP := map_bp]
  x[, c("map_chr", "map_bp") := NULL]
  data.table::setorder(x, .row_id)
  x[, .row_id := NULL]
  x
}

bim_ids <- fread(paste0(LD_BFILE, ".bim"), select = 2L, header = FALSE,
                 colClasses = "character", showProgress = FALSE)[[1L]]

clump_exposure <- function(trait, threshold_name, p_threshold) {
  cache_path <- file.path(instrument_cache_dir,
    paste0(safe_file(trait), "__", threshold_name, ".rds"))
  if (file.exists(cache_path)) {
    message("Reusing cached clumped instruments: ", cache_path)
    return(readRDS(cache_path))
  }
  x <- load_trait(trait)
  iv <- x[is.finite(P) & P < p_threshold]
  message(trait, " variants passing P < ", format(p_threshold, scientific = TRUE), ": ", nrow(iv))
  if (!nrow(iv)) stop("No variants pass threshold for ", trait)

  needs_mapping <- !grepl("^rs[0-9]+$", iv$SNP, ignore.case = TRUE)
  if (any(needs_mapping)) {
    if (!identical(TRAIT_COORD_BUILD[[trait]], "GRCh38")) {
      stop("Coordinate-to-rsID mapping is only enabled for GRCh38 inputs; ", trait,
           " is marked ", TRAIT_COORD_BUILD[[trait]])
    }
    mapped <- map_coordinates_to_rsid(iv[needs_mapping], trait)
    iv <- data.table::rbindlist(list(iv[!needs_mapping], mapped), fill = TRUE)
  }

  iv <- iv[!is.na(SNP) & grepl("^rs[0-9]+$", SNP, ignore.case = TRUE) & SNP %chin% bim_ids]
  message(trait, " candidate rsIDs present in PLINK panel: ", nrow(iv))
  if (!nrow(iv)) stop("No candidate rsIDs present in PLINK panel for ", trait)
  duplicate_ids <- iv[, .N, by = SNP][N > 1L, SNP]
  if (length(duplicate_ids)) iv <- iv[!SNP %chin% duplicate_ids]

  clumped <- ieugwasr::ld_clump_local(
    dat = data.frame(rsid = iv$SNP, pval = iv$P),
    clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = p_threshold,
    bfile = LD_BFILE, plink_bin = PLINK_BIN
  )
  if (is.null(clumped) || !nrow(clumped)) stop("No clumped instruments for ", trait)
  kept <- iv[SNP %chin% clumped$rsid]
  if (!nrow(kept)) stop("No instruments remain after clumping for ", trait)
  saveRDS(kept, cache_path)
  message("Independent instruments after clumping for ", trait, ": ", nrow(kept))
  kept
}

to_tsmr <- function(x, type, trait_name) {
  if (!nrow(x)) stop("No variants to format for ", trait_name, " (", type, ").")
  d <- as.data.frame(x)
  d$trait_name <- trait_name
  args <- list(
    dat = d, type = type, snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "A1", other_allele_col = "A2",
    pval_col = "P", phenotype_col = "trait_name"
  )
  if ("EAF" %in% names(d) && any(is.finite(d$EAF))) args$eaf_col <- "EAF"
  if ("N" %in% names(d) && any(is.finite(d$N))) args$samplesize_col <- "N"
  do.call(TwoSampleMR::format_data, args)
}

match_outcome_rows <- function(exp_iv, outcome, exposure_name, outcome_name) {
  direct_rows <- outcome[which(outcome[["SNP"]] %chin% exp_iv[["SNP"]])]
  coordinate_rows <- data.table()

  # Only coordinate-match against GRCh38 positions. For rsID exposures,
  # retrieve their GRCh38 locations from dbSNP; this is a lookup, not liftover.
  if (identical(TRAIT_COORD_BUILD[[outcome_name]], "GRCh38") &&
      all(c("CHR", "BP", "A1", "A2") %in% names(outcome))) {
    e <- add_grch38_coordinates_from_rsid(exp_iv)
    if (all(c("CHR", "BP") %in% names(e)) && nrow(e)) {
      e[, `:=`(
        match_chr = toupper(sub("^chr", "", as.character(CHR), ignore.case = TRUE)),
        match_bp = as.numeric(BP),
        allele_lo = pmin(toupper(A1), toupper(A2)),
        allele_hi = pmax(toupper(A1), toupper(A2)),
        exposure_snp = as.character(SNP)
      )]
      o <- data.table::copy(outcome)
      o[, `:=`(
        match_chr = toupper(sub("^chr", "", as.character(CHR), ignore.case = TRUE)),
        match_bp = as.numeric(BP),
        allele_lo = pmin(toupper(A1), toupper(A2)),
        allele_hi = pmax(toupper(A1), toupper(A2))
      )]
      exposure_keys <- unique(e[, .(match_chr, match_bp, allele_lo, allele_hi, exposure_snp)])
      coordinate_rows <- merge(
        o, exposure_keys,
        by = c("match_chr", "match_bp", "allele_lo", "allele_hi"),
        all = FALSE, allow.cartesian = FALSE
      )
      if (nrow(coordinate_rows)) {
        coordinate_rows[, SNP := exposure_snp]
        coordinate_rows[, c("match_chr", "match_bp", "allele_lo", "allele_hi", "exposure_snp") := NULL]
      }
    }
  }

  if (nrow(direct_rows) && nrow(coordinate_rows)) {
    coordinate_rows <- coordinate_rows[!SNP %chin% direct_rows[["SNP"]]]
  }
  rows <- rbindlist(list(direct_rows, coordinate_rows), fill = TRUE)
  if (nrow(rows)) rows <- unique(rows, by = "SNP")
  message(exposure_name, " -> ", outcome_name, ": ", nrow(direct_rows),
          " direct rsID matches; ", nrow(coordinate_rows), " added GRCh38 coordinate matches.")
  rows
}

method_ids <- c(
  "Inverse variance weighted" = "mr_ivw",
  "MR Egger" = "mr_egger_regression",
  "Weighted median" = "mr_weighted_median",
  "Weighted mode" = "mr_weighted_mode"
)

run_one_pair <- function(exposure_name, outcome_name, threshold_name, p_threshold) {
  pair_label <- paste(exposure_name, "to", outcome_name, threshold_name, sep = "__")
  cache_path <- file.path(cache_dir, paste0(safe_file(pair_label), ".rds"))
  if (file.exists(cache_path)) {
    message("Reusing cached pair: ", pair_label)
    return(readRDS(cache_path))
  }

  message("\n--- ", pair_label, " ---")
  answer <- tryCatch({
    exp_iv <- clump_exposure(exposure_name, threshold_name, p_threshold)
    outcome <- load_trait(outcome_name)
    outcome_rows <- match_outcome_rows(exp_iv, outcome, exposure_name, outcome_name)
    if (!nrow(outcome_rows)) {
      stop("No outcome variants matched clumped exposure instruments.")
    }

    exp_dat <- to_tsmr(exp_iv, "exposure", exposure_name)
    out_dat <- to_tsmr(outcome_rows, "outcome", outcome_name)
    harm <- TwoSampleMR::harmonise_data(exp_dat, out_dat, action = 3)
    harm <- subset(harm, mr_keep)
    n_harmonized <- nrow(harm)
    if (!n_harmonized) stop("No variants remained after allele harmonization.")

    main_rows <- list()
    n_by_method <- if (n_harmonized == 1L) {
      c("Wald ratio" = "mr_wald_ratio")
    } else {
      method_ids
    }
    for (method_label in names(n_by_method)) {
      method_id <- unname(n_by_method[[method_label]])
      if (identical(method_id, "mr_egger_regression") && n_harmonized < 3L) {
        main_rows[[method_label]] <- data.table(
          exposure = exposure_name, outcome = outcome_name,
          threshold = threshold_name, p_threshold = p_threshold,
          nsnp = n_harmonized, method = method_label,
          b = NA_real_, se = NA_real_, pval = NA_real_,
          result_status = "not_estimable: fewer than 3 instruments"
        )
        next
      }
      fit <- tryCatch(
        TwoSampleMR::mr(harm, method_list = method_id),
        error = function(e) e
      )
      if (inherits(fit, "error") || !nrow(fit)) {
        err <- if (inherits(fit, "error")) conditionMessage(fit) else "method returned no result"
        main_rows[[method_label]] <- data.table(
          exposure = exposure_name, outcome = outcome_name,
          threshold = threshold_name, p_threshold = p_threshold,
          nsnp = n_harmonized, method = method_label,
          b = NA_real_, se = NA_real_, pval = NA_real_, result_status = err
        )
      } else {
        r <- as.data.table(fit)
        r[, `:=`(exposure = exposure_name, outcome = outcome_name,
                 threshold = threshold_name, p_threshold = p_threshold,
                 result_status = "completed")]
        main_rows[[method_label]] <- r
      }
    }
    main <- rbindlist(main_rows, fill = TRUE, use.names = TRUE)

    sensitivity_rows <- list()
    if (n_harmonized >= 2L) {
      het <- tryCatch(TwoSampleMR::mr_heterogeneity(harm), error = function(e) e)
      if (!inherits(het, "error") && nrow(het)) {
        h <- as.data.table(het)
        h[, `:=`(exposure = exposure_name, outcome = outcome_name,
                 threshold = threshold_name, p_threshold = p_threshold,
                 sensitivity = "Cochran Q", result_status = "completed")]
        sensitivity_rows[["heterogeneity"]] <- h
      } else {
        sensitivity_rows[["heterogeneity"]] <- data.table(
          exposure = exposure_name, outcome = outcome_name,
          threshold = threshold_name, p_threshold = p_threshold,
          sensitivity = "Cochran Q", result_status =
            if (inherits(het, "error")) conditionMessage(het) else "no result"
        )
      }

      loo <- tryCatch(TwoSampleMR::mr_leaveoneout(harm), error = function(e) e)
      if (!inherits(loo, "error") && nrow(loo)) {
        l <- as.data.table(loo)
        l[, `:=`(exposure = exposure_name, outcome = outcome_name,
                 threshold = threshold_name, p_threshold = p_threshold,
                 sensitivity = "Leave one out IVW", result_status = "completed")]
        sensitivity_rows[["loo"]] <- l
      } else {
        sensitivity_rows[["loo"]] <- data.table(
          exposure = exposure_name, outcome = outcome_name,
          threshold = threshold_name, p_threshold = p_threshold,
          sensitivity = "Leave one out IVW", result_status =
            if (inherits(loo, "error")) conditionMessage(loo) else "no result"
        )
      }
    }

    if (n_harmonized >= 3L) {
      egger <- tryCatch(TwoSampleMR::mr_pleiotropy_test(harm), error = function(e) e)
      if (!inherits(egger, "error") && nrow(egger)) {
        p <- as.data.table(egger)
        p[, `:=`(exposure = exposure_name, outcome = outcome_name,
                 threshold = threshold_name, p_threshold = p_threshold,
                 sensitivity = "Egger intercept", result_status = "completed")]
        sensitivity_rows[["egger"]] <- p
      } else {
        sensitivity_rows[["egger"]] <- data.table(
          exposure = exposure_name, outcome = outcome_name,
          threshold = threshold_name, p_threshold = p_threshold,
          sensitivity = "Egger intercept", result_status =
            if (inherits(egger, "error")) conditionMessage(egger) else "no result"
        )
      }
    } else {
      sensitivity_rows[["egger"]] <- data.table(
        exposure = exposure_name, outcome = outcome_name,
        threshold = threshold_name, p_threshold = p_threshold,
        sensitivity = "Egger intercept",
        result_status = "not_estimable: fewer than 3 instruments"
      )
    }

    list(
      status = data.table(
        exposure = exposure_name, outcome = outcome_name,
        threshold = threshold_name, p_threshold = p_threshold,
        coordinate_build_exposure = TRAIT_COORD_BUILD[[exposure_name]],
        coordinate_build_outcome = TRAIT_COORD_BUILD[[outcome_name]],
        n_clumped = nrow(exp_iv), n_outcome_matches = nrow(outcome_rows),
        n_harmonized = n_harmonized, run_status = "completed", error = NA_character_
      ),
      main = main,
      sensitivity = rbindlist(sensitivity_rows, fill = TRUE, use.names = TRUE),
      harmonized = as.data.table(harm)[, `:=`(
        exposure = exposure_name,
        outcome = outcome_name,
        threshold = threshold_name,
        p_threshold = p_threshold
      )]
    )
  }, error = function(e) {
    message("FAILED ", pair_label, ": ", conditionMessage(e))
    list(
      status = data.table(
        exposure = exposure_name, outcome = outcome_name,
        threshold = threshold_name, p_threshold = p_threshold,
        coordinate_build_exposure = TRAIT_COORD_BUILD[[exposure_name]],
        coordinate_build_outcome = TRAIT_COORD_BUILD[[outcome_name]],
        n_clumped = NA_integer_, n_outcome_matches = NA_integer_,
        n_harmonized = NA_integer_, run_status = "failed",
        error = conditionMessage(e)
      ),
      main = data.table(), sensitivity = data.table(), harmonized = data.table()
    )
  })

  saveRDS(answer, cache_path)
  answer
}

# Run every direction between the loneliness datasets and the three core traits,
# plus both directions among the core traits. This gives 24 directed comparisons.
loneliness_core_edges <- CJ(
  exposure = LONELINESS_TRAITS,
  outcome = CORE_TRAITS,
  sorted = FALSE
)
core_edges <- as.data.table(t(combn(CORE_TRAITS, 2L)))
setnames(core_edges, c("exposure", "outcome"))
edges <- rbindlist(
  list(loneliness_core_edges, core_edges),
  use.names = TRUE
)
directed <- rbindlist(list(
  copy(edges), edges[, .(exposure = outcome, outcome = exposure)]
), use.names = TRUE)

all_results <- list()
total <- nrow(directed) * length(THRESHOLDS)
run_i <- 0L
for (threshold_name in names(THRESHOLDS)) {
  p_threshold <- unname(THRESHOLDS[[threshold_name]])
  for (i in seq_len(nrow(directed))) {
    run_i <- run_i + 1L
    message("\nProgress ", run_i, "/", total)
    all_results[[paste(threshold_name, i, sep = "_")]] <- run_one_pair(
      directed$exposure[[i]], directed$outcome[[i]], threshold_name, p_threshold
    )

    main_table <- rbindlist(lapply(all_results, `[[`, "main"), fill = TRUE, use.names = TRUE)
    sensitivity_table <- rbindlist(lapply(all_results, `[[`, "sensitivity"), fill = TRUE, use.names = TRUE)
    status_table <- rbindlist(lapply(all_results, `[[`, "status"), fill = TRUE, use.names = TRUE)
    fwrite(main_table, file.path(out_dir, "twosamplemr_main_results.csv"))
    fwrite(sensitivity_table, file.path(out_dir, "twosamplemr_sensitivity_results.csv"))
    fwrite(status_table, file.path(out_dir, "twosamplemr_status.csv"))
  }
}

saveRDS(all_results, file.path(out_dir, "twosamplemr_all_pair_results.rds"),
        compress = FALSE)
message("All requested TwoSampleMR analyses finished. Results: ", out_dir)
