#!/usr/bin/env Rscript

# Primary sumstats cleaner: reads the reference sheet, standardizes each file,
# and writes a manifest for the LDSC step that follows.
#
# Usage: clean_sumstats.R <mapping CSV/URL> <output dir> [already_munged.csv]
#                         [analysis] [qc_passed.tsv]
#
# Environment variables:
#   CLEAN_CORES     number of parallel workers across phenotypes (default 1).
#                   Each worker loads its own input file, so keep this modest
#                   for very large sources.
#   CLEAN_CHUNK_ROWS  rows per batch in the streamed coordinate-mapping path
#                   (default 250000).
#   CLEAN_BATCH_GB  input files larger than this (GB) that need coordinate
#                   rsID mapping are processed through the resumable batched
#                   path instead of loading everything into memory (default 2).
#
# Skip logic (a phenotype is never reprocessed when any of these apply):
#   1. Listed in the optional qc_passed.tsv registry (phenotypes whose cleaned
#      output already passed the deep QC checks).
#   2. Listed in the already_munged.csv registry.
#   3. A complete cleaned output file already exists (resumable reruns).
# Per-trait failures are caught and written to clean_failures.tsv so one bad
# source cannot abort a long batch run; failed traits are retried on rerun.
#
# Version 2 changes (after first QC pass):
#   - Column matching picks the alias-matching column with the most usable
#     values, so an empty "A1" column no longer shadows a populated "Allele1"
#     (BMI GIANT) or an empty "A2" a populated "other_allele" (Depression PGC).
#   - chr:pos_ref_alt / chr_pos_ref_alt style SNP IDs (Crohn, GlycA, Tyrosine,
#     Hypothyroidism) are parsed into CHR/BP/A1/A2 and routed through
#     coordinate-based rsID mapping.
#   - Coordinate mapping now triggers when fewer than 50% of SNP IDs are
#     rsIDs, not only when the SNP column is entirely NA.
#   - CHR values are normalized ("chr10" -> "10"; GSCAN files).
#   - Alleles are uppercased (SCALLOP/omicscience files) and rows with
#     non-ACGT allele codes (daner I/D indel codes; Bipolar PGC) are dropped.
#   - Rows with FRQ outside (0, 1) are dropped (Atopic dermatitis, T2D).
#   - Warns when median N disagrees with median N_CASE + N_CONTROL by >25%.
#   - Google Sheets edit URLs are converted to CSV export URLs automatically.
#   - Optional parallel processing across phenotypes (CLEAN_CORES).

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Command-line arguments and settings
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: clean_sumstats.R <mapping CSV or URL> <output directory> ",
       "[already_munged.csv] [analysis] [qc_passed.tsv]")
}
mapping_path <- args[[1L]]
out_dir <- args[[2L]]
already_munged_path <- if (length(args) >= 3L) args[[3L]] else ""
analysis <- if (length(args) >= 4L) args[[4L]] else ""
qc_passed_path <- if (length(args) >= 5L) args[[5L]] else ""
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- suppressWarnings(as.integer(Sys.getenv("CLEAN_CORES", "1")))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
batch_threshold_gb <- suppressWarnings(as.numeric(Sys.getenv("CLEAN_BATCH_GB", "2")))
if (is.na(batch_threshold_gb) || batch_threshold_gb < 0) batch_threshold_gb <- 5

message("Primary cleaner v2: best-column matching, pseudo-ID parsing, rsID-fraction ",
        "mapping trigger, allele/FRQ cleaning, resumable batches, CLEAN_CORES=", n_cores)

# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
safe_name <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
as_num <- function(x) {
  if (is.null(x)) return(NULL)
  suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
}
fill_value <- function(x, n, type = c("character", "numeric")) {
  type <- match.arg(type)
  if (is.null(x)) return(if (type == "character") rep(NA_character_, n) else rep(NA_real_, n))
  if (type == "character") as.character(x) else as_num(x)
}
first_match <- function(headers, alias_names) {
  hits <- find_col_matches(headers, alias_names)
  if (length(hits)) hits[[1L]] else NA_character_
}

# Google Sheets edit/share URLs must be converted to their CSV export endpoint
# before fread can download them.
as_export_url <- function(path) {
  if (!grepl("docs\\.google\\.com/spreadsheets", path)) return(path)
  if (grepl("format=csv", path, fixed = TRUE)) return(path)
  sheet_id <- sub(".*spreadsheets/d/([^/]+).*", "\\1", path)
  gid <- if (grepl("gid=[0-9]+", path)) sub(".*gid=([0-9]+).*", "\\1", path) else "0"
  sprintf("https://docs.google.com/spreadsheets/d/%s/export?format=csv&gid=%s", sheet_id, gid)
}

# PGC daner files begin with ## metadata lines; find the header row without
# decoding arbitrary leading bytes as UTF-8.
header_skip <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, n = 200L, warn = FALSE, encoding = "bytes")
  ascii <- vapply(lines, function(line) {
    bytes <- as.integer(charToRaw(line))
    rawToChar(as.raw(bytes[bytes < 128L]))
  }, character(1))
  first_data_line <- which(nzchar(trimws(ascii)) & !grepl("^##", trimws(ascii)))[1L]
  if (is.na(first_data_line)) stop("Could not find a header row in: ", path)
  first_data_line - 1L
}

# ---------------------------------------------------------------------------
# Column matching
# ---------------------------------------------------------------------------

# All header names matching any alias, in alias-preference order.
find_col_matches <- function(headers, alias_names) {
  hits <- match(norm(alias_names), norm(headers), nomatch = 0L)
  hits <- hits[hits > 0L]
  if (length(hits)) headers[hits] else character()
}

# Among all alias-matching columns, take the one with the most usable values.
# This prevents an empty column that happens to carry a preferred alias
# (for example an all-NA "A1") from shadowing a populated alternative
# ("Allele1", "effect_allele", ...).
best_column <- function(headers, values, alias_names) {
  candidates <- find_col_matches(headers, alias_names)
  if (!length(candidates)) return(NULL)
  scores <- vapply(candidates, function(column) {
    x <- values[[column]]
    sum(!is.na(x) & nzchar(trimws(as.character(x))))
  }, integer(1))
  candidates[[which.max(scores)]]
}

# Common column names used by the included GWAS sources.
aliases <- list(
  SNP = c("rsids", "rsid", "rs_id", "SNP", "SNPID", "MarkerName", "variant_id", "variant", "ID"),
  CHR = c("#chrom", "CHR", "chromosome", "chromsome", "chrom", "chr"),
  BP = c("POS", "BP", "base_pair_location", "position", "pos"),
  A1 = c("effect_allele", "A1", "Allele1", "Tested_Allele", "ALT", "inc_allele", "EA"),
  A2 = c("other_allele", "A2", "Allele2", "Other_Allele", "NEA", "REF", "non_effect_allele", "OA"),
  BETA = c("beta", "BETA", "Effect", "B", "log_odds"),
  OR = c("OR", "odds_ratio"),
  SE = c("standard_error", "SE", "sebeta", "StdErr", "stderr"),
  P = c("p_value", "P", "pval", "P.value", "Pvalue"),
  LOGP = c("neg_log_10_p_value", "mlogp", "minus_log10_p", "neg_log10_p", "LOG10P"),
  Z = c("Z", "Zscore", "Z_score", "z_score", "tstat"),
  DIRECTION = c("direction", "effect_direction", "beta_direction"),
  FRQ = c("effect_allele_frequency", "EAF", "Freq_Tested_Allele_in_HRS", "af_alt",
          "AF_Allele2", "FRQ", "af"),
  INFO = c("INFO", "imputation_info", "IMPINFO"),
  N = c("N", "n", "sample_size", "total_n", "nCompleteSamples", "Neff"),
  N_CASE = c("N_CASE", "N_CASES", "N_CAS", "Nca", "cases_n"),
  N_CONTROL = c("N_CONTROL", "N_CONTROLS", "N_CON", "Nco", "controls_n")
)

# ---------------------------------------------------------------------------
# Variant-level cleaning (shared by the standard and batched paths)
# ---------------------------------------------------------------------------

# Parse 'chr:pos_ref_alt' / 'chr_pos_ref_alt' style IDs (Crohn disease uses
# colons, GWAS Catalog protein files use underscores).
parse_pseudo_id <- function(ids) {
  strcapture(
    "^(?:chr)?([0-9XYM]+)[:_]([0-9]+)[:_]([A-Za-z]+)[:_]([A-Za-z]+)$",
    as.character(ids),
    data.frame(chr = character(), bp = character(), a1 = character(), a2 = character())
  )
}

# Normalize CHR, uppercase alleles, parse pseudo-IDs into coordinates, and
# drop rows that cannot be used for LDSC munging.
clean_variants <- function(out, trait) {
  out[, CHR := toupper(sub("^chr", "", trimws(as.character(CHR)), ignore.case = TRUE))]
  out[CHR %in% c("23"), CHR := "X"]
  out[CHR %in% c("24"), CHR := "Y"]
  out[CHR %in% c("25", "M"), CHR := "MT"]

  out[, A1 := toupper(trimws(as.character(A1)))]
  out[, A2 := toupper(trimws(as.character(A2)))]

  # Recover coordinates and alleles from pseudo-IDs, then clear the ID so the
  # coordinate-based rsID mapping replaces it.
  pseudo <- !is.na(out$SNP) & !grepl("^rs[0-9]+$", out$SNP)
  if (any(pseudo)) {
    parsed <- parse_pseudo_id(out$SNP[pseudo])
    ok <- !is.na(parsed$chr) & !is.na(parsed$bp)
    if (any(ok)) {
      message(trait, ": parsing ", sum(ok), " chr:pos_alleles style ID(s) into coordinates")
      idx <- which(pseudo)[ok]
      chr_v <- as.character(out$CHR)
      bp_v <- out$BP
      a1_v <- out$A1
      a2_v <- out$A2
      chr_v[idx] <- ifelse(is.na(chr_v[idx]) | !nzchar(chr_v[idx]), parsed$chr[ok], chr_v[idx])
      bp_v[idx] <- ifelse(!is.finite(bp_v[idx]),
                          suppressWarnings(as.numeric(parsed$bp[ok])), bp_v[idx])
      a1_v[idx] <- ifelse(is.na(a1_v[idx]) | !nzchar(a1_v[idx]), toupper(parsed$a1[ok]), a1_v[idx])
      a2_v[idx] <- ifelse(is.na(a2_v[idx]) | !nzchar(a2_v[idx]), toupper(parsed$a2[ok]), a2_v[idx])
      out[, CHR := chr_v]
      out[, BP := bp_v]
      out[, A1 := a1_v]
      out[, A2 := a2_v]
      out$SNP[idx] <- NA_character_
    }
  }

  # daner I/D indel codes and any other non-ACGT allele strings cannot be
  # munged; multi-base alleles (for example "AG") are kept.
  bad_allele <- (!is.na(out$A1) & !grepl("^[ACGT]+$", out$A1)) |
    (!is.na(out$A2) & !grepl("^[ACGT]+$", out$A2))
  if (any(bad_allele)) {
    message("Dropping ", sum(bad_allele), " row(s) with non-ACGT allele codes for ", trait)
    out <- out[!bad_allele]
  }

  # Allele frequencies must be inside (0, 1); FRQ = 1 rows are malformed.
  if ("FRQ" %in% names(out)) {
    bad_frq <- !is.na(out$FRQ) & (out$FRQ <= 0 | out$FRQ >= 1)
    if (any(bad_frq)) {
      message("Dropping ", sum(bad_frq), " row(s) with FRQ outside (0, 1) for ", trait)
      out <- out[!bad_frq]
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Coordinate-based rsID mapping
# ---------------------------------------------------------------------------

map_coordinates_to_rsids <- function(out, trait, src_text) {
  coordinate_build <- if (grepl("finngen", src_text, ignore.case = TRUE)) {
    "GRCh38"
  } else if (grepl("GWAS Catalog|DIAGRAM", src_text, ignore.case = TRUE)) {
    "GRCh37"
  } else {
    stop("Pseudo-rsIDs or missing rsIDs for ", trait, " but source/cohort ('",
         src_text, "') does not identify a genome build for coordinate lookup.")
  }
  required_reference_packages <- if (coordinate_build == "GRCh37") {
    c("MungeSumstats", "SNPlocs.Hsapiens.dbSNP155.GRCh37",
      "BSgenome.Hsapiens.1000genomes.hs37d5")
  } else {
    c("MungeSumstats", "SNPlocs.Hsapiens.dbSNP155.GRCh38",
      "BSgenome.Hsapiens.NCBI.GRCh38")
  }
  missing_reference_packages <- required_reference_packages[
    !vapply(required_reference_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_reference_packages)) {
    stop("Cannot map coordinates to rsIDs for ", trait, ". Missing package(s): ",
         paste(missing_reference_packages, collapse = ", "),
         ". Install them in R_LIBS_USER before running this job.")
  }
  for (package in required_reference_packages) {
    message(package, " version: ", as.character(utils::packageVersion(package)))
  }
  message("Mapping rsIDs from ", coordinate_build, " coordinates for ", trait)
  rows_with_coordinates <- !is.na(out$CHR) & is.finite(out$BP) & out$BP > 0
  if (sum(!rows_with_coordinates)) {
    message("Dropping ", sum(!rows_with_coordinates), " row(s) without usable coordinates for ", trait)
    out <- out[rows_with_coordinates]
  }
  out[, SNP := NULL]
  mapped_path <- tempfile(tmpdir = out_dir, fileext = ".tsv.gz")
  mapped <- MungeSumstats::format_sumstats(
    path = out,
    ref_genome = coordinate_build,
    snp_ids_are_rs_ids = FALSE,
    dbSNP = 155,
    return_data = TRUE,
    return_format = "data.table",
    save_path = mapped_path,
    force_new = TRUE,
    log_folder_ind = FALSE,
    sort_coordinates = FALSE
  )
  unlink(mapped_path)
  mapped <- as.data.table(mapped)
  if (!"SNP" %in% names(mapped) || all(is.na(mapped$SNP))) {
    stop("Could not map coordinates to rsIDs for ", trait,
         ". Check that the file build matches ", coordinate_build,
         " and that its matching SNPlocs reference is installed.")
  }
  mapped
}

# ---------------------------------------------------------------------------
# Reference sheet and skip registries
# ---------------------------------------------------------------------------

ref <- fread(as_export_url(mapping_path), na.strings = c("", "NA", "N/A"), showProgress = FALSE)
ph_col <- first_match(names(ref), "phenotype")
path_col <- first_match(names(ref), "file location")
if (is.na(ph_col) || is.na(path_col)) {
  stop("The reference sheet must have 'phenotype' and 'file location' columns.")
}
source_col <- first_match(names(ref), "source")
cohort_col <- first_match(names(ref), "cohort")
analysis_col <- first_match(names(ref), "analysis")
n_col <- first_match(names(ref), c("n", "sample size", "total n"))
n_case_col <- first_match(names(ref), c("n case", "n cases", "n_cas"))
n_control_col <- first_match(names(ref), c("n control", "n controls", "n_con"))

if (nzchar(analysis)) {
  if (is.na(analysis_col)) stop("The reference sheet needs an 'analysis' column for this run.")
  analysis_values <- toupper(trimws(as.character(ref[[analysis_col]])))
  requested_analysis <- toupper(trimws(analysis))
  ref <- ref[which(!is.na(analysis_values) & analysis_values == requested_analysis)]
  message("Using ", nrow(ref), " row(s) marked analysis = ", analysis)
}

manifest <- data.table(
  phenotype = as.character(ref[[ph_col]]),
  analysis = if (!is.na(analysis_col)) as.character(ref[[analysis_col]]) else "",
  source = if (!is.na(source_col)) as.character(ref[[source_col]]) else "",
  cohort = if (!is.na(cohort_col)) as.character(ref[[cohort_col]]) else "",
  input_path = trimws(as.character(ref[[path_col]]))
)
manifest <- manifest[!is.na(phenotype) & nzchar(phenotype) & !is.na(input_path) & nzchar(input_path)]
if (nzchar(analysis) &&
    any(toupper(trimws(manifest$analysis)) != requested_analysis, na.rm = TRUE)) {
  stop("Internal filter error: non-", analysis, " rows entered the primary manifest.")
}
fwrite(manifest, file.path(out_dir, "primary_manifest.tsv"), sep = "\t", na = "")
message("Primary manifest written: ", file.path(out_dir, "primary_manifest.tsv"),
        " (", nrow(manifest), " row(s); analysis = ", analysis, ")")

# Registry 1: phenotypes whose cleaned output already passed deep QC.
# sep is pinned to tab because phenotype names contain spaces.
qc_passed <- character()
if (nzchar(qc_passed_path)) {
  qc_passed <- as.character(fread(qc_passed_path, sep = "\t", showProgress = FALSE)$phenotype)
  message("QC-passed registry: ", length(qc_passed), " phenotype(s) will be skipped")
}

# Registry 2: phenotypes that are already LDSC-munged downstream.
already_munged <- character()
if (nzchar(already_munged_path)) {
  previous <- fread(already_munged_path, na.strings = c("", "NA", "N/A"), showProgress = FALSE)
  already_munged <- as.character(previous[[first_match(names(previous), "phenotype")]])
}

# A FinnGen R13 URL takes precedence over an older R12 pre-munged file with
# the same exact label.
r13 <- grepl("finngen-public-data-r13", manifest$input_path, ignore.case = TRUE)
already_munged <- setdiff(already_munged, manifest$phenotype[r13])

# Download URL inputs once into scratch; use local file paths unchanged.
resolve_path <- function(path, phenotype) {
  if (!grepl("^https?://", path, ignore.case = TRUE)) return(path)
  cache_dir <- file.path(out_dir, "raw_downloads")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  filename <- basename(sub("[?#].*$", "", path))
  target <- file.path(cache_dir, paste0(safe_name(phenotype), "__", filename))
  if (!file.exists(target)) {
    temporary <- tempfile(tmpdir = cache_dir)
    status <- system2("curl", c("--location", "--fail", "--silent", "--show-error",
                                "--retry", "3", "--output", temporary, path))
    if (!identical(as.integer(status), 0L) || !file.exists(temporary) ||
        file.info(temporary)$size == 0) {
      unlink(temporary)
      stop("Download failed: ", path)
    }
    if (!file.rename(temporary, target)) stop("Could not cache download: ", target)
  }
  target
}

# ---------------------------------------------------------------------------
# Batched coordinate mapping for very large GRCh37 files without rsIDs
# (streams the input in chunks so the whole file never resides in memory).
# ---------------------------------------------------------------------------

coordinate_batch_rows <- suppressWarnings(as.integer(Sys.getenv("CLEAN_CHUNK_ROWS", "250000")))
if (is.na(coordinate_batch_rows) || coordinate_batch_rows < 1L) {
  stop("CLEAN_CHUNK_ROWS must be a positive integer.")
}

clean_coordinates_in_batches <- function(path, skip, file_headers,
                                         read_headers, trait, out_path, sheet_n) {
  required_packages <- c(
    "SNPlocs.Hsapiens.dbSNP155.GRCh37",
    "BSgenome.Hsapiens.1000genomes.hs37d5",
    "BSgenome",
    "GenomicRanges",
    "IRanges",
    "S4Vectors"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Coordinate-based rsID lookup needs package(s): ",
         paste(missing_packages, collapse = ", "),
         ". Install them in R_LIBS_USER before running this job.")
  }
  snplocs_package <- "SNPlocs.Hsapiens.dbSNP155.GRCh37"
  snplocs <- getExportedValue(snplocs_package, snplocs_package)

  info <- file.info(path)
  signature_fields <- list(path = normalizePath(path), size = info$size,
                           mtime = as.numeric(info$mtime), headers = file_headers,
                           selected = read_headers, rows_per_batch = coordinate_batch_rows,
                           sample_size = sheet_n, build = "GRCh37",
                           mungesumstats = as.character(utils::packageVersion("MungeSumstats")))
  signature_path <- tempfile(tmpdir = out_dir)
  on.exit(unlink(signature_path), add = TRUE)
  trait_cache_dir <- file.path(out_dir, "batch_cache", safe_name(trait))
  dir.create(trait_cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Look for caches made by earlier cleaner versions before starting a new
  # folder, so completed lymphocyte batches from v1 runs are still reused. The
  # rest of the signature must still match exactly (input metadata, selected
  # columns, sample size, rows per batch).
  signature_versions <- 2:4
  cache_keys <- vapply(signature_versions, function(version) {
    fields <- signature_fields
    if (version <= 3L) {
      signature <- c(list(version = as.integer(version)), fields)
    } else {
      fields$cleaner <- "v2"
      signature <- c(list(version = as.integer(version)), fields)
    }
    writeBin(serialize(signature, NULL), signature_path)
    unname(tools::md5sum(signature_path))
  }, character(1))
  candidate_dirs <- file.path(trait_cache_dir, cache_keys)
  candidate_exists <- dir.exists(candidate_dirs)
  completed_counts <- vapply(candidate_dirs, function(candidate) {
    if (!dir.exists(candidate)) return(0L)
    length(Sys.glob(file.path(candidate, "chunk_*.tsv.gz"))) +
      length(Sys.glob(file.path(candidate, "chunk_*.empty")))
  }, integer(1))
  reusable <- which(candidate_exists & completed_counts > 0L)
  if (length(reusable)) {
    best <- reusable[order(completed_counts[reusable],
                           signature_versions[reusable] == 4L, decreasing = TRUE)][[1L]]
    cache_key <- cache_keys[[best]]
    message("Reusing compatible cache folder for ", trait, ": ",
            candidate_dirs[[best]], " (", completed_counts[[best]], " completed batches)")
  } else {
    cache_key <- cache_keys[[which(signature_versions == 4L)]]
    message("No compatible completed cache found for ", trait,
            "; starting a new cache at ", file.path(trait_cache_dir, cache_key))
  }
  cache_dir <- file.path(trait_cache_dir, cache_key)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  input_con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(input_con), add = TRUE)
  if (skip > 0L) readLines(input_con, n = skip, warn = FALSE)
  header_line <- readLines(input_con, n = 1L, warn = FALSE)
  if (!length(header_line)) stop("Could not read the header for ", trait)
  raw_header <- strsplit(sub("\\r$", "", header_line), "\t", fixed = TRUE)[[1L]]
  raw_header[[1L]] <- sub("^\\ufeff", "", raw_header[[1L]])
  raw_header <- trimws(raw_header)
  selected_indices <- match(norm(read_headers), norm(raw_header))
  if (anyNA(selected_indices)) {
    stop("Could not match selected fields to the streamed header for ", trait,
         ". Expected: ", paste(read_headers, collapse = ", "),
         "; observed: ", paste(raw_header, collapse = ", "))
  }
  chunk_files <- character()
  chunk_number <- 0L
  total_rows <- 0
  repeat {
    lines <- readLines(input_con, n = coordinate_batch_rows, warn = FALSE)
    if (!length(lines)) break
    chunk_number <- chunk_number + 1L
    total_rows <- total_rows + length(lines)
    chunk_path <- file.path(cache_dir, sprintf("chunk_%06d.tsv.gz", chunk_number))
    empty_path <- file.path(cache_dir, sprintf("chunk_%06d.empty", chunk_number))
    chunk_ok <- file.exists(chunk_path) && file.info(chunk_path)$size > 0 &&
      identical(as.integer(system2("gzip", c("-t", chunk_path), stdout = FALSE, stderr = FALSE)), 0L)

    if (chunk_ok || file.exists(empty_path)) {
      message("Reusing cached coordinate-lookup batch ", chunk_number, " for ", trait)
    } else {
      dt <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE,
                  col.names = raw_header[selected_indices], select = selected_indices,
                  check.names = FALSE, showProgress = FALSE)
      n <- nrow(dt)
      get_value <- function(key) {
        column <- best_column(names(dt), dt, aliases[[key]])
        if (is.null(column)) NULL else dt[[column]]
      }
      beta <- as_num(get_value("BETA"))
      odds_ratio <- as_num(get_value("OR"))
      if (is.null(beta) && !is.null(odds_ratio)) beta <- log(odds_ratio)
      se <- as_num(get_value("SE"))
      p <- as_num(get_value("P"))
      logp <- as_num(get_value("LOGP"))
      z <- as_num(get_value("Z"))
      if (is.null(p) && !is.null(logp)) p <- 10^(-logp)
      if (is.null(z) && !is.null(beta) && !is.null(se)) z <- beta / se
      if (is.null(z) && !is.null(beta) && !is.null(p)) {
        z <- sign(beta) * qnorm(pmax(p / 2, .Machine$double.xmin), lower.tail = FALSE)
      }
      if (is.null(p) && !is.null(z)) p <- 2 * pnorm(-abs(z))
      if (!is.null(p)) {
        p[!is.na(p) & p <= 0] <- .Machine$double.xmin
        p[!is.na(p) & p > 1] <- NA_real_
      }
      if (!is.null(z)) z[!is.finite(z)] <- NA_real_

      sheet_n_value <- as_num(sheet_n)
      if (length(sheet_n_value) != 1L || is.na(sheet_n_value)) sheet_n_value <- NA_real_
      sample_n <- as_num(get_value("N"))
      if (is.null(sample_n)) sample_n <- rep(sheet_n_value, n)
      else sample_n[!is.finite(sample_n)] <- sheet_n_value

      out <- data.table(
        SNP = fill_value(get_value("SNP"), n),
        CHR = fill_value(get_value("CHR"), n),
        BP = fill_value(get_value("BP"), n, "numeric"),
        A1 = fill_value(get_value("A1"), n),
        A2 = fill_value(get_value("A2"), n),
        BETA = fill_value(beta, n, "numeric"),
        SE = fill_value(se, n, "numeric"),
        P = fill_value(p, n, "numeric"),
        Z = fill_value(z, n, "numeric"),
        FRQ = fill_value(get_value("FRQ"), n, "numeric"),
        INFO = fill_value(get_value("INFO"), n, "numeric"),
        N = fill_value(sample_n, n, "numeric")
      )
      out <- clean_variants(out, trait)
      out <- out[!is.na(CHR) & nzchar(CHR) & is.finite(BP) & BP > 0 &
                   !is.na(A1) & nzchar(A1) & !is.na(A2) & nzchar(A2) &
                   is.finite(BETA) & is.finite(SE) & SE > 0 &
                   is.finite(P) & P > 0 & P <= 1]
      missing_id <- is.na(out$SNP) | !grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE)
      rm(dt, lines)

      # Match only the distinct chromosome/position pairs in this chunk. The
      # dbSNP package returns rsIDs in RefSNP_id; positions with multiple
      # distinct rsIDs are treated as ambiguous and left unmapped.
      out[, lookup_chr := CHR]
      out[, lookup_bp := suppressWarnings(as.integer(BP))]
      out[, mapped_snp := NA_character_]
      query_loci <- unique(out[missing_id & !is.na(lookup_chr) &
                                 lookup_chr %in% c(as.character(1:22), "X", "Y", "MT") &
                                 !is.na(lookup_bp) & lookup_bp > 0,
                               .(lookup_chr, lookup_bp)])
      if (nrow(query_loci)) {
        query_ranges <- GenomicRanges::GRanges(
          seqnames = query_loci$lookup_chr,
          ranges = IRanges::IRanges(start = query_loci$lookup_bp, width = 1L)
        )
        hits <- BSgenome::snpsByOverlaps(snplocs, query_ranges)
        hit_ids <- as.character(S4Vectors::mcols(hits)$RefSNP_id)
        hit_table <- data.table(
          lookup_chr = as.character(GenomicRanges::seqnames(hits)),
          lookup_bp = as.integer(IRanges::start(hits)),
          SNP = hit_ids
        )[!is.na(SNP) & nzchar(SNP)]
        locus_map <- hit_table[, .(
          mapped_snp = if (uniqueN(SNP) == 1L) unique(SNP)[[1L]] else NA_character_,
          ambiguous = uniqueN(SNP) > 1L
        ), by = .(lookup_chr, lookup_bp)]
        out[locus_map, on = .(lookup_chr, lookup_bp), mapped_snp := i.mapped_snp]
        mapped_rows <- sum(missing_id & !is.na(out$mapped_snp))
        ambiguous_loci <- sum(locus_map$ambiguous)
        if (mapped_rows || ambiguous_loci) {
          message("Direct GRCh37 rsID lookup for ", trait, " batch ", chunk_number,
                  ": mapped ", mapped_rows, " rows; ", ambiguous_loci,
                  " coordinate(s) had multiple rsIDs and were left unmapped")
        }
        out[missing_id & !is.na(mapped_snp), SNP := mapped_snp]
        rm(query_ranges, hits, hit_ids, hit_table, locus_map)
      }
      out[, c("lookup_chr", "lookup_bp", "mapped_snp") := NULL]
      out <- out[!is.na(SNP) & nzchar(SNP)]
      if (!"Z" %in% names(out) && all(c("BETA", "SE") %in% names(out))) out[, Z := BETA / SE]
      if (!"N_CASE" %in% names(out)) out[, N_CASE := NA_real_]
      if (!"N_CONTROL" %in% names(out)) out[, N_CONTROL := NA_real_]
      out <- out[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, Z,
                     FRQ, INFO, N, N_CASE, N_CONTROL)]
      if (nrow(out)) {
        fwrite(out, chunk_path, sep = "\t", na = "NA")
        if (file.exists(empty_path)) unlink(empty_path)
      } else {
        file.create(empty_path)
        if (file.exists(chunk_path)) unlink(chunk_path)
      }
      rm(out)
      # Fully collect between batches; repeated reference lookups can otherwise
      # leave enough temporary objects resident to exhaust the Slurm memory.
      invisible(gc())
    }
    if (file.exists(chunk_path)) chunk_files <- c(chunk_files, chunk_path)
    rm(lines)
    if (chunk_number %% 10L == 0L) {
      message("Processed coordinate-lookup batch ", chunk_number, " for ", trait,
              " (", total_rows, " source rows)")
    }
  }

  if (!length(chunk_files)) stop("No rsID-mapped rows were produced for ", trait)
  output_con <- gzfile(out_path, open = "wt")
  on.exit(try(close(output_con), silent = TRUE), add = TRUE)
  header_written <- FALSE
  for (chunk_path in chunk_files) {
    chunk_con <- gzfile(chunk_path, open = "rt")
    chunk_lines <- readLines(chunk_con, warn = FALSE)
    close(chunk_con)
    if (length(chunk_lines)) {
      if (header_written) chunk_lines <- chunk_lines[-1L]
      if (length(chunk_lines)) {
        writeLines(chunk_lines, output_con)
        header_written <- TRUE
      }
    }
    rm(chunk_lines)
  }
  close(output_con)
  output_headers <- names(fread(out_path, nrows = 0L, showProgress = FALSE))
  if (!all(c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Z") %in% output_headers)) {
    stop("Batched output is missing required LDSC fields for ", trait)
  }
  message("Saved batched GRCh37 coordinate-mapped output for ", trait,
          " (", total_rows, " source rows; batch cache: ", cache_dir, ")")
  invisible(out_path)
}

# ---------------------------------------------------------------------------
# Content-level validation of existing cleaned outputs
# ---------------------------------------------------------------------------

# v1 outputs carry the correct headers but can contain NA alleles, pseudo-ID
# SNP strings, "chr"-prefixed CHR, out-of-range FRQ, or inconsistent Z values.
# An existing file is only reused when its content passes these checks;
# otherwise it is reprocessed in full.
validate_existing_output <- function(path, trait) {
  required <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Z", "FRQ")
  headers <- names(fread(path, nrows = 0L, showProgress = FALSE))
  dt <- fread(path, select = intersect(required, headers), showProgress = FALSE)
  problems <- character()

  # Missing rsIDs are acceptable (munging drops them); pseudo-ID strings such
  # as chr:pos_alleles are not -- they mean the file predates v2 cleaning.
  snp_nonna <- dt$SNP[!is.na(dt$SNP)]
  pct_pseudo <- if (length(snp_nonna)) {
    100 * mean(!grepl("^rs[0-9]+$", snp_nonna))
  } else 100
  if (pct_pseudo > 1) {
    problems <- c(problems, sprintf("%.1f%% of SNP IDs are pseudo-IDs, not rsIDs", pct_pseudo))
  }

  allele_ok <- (!is.na(dt$A1) & grepl("^[ACGT]+$", dt$A1)) &
    (!is.na(dt$A2) & grepl("^[ACGT]+$", dt$A2))
  pct_allele_ok <- 100 * mean(allele_ok)
  if (pct_allele_ok < 99) {
    problems <- c(problems, sprintf("only %.1f%% of rows have usable A1/A2 alleles", pct_allele_ok))
  }

  chr_ok <- !grepl("^chr", dt$CHR, ignore.case = TRUE) &
    toupper(dt$CHR) %in% c(as.character(1:22), "X", "Y", "MT")
  pct_chr_ok <- 100 * mean(chr_ok, na.rm = TRUE)
  if (pct_chr_ok < 99) {
    problems <- c(problems, sprintf("only %.1f%% of rows have valid CHR values", pct_chr_ok))
  }

  if ("FRQ" %in% names(dt)) {
    n_bad_frq <- sum(!is.na(dt$FRQ) & (dt$FRQ <= 0 | dt$FRQ >= 1))
    if (n_bad_frq > 0) {
      problems <- c(problems, sprintf("%d row(s) with FRQ outside (0, 1)", n_bad_frq))
    }
  }

  if (all(c("Z", "BETA", "SE") %in% names(dt))) {
    check <- is.finite(dt$Z) & is.finite(dt$BETA) & is.finite(dt$SE) & dt$SE > 0
    if (sum(check) > 0) {
      z_ref <- dt$BETA[check] / dt$SE[check]
      pct_inconsistent <- 100 * mean(abs(dt$Z[check] - z_ref) > 0.05 * pmax(abs(z_ref), 1))
      if (pct_inconsistent > 5) {
        problems <- c(problems, sprintf("%.1f%% of rows have Z inconsistent with BETA/SE",
                                        pct_inconsistent))
      }
    }
  }

  n_dup <- length(dt$SNP) - uniqueN(dt$SNP)
  if (n_dup > 0.01 * length(dt$SNP)) {
    problems <- c(problems, sprintf("%d duplicate SNP IDs", n_dup))
  }

  if (length(problems)) {
    message("Existing output for ", trait, " failed validation: ",
            paste(problems, collapse = "; "))
    rm(dt); gc(FALSE)
    invisible(FALSE)
  } else {
    rm(dt); gc(FALSE)
    invisible(TRUE)
  }
}

# ---------------------------------------------------------------------------
# Per-trait cleaning (one call per phenotype; safe to run in parallel workers)
# ---------------------------------------------------------------------------

clean_one_trait <- function(i) {
  trait <- manifest$phenotype[[i]]
  out_path <- file.path(out_dir, sprintf("%03d_%s_clean.tsv.gz", i, safe_name(trait)))
  result <- data.table(
    phenotype = trait, source = manifest$source[[i]], input_path = manifest$input_path[[i]],
    clean_path = "", status = "", note = "", source_headers = ""
  )

  # Skip 1: output already passed deep QC.
  if (trait %in% qc_passed) {
    result[, status := "skipped_qc_passed"]
    message("QC-passed; skipping: ", trait)
    return(result)
  }
  # Skip 2: already LDSC-munged downstream.
  if (trait %in% already_munged) {
    result[, status := "skipped_already_munged"]
    message("Already LDSC-munged; skipping: ", trait)
    return(result)
  }
  # Skip 3: a complete cleaned output already exists (resumable reruns).
  # Headers alone are not enough -- v1 outputs have the right headers but can
  # be content-broken -- so the file must also pass content validation.
  output_gzip_ok <- file.exists(out_path) && !is.na(file.info(out_path)$size) &&
    file.info(out_path)$size > 0 &&
    identical(as.integer(system2("gzip", c("-t", out_path), stdout = FALSE, stderr = FALSE)), 0L)
  if (output_gzip_ok) {
    output_headers <- tryCatch(
      names(fread(out_path, nrows = 0L, showProgress = FALSE)),
      error = function(e) character()
    )
    if (all(c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Z") %in% output_headers)) {
      if (validate_existing_output(out_path, trait)) {
        result[, status := "reused_existing"]
        result[, clean_path := normalizePath(out_path, mustWork = TRUE)]
        message("Existing output passed validation; skipping: ", trait)
        return(result)
      }
      message("Existing output is content-broken; reprocessing: ", trait)
    } else {
      message("Existing output is empty or incomplete; recreating: ", trait)
    }
  }

  tryCatch({
    input_path <- resolve_path(manifest$input_path[[i]], trait)
    message("Input for ", trait, ": ", input_path)
    if (!file.exists(input_path)) stop("Input file not found: ", input_path)

    skip <- header_skip(input_path)
    file_headers <- names(fread(input_path, nrows = 0L, skip = skip, showProgress = FALSE))
    result[, source_headers := paste(file_headers, collapse = " | ")]

    # If a PanUKB file includes ancestry-specific fields, use its EUR columns.
    is_panukb <- grepl("pan[ _-]?ukb",
                       paste(manifest$source[[i]], manifest$cohort[[i]]), ignore.case = TRUE)
    eur_headers <- grep("(_EUR|\\.EUR|^EUR_)", file_headers, value = TRUE, ignore.case = TRUE)
    recognized <- norm(file_headers) %in% norm(unlist(aliases))
    read_headers <- unique(c(file_headers[recognized], eur_headers))
    if (!length(read_headers)) stop("No recognized sumstats columns in ", input_path)

    # Large GRCh37 GWAS Catalog files that need coordinate-based rsID mapping
    # (no rsID column, or a column of chr:pos_alleles pseudo-IDs) are streamed
    # through rsID lookup in batches so they never load fully; feeding 100M+
    # rows to format_sumstats in memory would exhaust the node.
    has_snp_column <- length(find_col_matches(file_headers, aliases$SNP)) > 0L
    snp_is_pseudo <- FALSE
    if (has_snp_column) {
      snp_col_name <- find_col_matches(file_headers, aliases$SNP)[[1L]]
      peek <- fread(input_path, nrows = 1000, select = snp_col_name, skip = skip,
                    showProgress = FALSE)
      snp_is_pseudo <- mean(grepl("^rs[0-9]+$", peek[[1L]]), na.rm = TRUE) < 0.5
      rm(peek)
    }
    file_size_gb <- if (file.exists(input_path)) file.info(input_path)$size / 1024^3 else 0
    use_batched <- grepl("GWAS Catalog", manifest$source[[i]], ignore.case = TRUE) &&
      (!has_snp_column || snp_is_pseudo) &&
      (identical(trait, "Lymphocyte Count") || file_size_gb > batch_threshold_gb)
    if (use_batched) {
      source_has_n <- !is.na(first_match(file_headers, aliases$N))
      sheet_n <- if (is.na(n_col)) NA_real_ else as_num(ref[[n_col]][[i]])
      if (!source_has_n && (is.null(sheet_n) || !length(sheet_n) || is.na(sheet_n[[1L]]))) {
        stop("Lymphocyte Count needs a sample size in the source file or reference sheet before LDSC munging.")
      }
      clean_coordinates_in_batches(
        path = input_path, skip = skip, file_headers = file_headers,
        read_headers = read_headers, trait = trait, out_path = out_path,
        sheet_n = sheet_n
      )
      result[, status := "cleaned_batched"]
      result[, clean_path := normalizePath(out_path, mustWork = TRUE)]
      message("Saved cleaned sumstats: ", trait)
      return(result)
    }

    dt <- fread(input_path, select = read_headers, skip = skip, showProgress = TRUE)

    get_value <- function(key) {
      column <- best_column(names(dt), dt, aliases[[key]])
      if (is_panukb && length(eur_headers)) {
        eur_names <- c(paste0(aliases[[key]], "_EUR"), paste0(aliases[[key]], ".EUR"),
                       paste0("EUR_", aliases[[key]]))
        eur_column <- best_column(names(dt), dt, eur_names)
        if (!is.null(eur_column)) column <- eur_column
      }
      if (is.null(column)) NULL else dt[[column]]
    }

    a1 <- get_value("A1")
    a2 <- get_value("A2")
    # FinnGen uses ALT as the effect allele and REF as the other allele.
    if (all(c("ALT", "REF") %in% names(dt)) && is.na(first_match(names(dt), c("A1", "effect_allele", "Allele1")))) {
      a1 <- dt[[find_col_matches(names(dt), "ALT")[[1L]]]]
      a2 <- dt[[find_col_matches(names(dt), "REF")[[1L]]]]
    }

    beta <- as_num(get_value("BETA"))
    odds_ratio <- as_num(get_value("OR"))
    if (is.null(beta) && !is.null(odds_ratio)) beta <- log(odds_ratio)
    se <- as_num(get_value("SE"))
    p <- as_num(get_value("P"))
    logp <- as_num(get_value("LOGP"))
    z <- as_num(get_value("Z"))
    direction <- get_value("DIRECTION")
    if (is.null(p) && !is.null(logp)) p <- 10^(-logp)
    n_stats <- nrow(dt)
    if (is.null(p)) p <- rep(NA_real_, n_stats)
    if (is.null(z)) z <- rep(NA_real_, n_stats)
    if (is.null(beta)) beta <- rep(NA_real_, n_stats)
    if (is.null(se)) se <- rep(NA_real_, n_stats)

    # Fill Z from available signed statistics. Some files include a Z column
    # that is present but empty, so fill missing entries rather than checking
    # only whether the column exists.
    missing_z <- !is.finite(z)
    z_from_beta_se <- missing_z & is.finite(beta) & is.finite(se) & se > 0
    z[z_from_beta_se] <- beta[z_from_beta_se] / se[z_from_beta_se]
    missing_z <- !is.finite(z)
    z_from_beta_p <- missing_z & is.finite(beta) & is.finite(p) & p > 0 & p <= 1
    z[z_from_beta_p] <- sign(beta[z_from_beta_p]) *
      qnorm(pmax(p[z_from_beta_p] / 2, .Machine$double.xmin), lower.tail = FALSE)

    # A direction string is usable only when all its non-missing signs agree.
    # Mixed strings (for example '++-') are ambiguous and remain missing.
    if (!is.null(direction)) {
      direction <- trimws(as.character(direction))
      direction_sign <- rep(NA_real_, length(direction))
      direction_sign[grepl("^\\++$", direction)] <- 1
      direction_sign[grepl("^-+$", direction)] <- -1
      missing_z <- !is.finite(z)
      z_from_direction <- missing_z & is.finite(p) & p > 0 & p <= 1 &
        !is.na(direction_sign)
      z[z_from_direction] <- direction_sign[z_from_direction] *
        qnorm(pmax(p[z_from_direction] / 2, .Machine$double.xmin), lower.tail = FALSE)
    }
    if (is.null(p) && !is.null(z)) p <- 2 * pnorm(-abs(z))

    add_sheet_sample_size <- function(values, column) {
      if (is.na(column)) return(values)
      sheet_n <- as_num(ref[[column]][[i]])
      if (is.null(values)) return(rep(sheet_n, nrow(dt)))
      values[!is.finite(values)] <- sheet_n
      values
    }
    sample_n <- add_sheet_sample_size(as_num(get_value("N")), n_col)
    n_case <- add_sheet_sample_size(as_num(get_value("N_CASE")), n_case_col)
    n_control <- add_sheet_sample_size(as_num(get_value("N_CONTROL")), n_control_col)

    if (!is.null(p)) {
      p[!is.na(p) & p <= 0] <- .Machine$double.xmin
      p[!is.na(p) & p > 1] <- NA_real_
    }
    if (!is.null(z)) z[!is.finite(z)] <- NA_real_
    n_rows <- nrow(dt)
    out <- as.data.table(setNames(
      list(fill_value(get_value("SNP"), n_rows), fill_value(get_value("CHR"), n_rows),
           fill_value(get_value("BP"), n_rows, "numeric"), fill_value(a1, n_rows),
           fill_value(a2, n_rows), fill_value(beta, n_rows, "numeric"),
           fill_value(se, n_rows, "numeric"), fill_value(p, n_rows, "numeric"),
           fill_value(z, n_rows, "numeric")),
      c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Z")
    ))
    optional <- list(FRQ = get_value("FRQ"), INFO = get_value("INFO"), N = sample_n,
                     N_CASE = n_case, N_CONTROL = n_control)
    for (column in names(optional)) {
      value <- fill_value(optional[[column]], n_rows, "numeric")
      if (any(is.finite(value))) set(out, j = column, value = value)
    }

    # Normalize CHR/alleles, parse pseudo-IDs, drop unusable rows.
    out <- clean_variants(out, trait)

    # An allele column that matched but is entirely empty means the alias list
    # did not fit this source; fail loudly instead of writing a broken file.
    if (nrow(out) && (all(is.na(out$A1)) || all(is.na(out$A2)))) {
      stop("Effect or other allele column is entirely missing for ", trait,
           ". Matched headers were: ", paste(names(dt), collapse = ", "),
           ". Check sumstats_headers.tsv for the actual source headers and add ",
           "the correct alias to the alias list.")
    }
    if (!nrow(out)) {
      stop("No usable rows remain after cleaning for ", trait,
           ". Check the allele, FRQ, and coordinate columns in sumstats_headers.tsv.")
    }

    # Map coordinates to rsIDs whenever the SNP column is mostly not rsIDs
    # (missing, chr:pos_alleles pseudo-IDs, or heavily malformed).
    src_text <- paste(trait, manifest$source[[i]], manifest$cohort[[i]])
    pct_rs <- mean(grepl("^rs[0-9]+$", out$SNP), na.rm = TRUE)
    has_coordinates <- any(!is.na(out$CHR) & nzchar(out$CHR) & is.finite(out$BP) & out$BP > 0)
    if (pct_rs < 0.5) {
      if (!has_coordinates) {
        stop("Fewer than 50% of SNP IDs are rsIDs for ", trait,
             " and no usable CHR/BP coordinates are available for lookup. ",
             "Parsed headers: ", paste(names(dt), collapse = ", "), ".")
      }
      out <- map_coordinates_to_rsids(out, trait, src_text)
    }

    # Drop any residual rows that never received a usable rsID.
    bad_id <- is.na(out$SNP) | !grepl("^rs[0-9]+$", out$SNP)
    if (any(bad_id)) {
      message("Dropping ", sum(bad_id), " row(s) without usable rsIDs for ", trait)
      out <- out[!bad_id]
    }

    if (all(is.na(out$P)) || all(is.na(out$Z))) {
      n_valid_p <- sum(is.finite(out$P) & out$P > 0 & out$P <= 1)
      n_valid_z <- sum(is.finite(out$Z))
      direction_values <- if (is.null(direction)) "<no direction column>" else
        paste(utils::head(unique(na.omit(as.character(direction))), 8L), collapse = ", ")
      stop("Need P and signed effect information for ", trait,
           ". Usable P rows: ", n_valid_p, "; signed Z rows: ", n_valid_z,
           "; direction examples: ", direction_values,
           ". Check the actual source headers and direction encoding in sumstats_headers.tsv.")
    }
    if (all(is.na(out$N)) && (all(is.na(out$N_CASE)) || all(is.na(out$N_CONTROL)))) {
      stop("No sample-size information found for ", trait, ". Add N to the source or reference sheet.")
    }

    # Sample-size sanity: N should agree with N_CASE + N_CONTROL when both
    # exist -- unless N is an effective-sample-size column (e.g., DIAGRAM ships
    # Neff), which is deliberately smaller than N_CASE + N_CONTROL and is
    # exactly what LDSC wants.
    if (all(c("N", "N_CASE", "N_CONTROL") %in% names(out))) {
      n_med <- median(out$N, na.rm = TRUE)
      cc_med <- median(out$N_CASE + out$N_CONTROL, na.rm = TRUE)
      neff_ok <- out$N_CASE > 0 & out$N_CONTROL > 0
      neff_med <- if (any(neff_ok)) {
        median(4 / (1 / out$N_CASE[neff_ok] + 1 / out$N_CONTROL[neff_ok]), na.rm = TRUE)
      } else NA_real_
      is_neff_column <- is.finite(neff_med) && neff_med > 0 &&
        abs(n_med - neff_med) <= 0.25 * neff_med
      if (is.finite(n_med) && is.finite(cc_med) && cc_med > 0 &&
          abs(n_med - cc_med) > 0.25 * cc_med && !is_neff_column) {
        message("WARNING: ", trait, ": median N (", round(n_med),
                ") differs from median N_CASE + N_CONTROL (", round(cc_med),
                ") by more than 25%. Reconcile the sample-size columns before munging.")
      }
    }

    fwrite(out, out_path, sep = "\t", na = "NA")
    result[, status := "cleaned"]
    result[, clean_path := normalizePath(out_path, mustWork = TRUE)]
    message("Saved cleaned sumstats: ", trait)
    rm(dt, out)
    gc(FALSE)
  }, error = function(e) {
    result[, status := "failed"]
    result[, note := conditionMessage(e)]
    message("FAILED: ", trait, " -- ", conditionMessage(e))
  })
  result
}

# ---------------------------------------------------------------------------
# Run across phenotypes (parallel when CLEAN_CORES > 1) and write reports
# ---------------------------------------------------------------------------

if (n_cores == 1L) {
  results <- rbindlist(lapply(seq_len(nrow(manifest)), clean_one_trait), fill = TRUE)
} else {
  # One forked worker per trait; preschedule off so each worker holds only its
  # own input file in memory.
  results <- rbindlist(
    parallel::mclapply(seq_len(nrow(manifest)), clean_one_trait,
                       mc.cores = n_cores, mc.preschedule = FALSE),
    fill = TRUE
  )
}

fwrite(results, file.path(out_dir, "clean_status.tsv"), sep = "\t", na = "")
fwrite(results[status %in% c("cleaned", "cleaned_batched", "reused_existing") &
                 nzchar(clean_path),
               .(phenotype, source, input_path, clean_path)],
       file.path(out_dir, "clean_manifest.tsv"), sep = "\t", na = "")
fwrite(results[nzchar(source_headers),
               .(phenotype, source, file_path = input_path, columns = source_headers)],
       file.path(out_dir, "sumstats_headers.tsv"), sep = "\t", na = "")
failures <- results[status == "failed"]
if (nrow(failures)) {
  fwrite(failures, file.path(out_dir, "clean_failures.tsv"), sep = "\t", na = "")
  message(nrow(failures), " trait(s) failed; see ", file.path(out_dir, "clean_failures.tsv"),
          ". Failed traits are retried automatically on the next run.")
}
message("Done. Status: ", paste(results[, .N, by = status][, paste0(status, "=", N)],
                                collapse = ", "), ". Reports are in ", out_dir)

