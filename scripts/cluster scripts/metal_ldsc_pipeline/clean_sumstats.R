#!/usr/bin/env Rscript

# Read the selected rows from the reference sheet, standardize each file, and
# write a simple manifest for the LDSC step that follows.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: clean_sumstats.R <mapping CSV or URL> <output directory> [already_munged.csv] [analysis]")
}

mapping_path <- args[[1L]]
out_dir <- args[[2L]]
already_munged_path <- if (length(args) >= 3L) args[[3L]] else ""
analysis <- if (length(args) >= 4L) args[[4L]] else ""
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Primary cleaner version: resumable GRCh37 batches, signed-Z fallback, and partial-coordinate rsID mapping")

# Small helpers keep source-format differences in one place.
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
find_col <- function(headers, aliases) {
  hit <- match(norm(aliases), norm(headers), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) headers[[hit[[1L]]]] else NULL
}

# The reference sheet uses exact phenotype names to distinguish cohorts.
ref <- fread(mapping_path, na.strings = c("", "NA", "N/A"), showProgress = FALSE)
ph_col <- find_col(names(ref), "phenotype")
path_col <- find_col(names(ref), "file location")
if (is.null(ph_col) || is.null(path_col)) {
  stop("The reference sheet must have 'phenotype' and 'file location' columns.")
}
source_col <- find_col(names(ref), "source")
cohort_col <- find_col(names(ref), "cohort")
analysis_col <- find_col(names(ref), "analysis")
n_col <- find_col(names(ref), c("n", "sample size", "total n"))
n_case_col <- find_col(names(ref), c("n case", "n cases", "n_cas"))
n_control_col <- find_col(names(ref), c("n control", "n controls", "n_con"))

if (nzchar(analysis)) {
  if (is.null(analysis_col)) stop("The reference sheet needs an 'analysis' column for this run.")
  analysis_values <- toupper(trimws(as.character(ref[[analysis_col]])))
  requested_analysis <- toupper(trimws(analysis))
  keep_rows <- which(!is.na(analysis_values) & analysis_values == requested_analysis)
  ref <- ref[keep_rows]
  message("Using ", nrow(ref), " row(s) marked analysis = ", analysis)
}

# Save the requested phenotype list even when some files are already munged.
manifest <- data.table(
  phenotype = as.character(ref[[ph_col]]),
  analysis = if (!is.null(analysis_col)) as.character(ref[[analysis_col]]) else "",
  source = if (!is.null(source_col)) as.character(ref[[source_col]]) else "",
  cohort = if (!is.null(cohort_col)) as.character(ref[[cohort_col]]) else "",
  input_path = trimws(as.character(ref[[path_col]]))
)
manifest <- manifest[!is.na(phenotype) & nzchar(phenotype) & !is.na(input_path) & nzchar(input_path)]
if (nzchar(analysis) &&
    any(toupper(trimws(manifest$analysis)) != requested_analysis, na.rm = TRUE)) {
  stop("Internal filter error: non-", analysis, " rows entered the primary manifest.")
}
fwrite(manifest, file.path(out_dir, "primary_manifest.tsv"), sep = "\t", na = "")
message(
  "Primary manifest written: ", file.path(out_dir, "primary_manifest.tsv"),
  " (", nrow(manifest), " row(s); analysis = ", analysis, ")"
)

# Exact phenotype matches in this registry already have LDSC-ready files.
already_munged <- character()
if (nzchar(already_munged_path)) {
  previous <- fread(already_munged_path, na.strings = c("", "NA", "N/A"), showProgress = FALSE)
  previous_ph_col <- find_col(names(previous), "phenotype")
  already_munged <- as.character(previous[[previous_ph_col]])
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
    if (!identical(as.integer(status), 0L) || !file.exists(temporary) || file.info(temporary)$size == 0) {
      unlink(temporary)
      stop("Download failed: ", path)
    }
    if (!file.rename(temporary, target)) stop("Could not cache download: ", target)
  }
  target
}

# PGC files can begin with ## metadata. Read past it without decoding arbitrary
# leading bytes as UTF-8.
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

# Common column names used by the included GWAS sources.
aliases <- list(
  SNP = c("rsids", "rsid", "rs_id", "SNP", "SNPID", "MarkerName", "variant_id", "variant", "ID"),
  CHR = c("#chrom", "CHR", "chromosome", "chromsome", "chrom", "chr"),
  BP = c("POS", "BP", "base_pair_location", "position", "pos"),
  A1 = c("effect_allele", "A1", "Allele1", "ALT", "inc_allele", "EA"),
  A2 = c("other_allele", "A2", "Allele2", "REF", "non_effect_allele", "OA"),
  BETA = c("beta", "BETA", "Effect", "B", "log_odds"),
  OR = c("OR", "odds_ratio"),
  SE = c("standard_error", "SE", "sebeta", "StdErr", "stderr"),
  P = c("p_value", "P", "pval", "P.value", "Pvalue"),
  LOGP = c("neg_log_10_p_value", "mlogp", "minus_log10_p", "neg_log10_p", "LOG10P"),
  Z = c("Z", "Zscore", "Z_score", "z_score", "tstat"),
  DIRECTION = c("direction", "effect_direction", "beta_direction"),
  FRQ = c("effect_allele_frequency", "EAF", "af_alt", "AF_Allele2", "FRQ", "af"),
  INFO = c("INFO", "imputation_info"),
  N = c("N", "n", "sample_size", "total_n", "nCompleteSamples", "Neff"),
  N_CASE = c("N_CASE", "N_CASES", "N_CAS", "Nca", "cases_n"),
  N_CONTROL = c("N_CONTROL", "N_CONTROLS", "N_CON", "Nco", "controls_n")
)

# Stream the very large GRCh37 GWAS Catalog lymphocyte file through rsID
# lookup in batches, caching finished batches so an interrupted job can resume.
coordinate_batch_rows <- suppressWarnings(as.integer(Sys.getenv("CLEAN_CHUNK_ROWS", "250000")))
if (is.na(coordinate_batch_rows) || coordinate_batch_rows < 1L) {
  stop("CLEAN_CHUNK_ROWS must be a positive integer.")
}

clean_gwas_catalog_coordinates_in_batches <- function(path, skip, file_headers,
                                                       read_headers, trait, out_path,
                                                       sheet_n) {
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
  # Query the installed GRCh37 dbSNP index directly for each batch. Calling
  # format_sumstats() repeatedly reloads/checks the reference and accumulated
  # enough memory to OOM this 50-million-row source.
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

  # Look for caches made by the earlier signature versions before starting a
  # new folder. The rest of the signature must still match exactly, including
  # input metadata, selected columns, sample size, and rows per batch.
  signature_versions <- 2:3
  cache_keys <- vapply(signature_versions, function(version) {
    signature <- c(list(version = as.integer(version)), signature_fields)
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
                           signature_versions[reusable] == 3L, decreasing = TRUE)][[1L]]
    cache_key <- cache_keys[[best]]
    message("Reusing compatible cache folder for ", trait, ": ",
            candidate_dirs[[best]], " (", completed_counts[[best]], " completed batches)")
  } else {
    cache_key <- cache_keys[[which(signature_versions == 3L)]]
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
        selected <- find_col(names(dt), aliases[[key]])
        if (is.null(selected)) NULL else dt[[selected]]
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
        SNP = rep(NA_character_, n),
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
      out <- out[!is.na(CHR) & nzchar(CHR) & is.finite(BP) & BP > 0 &
                   !is.na(A1) & nzchar(A1) & !is.na(A2) & nzchar(A2) &
                   is.finite(BETA) & is.finite(SE) & SE > 0 &
                   is.finite(P) & P > 0 & P <= 1]
      missing_id <- is.na(out$SNP) | !grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE)
      rm(dt, lines)

      # Match only the distinct chromosome/position pairs in this chunk. The
      # dbSNP package returns rsIDs in RefSNP_id; positions with multiple
      # distinct rsIDs are treated as ambiguous and left unmapped.
      out[, lookup_chr := toupper(sub("^chr", "", trimws(as.character(CHR)), ignore.case = TRUE))]
      out[lookup_chr == "M", lookup_chr := "MT"]
      out[lookup_chr == "23", lookup_chr := "X"]
      out[lookup_chr == "24", lookup_chr := "Y"]
      out[lookup_chr == "25", lookup_chr := "MT"]
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
  # close() below finalizes the gzip stream. Use a guarded cleanup for errors;
  # isOpen() can itself fail on an already-closed connection on some R builds.
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

cleaned <- list()
headers_report <- list()
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$phenotype[[i]]
  out_path <- file.path(out_dir, sprintf("%03d_%s_clean.tsv.gz", i, safe_name(trait)))

  # Reuse a complete-looking cleaned file from an earlier run. This lets a
  # rerun resume after a failed input without reprocessing finished traits.
  output_gzip_ok <- file.exists(out_path) && !is.na(file.info(out_path)$size) &&
    file.info(out_path)$size > 0 &&
    identical(as.integer(system2("gzip", c("-t", out_path), stdout = FALSE, stderr = FALSE)), 0L)
  if (output_gzip_ok) {
    output_headers <- tryCatch(
      names(fread(out_path, nrows = 0L, showProgress = FALSE)),
      error = function(e) character()
    )
    if (all(c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Z") %in% output_headers)) {
      cleaned[[length(cleaned) + 1L]] <- data.table(
        phenotype = trait, source = manifest$source[[i]], input_path = manifest$input_path[[i]],
        clean_path = normalizePath(out_path, mustWork = TRUE)
      )
      message("Completed output already exists; skipping: ", trait)
      next
    }
    message("Existing output is empty or incomplete; recreating: ", trait)
  }

  input_path <- resolve_path(manifest$input_path[[i]], trait)
  message("Input for ", trait, ": ", input_path)

  if (trait %in% already_munged) {
    message("Already LDSC-munged; leaving as-is: ", trait)
    next
  }
  if (!file.exists(input_path)) stop("Input file not found for ", trait, ": ", input_path)

  skip <- header_skip(input_path)
  file_headers <- names(fread(input_path, nrows = 0L, skip = skip, showProgress = FALSE))
  headers_report[[length(headers_report) + 1L]] <- data.table(
    phenotype = trait, source = manifest$source[[i]], file_path = manifest$input_path[[i]],
    columns = paste(file_headers, collapse = " | ")
  )

  # If a PanUKB file includes ancestry-specific fields, use its EUR columns.
  is_panukb <- grepl("pan[ _-]?ukb", paste(manifest$source[[i]], manifest$cohort[[i]]), ignore.case = TRUE)
  eur_headers <- grep("(_EUR|\\.EUR|^EUR_)", file_headers, value = TRUE, ignore.case = TRUE)
  recognized <- norm(file_headers) %in% norm(unlist(aliases))
  read_headers <- unique(c(file_headers[recognized], eur_headers))
  if (!length(read_headers)) stop("No recognized sumstats columns in ", input_path)

  # This GWAS Catalog file has coordinates but no rsID column. Map its known
  # GRCh37 coordinates in batches instead of loading all 50 million rows.
  is_chunked_lymphocyte <- identical(trait, "Lymphocyte Count") &&
    grepl("GWAS Catalog", manifest$source[[i]], ignore.case = TRUE) &&
    is.null(find_col(file_headers, aliases$SNP))
  if (is_chunked_lymphocyte) {
    source_has_n <- !is.null(find_col(file_headers, aliases$N))
    sheet_n <- if (is.null(n_col)) NA_real_ else as_num(ref[[n_col]][[i]])
    if (!source_has_n && (is.null(sheet_n) || !length(sheet_n) || is.na(sheet_n[[1L]]))) {
      stop("Lymphocyte Count needs a sample size in the source file or reference sheet before LDSC munging.")
    }
    clean_gwas_catalog_coordinates_in_batches(
      path = input_path, skip = skip, file_headers = file_headers,
      read_headers = read_headers, trait = trait, out_path = out_path,
      sheet_n = sheet_n
    )
    cleaned[[length(cleaned) + 1L]] <- data.table(
      phenotype = trait, source = manifest$source[[i]], input_path = manifest$input_path[[i]],
      clean_path = normalizePath(out_path, mustWork = TRUE)
    )
    next
  }

  dt <- fread(input_path, select = read_headers, skip = skip, showProgress = TRUE)

  get_value <- function(key) {
    options <- aliases[[key]]
    if (is_panukb && length(eur_headers)) {
      eur_names <- c(paste0(options, "_EUR"), paste0(options, ".EUR"), paste0("EUR_", options))
      selected <- find_col(names(dt), eur_names)
      if (!is.null(selected)) return(dt[[selected]])
    }
    selected <- find_col(names(dt), options)
    if (is.null(selected)) NULL else dt[[selected]]
  }

  a1 <- get_value("A1")
  a2 <- get_value("A2")
  # FinnGen uses ALT as the effect allele and REF as the other allele.
  if (all(c("ALT", "REF") %in% names(dt)) &&
      is.null(find_col(names(dt), c("A1", "effect_allele", "Allele1")))) {
    a1 <- dt[[find_col(names(dt), "ALT")]]
    a2 <- dt[[find_col(names(dt), "REF")]]
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
    if (is.null(column)) return(values)
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
  snp <- fill_value(get_value("SNP"), n_rows)
  chr <- fill_value(get_value("CHR"), n_rows)
  bp <- fill_value(get_value("BP"), n_rows, "numeric")
  values <- list(snp, chr, bp, fill_value(a1, n_rows), fill_value(a2, n_rows),
                 fill_value(beta, n_rows, "numeric"), fill_value(se, n_rows, "numeric"),
                 fill_value(p, n_rows, "numeric"), fill_value(z, n_rows, "numeric"))
  out <- as.data.table(setNames(
    values,
    c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "Z")
  ))
  optional <- list(FRQ = get_value("FRQ"), INFO = get_value("INFO"), N = sample_n,
                   N_CASE = n_case, N_CONTROL = n_control)
  for (column in names(optional)) {
    value <- fill_value(optional[[column]], n_rows, "numeric")
    if (any(is.finite(value))) set(out, j = column, value = value)
  }

  if (all(is.na(out$SNP)) && any(!is.na(out$CHR) & is.finite(out$BP) & out$BP > 0)) {
    is_finngen <- grepl(
      "finngen",
      paste(trait, manifest$source[[i]], manifest$cohort[[i]]),
      ignore.case = TRUE
    )
    is_lymphocyte_gwas_catalog <- identical(trait, "Lymphocyte Count") &&
      grepl("GWAS Catalog", manifest$source[[i]], ignore.case = TRUE)
    is_diagram <- grepl(
      "DIAGRAM",
      paste(manifest$source[[i]], manifest$cohort[[i]]),
      ignore.case = TRUE
    )
    coordinate_build <- if (is_finngen) {
      "GRCh38"
    } else if (is_lymphocyte_gwas_catalog || is_diagram) {
      "GRCh37"
    } else {
      stop("No variant ID column found for ", trait,
           ". Source/cohort were '", manifest$source[[i]], "' / '",
           manifest$cohort[[i]], "'. Coordinate-based rsID lookup needs a known genome build.")
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
      stop("No rsIDs found for ", trait, ". Missing package(s): ",
           paste(missing_reference_packages, collapse = ", "),
           ". Install them in R_LIBS_USER before running this job.")
    }
    for (package in required_reference_packages) {
      message(package, " version: ", as.character(utils::packageVersion(package)))
    }
    message("No rsID column for ", trait, "; looking up rsIDs from ", coordinate_build,
            " coordinates")
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
    out <- as.data.table(mapped)
    if (!"SNP" %in% names(out) || all(is.na(out$SNP))) {
      stop("Could not map coordinates to rsIDs for ", trait,
           ". Check that the file build matches ", coordinate_build,
           " and that its matching SNPlocs reference is installed.")
    }
  }
  if (all(is.na(out$SNP))) {
    id_column <- find_col(names(dt), aliases$SNP)
    chr_column <- find_col(names(dt), aliases$CHR)
    bp_column <- find_col(names(dt), aliases$BP)
    id_count <- if (is.null(id_column)) 0L else sum(!is.na(dt[[id_column]]) & nzchar(trimws(as.character(dt[[id_column]]))))
    coordinate_count <- if (is.null(chr_column) || is.null(bp_column)) 0L else
      sum(!is.na(dt[[chr_column]]) & nzchar(trimws(as.character(dt[[chr_column]]))) &
            is.finite(as_num(dt[[bp_column]])) & as_num(dt[[bp_column]]) > 0)
    stop("No usable variant IDs or coordinate-based rsID mapping for ", trait,
         ". Parsed ID column: ", if (is.null(id_column)) "<none>" else id_column,
         " (", id_count, " nonblank values); chromosome column: ",
         if (is.null(chr_column)) "<none>" else chr_column,
         "; position column: ", if (is.null(bp_column)) "<none>" else bp_column,
         " (", coordinate_count, " rows with both coordinates). Parsed headers: ",
         paste(names(dt), collapse = ", "),
         ". Check that these are the headers from the intended source file.")
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

  fwrite(out, out_path, sep = "\t", na = "NA")
  cleaned[[length(cleaned) + 1L]] <- data.table(
    phenotype = trait, source = manifest$source[[i]], input_path = manifest$input_path[[i]],
    clean_path = normalizePath(out_path, mustWork = TRUE)
  )
  message("Saved cleaned sumstats: ", trait)
  rm(dt, out)
  gc(FALSE)
}

if (!length(headers_report)) headers_report <- list(data.table(phenotype = character(), source = character(), file_path = character(), columns = character()))
if (!length(cleaned)) cleaned <- list(data.table(phenotype = character(), source = character(), input_path = character(), clean_path = character()))
fwrite(rbindlist(headers_report, fill = TRUE), file.path(out_dir, "sumstats_headers.tsv"), sep = "\t", na = "")
fwrite(rbindlist(cleaned, fill = TRUE), file.path(out_dir, "clean_manifest.tsv"), sep = "\t", na = "")
message("Done. Headers and cleaned-file manifest are in ", out_dir)
