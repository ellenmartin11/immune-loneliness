#!/usr/bin/env Rscript

############################
# METAL Preparation Script #
############################

# This script prepares sumstats for combination in METAL, and is implemented by submit_metal_ivw.sh

# Clean the spreadsheet's METAL rows, then combine the seven prespecified pairs with IVW.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: prepare_metal_ivw.R <mapping CSV or URL> <output directory> <METAL executable>")
}
mapping_path <- args[[1L]]
out_dir <- args[[2L]]
metal_bin <- args[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
input_dir <- file.path(out_dir, "clean_inputs")
download_dir <- file.path(out_dir, "raw_downloads")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
message("METAL prep version: chunked PanUKB liftover, resumable cache, direct chain download, and duplicate-rsID exclusion")

# Check both source and target reference resources before processing large files.
required_reference_packages <- c(
  "MungeSumstats",
  "SNPlocs.Hsapiens.dbSNP155.GRCh37",
  "BSgenome.Hsapiens.1000genomes.hs37d5",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38",
  "BSgenome.Hsapiens.NCBI.GRCh38"
)
missing_reference_packages <- required_reference_packages[
  !vapply(required_reference_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_reference_packages)) {
  stop("Missing MungeSumstats/reference package(s): ",
       paste(missing_reference_packages, collapse = ", "),
       ". Install them in R_LIBS_USER before submitting this job.")
}
for (package in required_reference_packages) {
  message(package, " version: ", as.character(utils::packageVersion(package)))
}

# Cleaning names
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
find_col <- function(nms, choices) {
  hit <- match(norm(choices), norm(nms), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) nms[[hit[[1L]]]] else NULL
}
safe_name <- function(x) gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", x))
as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))

# Keep only dbSNP rsIDs; never use chromosome-position IDs across genome builds.
normalize_rsid <- function(x) {
  x <- tolower(trimws(as.character(x)))
  # FinnGen rsids can contain more than one identifier; use the first valid rsID.
  match_pos <- regexpr("rs[0-9]+", x, ignore.case = TRUE, perl = TRUE)
  match_len <- attr(match_pos, "match.length")
  matched <- !is.na(match_pos) & match_pos > 0L
  hit <- rep(NA_character_, length(x))
  hit[matched] <- substring(x[matched], match_pos[matched],
                             match_pos[matched] + match_len[matched] - 1L)
  hit
}

# Download linked (non-downloaded) raw sumstats once; local file locations pass through unchanged.
resolve_path <- function(path, phenotype) {
  if (!grepl("^https?://", path, ignore.case = TRUE)) return(path)
  url_name <- basename(sub("[?#].*$", "", path))
  target <- file.path(download_dir, paste0(safe_name(phenotype), "__", url_name))
  if (file.exists(target) && file.info(target)$size > 0) return(target)
  
  # Caching so we don't redownload everytime we run the SLURM jobs
  tmp <- tempfile(tmpdir = download_dir)
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  if (!nzchar(Sys.which("curl"))) stop("curl is required to download: ", path)
  status <- system2("curl", c("--location", "--fail", "--silent", "--show-error",
                               "--retry", "3", "--retry-delay", "5", "--output", tmp, path))
  if (!identical(as.integer(status), 0L) || !file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Download failed: ", path)
  }
  if (grepl("\\.gz$", url_name, ignore.case = TRUE) &&
      !identical(as.integer(system2("gzip", c("-t", tmp), stdout = FALSE, stderr = FALSE)), 0L)) {
    stop("Downloaded file failed gzip integrity check: ", path)
  }
  if (!file.rename(tmp, target)) stop("Unable to cache downloaded input: ", target)
  target
}

ref <- fread(mapping_path, strip.white = FALSE, na.strings = c("", "NA", "N/A"), showProgress = FALSE)
ph_col <- find_col(names(ref), "phenotype")
analysis_col <- find_col(names(ref), "analysis")
source_col <- find_col(names(ref), "source")
cohort_col <- find_col(names(ref), "cohort")
path_col <- find_col(names(ref), "file location")
n_col <- find_col(names(ref), c("n", "sample size", "total n"))
if (any(vapply(list(ph_col, analysis_col, source_col, cohort_col, path_col, n_col), is.null, logical(1)))) {
  stop("Mapping sheet needs phenotype, analysis, cohort, source, n, and file location columns.")
}
# Should only run on my rows that specifically indicate are meant for METAL analysis (rather than Primary)
ref <- ref[as.character(ref[[analysis_col]]) == "METAL"]
if (!nrow(ref)) stop("No rows have analysis = METAL.")

# Each entry is an EXACT phenotype name from the sheet -- checked for spelling
pairs <- list(
  "Albumin METAL" = c("Albumin PanUKB", "Albumin FinnGen"),
  "C-reactive protein METAL" = c("C-reactive protein UKB", "C-reactive protein FinnGen"),
  "Eosinophil Count METAL" = c("Eosinophil Count PanUKB", "Eosinophil Count FinnGen"),
  "Monocytes METAL" = c("Monocytes PanUKB", "Monocytes FinnGen"),
  "Neutrophils METAL" = c("Neutrophils PanUKB", "Neutrophils FinnGen"),
  "Platelet Count METAL" = c("Platelet Count PanUKB", "Platelet Count FinnGen"),
  "Triglycerides METAL" = c("Triglycerides PanUKB", "Triglycerides FinnGen")
)
if (anyDuplicated(ref[[ph_col]])) stop("Duplicate exact phenotype labels among METAL rows.")
missing_traits <- setdiff(unique(unlist(pairs)), as.character(ref[[ph_col]]))
if (length(missing_traits)) stop("Missing METAL row(s) from the spreadsheet: ", paste(missing_traits, collapse = "; "))

# PanUKB files are streamed in batches so mapping and liftover do not hold the
# complete 29-million-row table and its intermediate copies in memory.
panukb_batch_rows <- suppressWarnings(as.integer(Sys.getenv("METAL_CHUNK_ROWS", "250000")))
if (is.na(panukb_batch_rows) || panukb_batch_rows < 1L) {
  stop("METAL_CHUNK_ROWS must be a positive integer.")
}
finngen_batch_rows <- suppressWarnings(as.integer(Sys.getenv("METAL_FINGEN_CHUNK_ROWS", "100000")))
if (is.na(finngen_batch_rows) || finngen_batch_rows < 1L) {
  stop("METAL_FINGEN_CHUNK_ROWS must be a positive integer.")
}

prepare_panukb_batch <- function(dt, n_value, phenotype, lookup_path, chain_file) {
  get_batch <- function(choices) {
    selected <- find_col(names(dt), choices)
    if (is.null(selected)) NULL else dt[[selected]]
  }
  n <- nrow(dt)
  id <- get_batch(c("rs_id", "rsids", "rsid", "dbsnp", "SNP", "MarkerName", "variant_id", "variant", "ID"))
  if (is.null(id)) id <- rep(NA_character_, n)
  id <- normalize_rsid(id)
  if (length(id) != n) stop("rsID normalization changed the batch row count for ", phenotype)
  chr <- get_batch(c("CHR", "chromosome", "chrom", "chr"))
  bp <- get_batch(c("BP", "POS", "base_pair_location", "pos", "position"))
  a1 <- get_batch(c("ALT", "effect_allele", "A1"))
  a2 <- get_batch(c("REF", "other_allele", "A2"))
  beta <- as_num(get_batch(c("beta_EUR", "BETA", "beta", "beta_meta", "Effect", "B", "log_odds")))
  se <- as_num(get_batch(c("se_EUR", "SE", "standard_error", "sebeta", "se_meta", "StdErr", "stderr")))
  p <- as_num(get_batch(c("pval_EUR", "P", "p_value", "pval", "pval_meta", "pval_metal", "P.value", "Pvalue")))
  if (is.null(chr) || is.null(bp) || is.null(a1) || is.null(a2) || is.null(beta) || is.null(se) || is.null(p)) {
    stop("Missing required PanUKB columns in batch for ", phenotype)
  }
  out <- data.table(SNP = id, CHR = as.character(chr), BP = as_num(bp),
                    A1 = as.character(a1), A2 = as.character(a2),
                    BETA = beta, SE = se, P = p, N = rep(n_value, n))
  if (nrow(out) != n) stop("PanUKB batch table row count changed for ", phenotype)
  out <- out[!is.na(A1) & nzchar(A1) & !is.na(A2) & nzchar(A2) &
               is.finite(BETA) & is.finite(SE) & SE > 0 &
               is.finite(P) & P > 0 & P <= 1 & !is.na(CHR) & nzchar(CHR) & is.finite(BP)]
  if (!nrow(out)) return(out[, .(SNP, A1, A2, BETA, SE, P, N)])

  missing_rsid <- is.na(out$SNP) | !grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE)
  if (any(missing_rsid)) {
    lookup_input <- data.table::copy(out[missing_rsid])
    lookup_input[, SNP := NULL]
    formatted <- MungeSumstats::format_sumstats(
      path = lookup_input, ref_genome = "GRCh37", dbSNP = 155,
      snp_ids_are_rs_ids = FALSE, return_data = TRUE,
      return_format = "data.table", save_path = lookup_path,
      force_new = TRUE, log_folder_ind = FALSE, sort_coordinates = FALSE
    )
    rm(lookup_input)
    if (!data.table::is.data.table(formatted) && !is.data.frame(formatted)) {
      stop("MungeSumstats did not return a table for PanUKB batch of ", phenotype)
    }
    mapped <- as.data.table(formatted)
    rm(formatted)
    required <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
    if (!all(required %in% names(mapped))) {
      stop("Unexpected coordinate-lookup columns for ", phenotype, ": ", paste(names(mapped), collapse = ", "))
    }
    valid_rsid <- !is.na(mapped$SNP) & grepl("^rs[0-9]+$", mapped$SNP, ignore.case = TRUE)
    mapped <- mapped[valid_rsid & !is.na(A1) & !is.na(A2) & is.finite(BETA) &
                       is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
    if (nrow(mapped)) mapped[, N := n_value]
    out <- data.table::rbindlist(list(out[!missing_rsid], mapped), fill = TRUE, use.names = TRUE)
    rm(mapped)
  }
  if (!nrow(out)) return(data.table(SNP = character(), A1 = character(), A2 = character(),
                                    BETA = numeric(), SE = numeric(), P = numeric(), N = numeric()))

  data.table::setcolorder(out, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"))
  lifted <- MungeSumstats::liftover(
    sumstats_dt = out, ref_genome = "GRCh37", convert_ref_genome = "GRCh38",
    chain_source = "ensembl", local_chain = chain_file, imputation_ind = FALSE
  )
  rm(out)
  lifted <- as.data.table(lifted)
  if (!all(c("SNP", "A1", "A2", "BETA", "SE", "P") %in% names(lifted))) {
    stop("Unexpected liftover output columns for ", phenotype)
  }
  if (!"N" %in% names(lifted)) lifted[, N := n_value]
  lifted <- lifted[!is.na(SNP) & grepl("^rs[0-9]+$", SNP, ignore.case = TRUE) &
                     !is.na(A1) & !is.na(A2) & is.finite(BETA) & is.finite(SE) & SE > 0 &
                     is.finite(P) & P > 0 & P <= 1]
  # A dbSNP rsID can label more than one allelic record. Exclude every
  # occurrence of such IDs rather than choosing an arbitrary effect estimate.
  duplicate_ids <- unique(lifted[duplicated(SNP) | duplicated(SNP, fromLast = TRUE), SNP])
  if (length(duplicate_ids)) {
    message("Excluding ", sum(lifted$SNP %in% duplicate_ids), " ambiguous rows across ",
            length(duplicate_ids), " duplicate rsIDs in PanUKB batch for ", phenotype,
            "; examples: ", paste(head(duplicate_ids, 5L), collapse = ", "))
    lifted <- lifted[!SNP %in% duplicate_ids]
  }
  lifted[, .(SNP, A1, A2, BETA, SE, P, N)]
}

prepare_panukb_in_batches <- function(path, header, select, n_value, phenotype, out_path) {
  info <- file.info(path)
  signature_fields <- list(path = normalizePath(path), size = info$size,
                           mtime = as.numeric(info$mtime), header = header, select = select,
                           n = n_value, batch_rows = panukb_batch_rows,
                           mungesumstats = as.character(utils::packageVersion("MungeSumstats")))
  signature_file <- tempfile()
  on.exit(unlink(signature_file), add = TRUE)
  phenotype_cache_dir <- file.path(out_dir, "batch_cache", safe_name(phenotype))
  dir.create(phenotype_cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Earlier script uploads used signature versions 2 and 3. Check those exact
  # cache keys before creating a new folder, so completed work survives script
  # replacement while still requiring the same input and processing settings.
  signature_versions <- 1:3
  cache_keys <- vapply(signature_versions, function(version) {
    writeBin(serialize(c(list(version = as.integer(version)), signature_fields), NULL),
             signature_file)
    unname(tools::md5sum(signature_file))
  }, character(1))
  candidate_dirs <- file.path(phenotype_cache_dir, cache_keys)
  candidate_exists <- dir.exists(candidate_dirs)
  completed_counts <- vapply(candidate_dirs, function(candidate) {
    if (!dir.exists(candidate)) return(0L)
    length(Sys.glob(file.path(candidate, "chunk_*.tsv.gz"))) +
      length(Sys.glob(file.path(candidate, "chunk_*.empty")))
  }, integer(1))
  reusable <- which(candidate_exists & completed_counts > 0L)
  if (length(reusable)) {
    # Prefer the current key on ties; otherwise continue from the cache with
    # the most completed batches among the matching signatures.
    best <- reusable[order(completed_counts[reusable],
                           signature_versions[reusable] == 2L, decreasing = TRUE)][[1L]]
    cache_key <- cache_keys[[best]]
    message("Reusing compatible cache folder for ", phenotype, ": ",
            candidate_dirs[[best]], " (", completed_counts[[best]], " completed batches)")
  } else {
    cache_key <- cache_keys[[which(signature_versions == 2L)]]
    message("No compatible completed cache found for ", phenotype,
            "; starting a new cache at ", file.path(phenotype_cache_dir, cache_key))
  }
  cache_dir <- file.path(phenotype_cache_dir, cache_key)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  chain_dir <- file.path(out_dir, "batch_cache", "reference")
  dir.create(chain_dir, recursive = TRUE, showWarnings = FALSE)
  chain_file <- file.path(chain_dir, "GRCh37_to_GRCh38.chain")
  compressed_chain <- paste0(chain_file, ".gz")
  chain_ok <- file.exists(chain_file) && file.info(chain_file)$size > 0 &&
    startsWith(readLines(chain_file, n = 1L, warn = FALSE), "chain ")
  if (!chain_ok) {
    if (file.exists(chain_file)) unlink(chain_file)
    compressed_ok <- file.exists(compressed_chain) && file.info(compressed_chain)$size > 0 &&
      identical(as.integer(system2("gzip", c("-t", compressed_chain), stdout = FALSE, stderr = FALSE)), 0L)
    if (!compressed_ok) {
      if (file.exists(compressed_chain)) unlink(compressed_chain)
      temporary_gzip <- tempfile(tmpdir = chain_dir, fileext = ".chain.gz")
      chain_url <- "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz"
      status <- system2("curl", c("--location", "--fail", "--silent", "--show-error",
                                   "--retry", "3", "--output", temporary_gzip, chain_url))
      compressed_ok <- identical(as.integer(status), 0L) && file.exists(temporary_gzip) &&
        file.info(temporary_gzip)$size > 0 &&
        identical(as.integer(system2("gzip", c("-t", temporary_gzip), stdout = FALSE, stderr = FALSE)), 0L)
      if (!compressed_ok) {
        unlink(temporary_gzip)
        stop("Could not download a valid GRCh37-to-GRCh38 Ensembl chain file. Check cluster network access and curl/gzip availability.")
      }
      if (!file.rename(temporary_gzip, compressed_chain)) {
        unlink(temporary_gzip)
        stop("Could not cache the compressed GRCh37-to-GRCh38 chain file in ", chain_dir)
      }
    }

    temporary_chain <- tempfile(tmpdir = chain_dir, fileext = ".chain")
    status <- system2("gzip", c("-dc", compressed_chain), stdout = temporary_chain)
    uncompressed_ok <- identical(as.integer(status), 0L) && file.exists(temporary_chain) &&
      file.info(temporary_chain)$size > 0 &&
      startsWith(readLines(temporary_chain, n = 1L, warn = FALSE), "chain ")
    if (!uncompressed_ok) {
      unlink(temporary_chain)
      stop("Could not create a valid uncompressed GRCh37-to-GRCh38 chain file from ", compressed_chain)
    }
    if (!file.rename(temporary_chain, chain_file)) {
      unlink(temporary_chain)
      stop("Could not cache the uncompressed GRCh37-to-GRCh38 chain file in ", chain_dir)
    }
    message("Downloaded/cached and unpacked GRCh37-to-GRCh38 chain file: ", chain_file)
  } else {
    message("Reusing cached uncompressed GRCh37-to-GRCh38 chain file: ", chain_file)
  }

  input_con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(input_con), add = TRUE)
  header_line <- readLines(input_con, n = 1L, warn = FALSE)
  if (!length(header_line)) stop("Could not read PanUKB header: ", path)
  input_header <- strsplit(sub("\\r$", "", header_line), "\t", fixed = TRUE)[[1L]]
  if (!identical(input_header, header)) stop("PanUKB header changed while opening: ", phenotype)
  select_idx <- match(select, input_header)
  if (anyNA(select_idx)) stop("Could not resolve selected PanUKB columns for ", phenotype)

  chunk_paths <- character()
  total_rows <- 0L
  chunk_number <- 0L
  lookup_path <- tempfile(tmpdir = out_dir, fileext = ".tsv.gz")
  on.exit(unlink(lookup_path), add = TRUE)
  repeat {
    lines <- readLines(input_con, n = panukb_batch_rows, warn = FALSE)
    if (!length(lines)) break
    chunk_number <- chunk_number + 1L
    total_rows <- total_rows + length(lines)
    chunk_path <- file.path(cache_dir, sprintf("chunk_%06d.tsv.gz", chunk_number))
    empty_path <- file.path(cache_dir, sprintf("chunk_%06d.empty", chunk_number))
    cached <- file.exists(empty_path) || (file.exists(chunk_path) && file.info(chunk_path)$size > 0 &&
      identical(as.integer(system2("gzip", c("-t", chunk_path), stdout = FALSE, stderr = FALSE)), 0L))
    if (cached) {
      message("Reusing cached PanUKB batch ", chunk_number, " for ", phenotype)
    } else {
      # fread assigns col.names after applying select, so provide the selected
      # names while keeping select_idx relative to the full source row.
      dt_chunk <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE,
                        col.names = input_header[select_idx], select = select_idx, showProgress = FALSE,
                        check.names = FALSE)
      ready <- prepare_panukb_batch(dt_chunk, n_value, phenotype, lookup_path, chain_file)
      rm(dt_chunk)
      if (nrow(ready)) {
        fwrite(ready, chunk_path, sep = "\t", na = "NA")
        if (file.exists(empty_path)) unlink(empty_path)
      } else {
        file.create(empty_path)
        if (file.exists(chunk_path)) unlink(chunk_path)
      }
      rm(ready)
    }
    if (file.exists(chunk_path)) chunk_paths <- c(chunk_paths, chunk_path)
    rm(lines)
    if (chunk_number %% 10L == 0L) {
      message("Processed/read ", chunk_number, " PanUKB batches for ", phenotype,
              " (", total_rows, " rows); completed batches are cached.")
    }
    invisible(gc(FALSE))
  }
  if (!length(chunk_paths)) stop("No usable PanUKB variants remained for ", phenotype)

  output_con <- gzfile(out_path, open = "wt")
  output_open <- TRUE
  on.exit(if (output_open) {
    try(close(output_con), silent = TRUE)
  }, add = TRUE)
  header_written <- FALSE
  variants <- 0L
  for (chunk_path in chunk_paths) {
    chunk_con <- gzfile(chunk_path, open = "rt")
    chunk_lines <- readLines(chunk_con, warn = FALSE)
    close(chunk_con)
    if (length(chunk_lines)) {
      if (header_written) chunk_lines <- chunk_lines[-1L]
      if (length(chunk_lines)) {
        writeLines(chunk_lines, output_con)
        variants <- variants + length(chunk_lines) - as.integer(!header_written)
        header_written <- TRUE
      }
    }
    rm(chunk_lines)
  }
  close(output_con)
  output_open <- FALSE
  if (!variants) stop("PanUKB batch preparation produced no usable rows for ", phenotype,
                      ". Check the selected source columns and coordinate/effect fields.")
  ids <- fread(out_path, select = "SNP", showProgress = FALSE)$SNP
  duplicate_ids <- unique(ids[duplicated(ids) | duplicated(ids, fromLast = TRUE)])
  duplicate_rows_removed <- 0L
  if (length(duplicate_ids)) {
    duplicate_rows_removed <- sum(ids %in% duplicate_ids)
    message("Excluding all ", duplicate_rows_removed, " rows for ", length(duplicate_ids),
            " rsIDs duplicated across PanUKB batches for ", phenotype,
            "; examples: ", paste(head(duplicate_ids, 5L), collapse = ", "))
    filtered_path <- paste0(out_path, ".deduplicating")
    filtered_con <- gzfile(filtered_path, open = "wt")
    header_written <- FALSE
    for (chunk_path in chunk_paths) {
      chunk_dt <- fread(chunk_path, showProgress = FALSE)
      chunk_dt <- chunk_dt[!SNP %chin% duplicate_ids]
      if (nrow(chunk_dt)) {
        write.table(chunk_dt, file = filtered_con, sep = "\t", quote = FALSE,
                    row.names = FALSE, col.names = !header_written, na = "NA")
        header_written <- TRUE
      }
      rm(chunk_dt)
    }
    close(filtered_con)
    if (!header_written) {
      unlink(filtered_path)
      stop("No unique PanUKB rsIDs remain after excluding duplicate identifiers for ", phenotype)
    }
    if (!file.rename(filtered_path, out_path)) {
      unlink(filtered_path)
      stop("Could not replace PanUKB output after removing duplicate rsIDs: ", out_path)
    }
    ids <- ids[!ids %chin% duplicate_ids]
  }
  list(path = out_path, variants = length(ids), unique_rsids = uniqueN(ids),
       duplicate_rows_removed = duplicate_rows_removed, rows_read = total_rows)
}

# Stream large FinnGen inputs too, especially lab-value files without rsIDs.
# Only one source chunk and its coordinate lookup are held in memory at a time.
prepare_finngen_in_batches <- function(path, header, select, n_value, phenotype, out_path) {
  info <- file.info(path)
  signature_fields <- list(path = normalizePath(path), size = info$size,
                           mtime = as.numeric(info$mtime), header = header,
                           select = select, n = n_value,
                           batch_rows = finngen_batch_rows,
                           mungesumstats = as.character(utils::packageVersion("MungeSumstats")))
  signature_path <- tempfile(tmpdir = out_dir)
  on.exit(unlink(signature_path), add = TRUE)
  trait_cache_dir <- file.path(out_dir, "batch_cache", safe_name(phenotype))
  dir.create(trait_cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Search prior cache-signature versions so rerunning the Slurm job can resume
  # completed chunks even if the script itself has been replaced.
  signature_versions <- 1:3
  cache_keys <- vapply(signature_versions, function(version) {
    writeBin(serialize(c(list(version = as.integer(version)), signature_fields), NULL),
             signature_path)
    unname(tools::md5sum(signature_path))
  }, character(1))
  candidate_dirs <- file.path(trait_cache_dir, cache_keys)
  completed_counts <- vapply(candidate_dirs, function(candidate) {
    if (!dir.exists(candidate)) return(0L)
    length(Sys.glob(file.path(candidate, "chunk_*.tsv.gz"))) +
      length(Sys.glob(file.path(candidate, "chunk_*.empty")))
  }, integer(1))
  reusable <- which(dir.exists(candidate_dirs) & completed_counts > 0L)
  if (length(reusable)) {
    best <- reusable[order(completed_counts[reusable],
                           signature_versions[reusable] == 1L, decreasing = TRUE)][[1L]]
    cache_key <- cache_keys[[best]]
    message("Reusing compatible FinnGen cache for ", phenotype, ": ",
            candidate_dirs[[best]], " (", completed_counts[[best]], " completed batches)")
  } else {
    cache_key <- cache_keys[[1L]]
    message("No completed FinnGen cache found for ", phenotype,
            "; starting at ", file.path(trait_cache_dir, cache_key))
  }
  cache_dir <- file.path(trait_cache_dir, cache_key)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  input_con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(try(close(input_con), silent = TRUE), add = TRUE)
  header_line <- readLines(input_con, n = 1L, warn = FALSE)
  if (!length(header_line)) stop("Could not read FinnGen header for ", phenotype)
  raw_header <- trimws(strsplit(sub("\\r$", "", header_line), "\t", fixed = TRUE)[[1L]])
  raw_header[[1L]] <- sub("^\\ufeff", "", raw_header[[1L]])
  selected_indices <- match(norm(select), norm(raw_header))
  if (anyNA(selected_indices)) {
    stop("FinnGen header changed while opening ", phenotype, ". Expected: ",
         paste(select, collapse = ", "), "; observed: ", paste(raw_header, collapse = ", "))
  }

  lookup_path <- tempfile(tmpdir = out_dir, fileext = ".tsv.gz")
  on.exit(unlink(lookup_path), add = TRUE)
  chunk_paths <- character()
  chunk_number <- 0L
  total_rows <- 0L
  repeat {
    lines <- readLines(input_con, n = finngen_batch_rows, warn = FALSE)
    if (!length(lines)) break
    chunk_number <- chunk_number + 1L
    total_rows <- total_rows + length(lines)
    chunk_path <- file.path(cache_dir, sprintf("chunk_%06d.tsv.gz", chunk_number))
    empty_path <- file.path(cache_dir, sprintf("chunk_%06d.empty", chunk_number))
    cached <- file.exists(empty_path) || (file.exists(chunk_path) &&
      file.info(chunk_path)$size > 0 &&
      identical(as.integer(system2("gzip", c("-t", chunk_path), stdout = FALSE, stderr = FALSE)), 0L))
    if (cached) {
      message("Reusing cached FinnGen batch ", chunk_number, " for ", phenotype)
    } else {
      dt <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE,
                  col.names = raw_header[selected_indices], select = selected_indices,
                  check.names = FALSE, showProgress = FALSE)
      n <- nrow(dt)
      get_batch <- function(choices) {
        column <- find_col(names(dt), choices)
        if (is.null(column)) NULL else dt[[column]]
      }
      chr <- get_batch(c("CHR", "chromosome", "chrom", "#chrom", "chr"))
      bp <- as_num(get_batch(c("BP", "POS", "base_pair_location", "pos", "position")))
      raw_id <- get_batch(c("rs_id", "rsids", "rsid", "dbsnp", "SNP", "MarkerName", "variant_id", "variant", "ID"))
      id <- if (is.null(raw_id)) rep(NA_character_, n) else normalize_rsid(raw_id)
      if (all(c("alt", "ref") %in% norm(names(dt)))) {
        a1 <- get_batch("ALT")
        a2 <- get_batch("REF")
      } else {
        a1 <- get_batch(c("A1", "effect_allele", "Allele1", "ALT", "inc_allele"))
        a2 <- get_batch(c("A2", "other_allele", "Allele2", "REF", "non_effect_allele"))
      }
      beta <- as_num(get_batch(c("BETA", "beta", "Effect", "B", "log_odds")))
      se <- as_num(get_batch(c("SE", "standard_error", "sebeta", "StdErr", "stderr")))
      p <- as_num(get_batch(c("P", "p_value", "pval", "P.value", "Pvalue")))
      if (is.null(beta) || is.null(se) || is.null(a1) || is.null(a2)) {
        stop("Missing effect, SE, or allele fields in FinnGen batch ", chunk_number,
             " for ", phenotype)
      }
      if (is.null(p)) p <- 2 * pnorm(-abs(beta / se))
      p[!is.na(p) & p <= 0] <- .Machine$double.xmin
      batch_n <- as_num(n_value)
      if (length(batch_n) != 1L || is.na(batch_n)) batch_n <- NA_real_
      out <- data.table(SNP = id,
                        CHR = if (is.null(chr)) rep(NA_character_, n) else as.character(chr),
                        BP = if (is.null(bp)) rep(NA_real_, n) else bp,
                        A1 = as.character(a1), A2 = as.character(a2),
                        BETA = beta, SE = se, P = p, N = rep(batch_n, n))
      out <- out[!is.na(A1) & nzchar(A1) & !is.na(A2) & nzchar(A2) &
                   is.finite(BETA) & is.finite(SE) & SE > 0 &
                   is.finite(P) & P > 0 & P <= 1]
      missing_id <- is.na(out$SNP) | !grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE)
      if (any(missing_id)) {
        lookup <- data.table::copy(out[missing_id])
        if (anyNA(lookup$CHR) || anyNA(lookup$BP) || any(!is.finite(lookup$BP))) {
          stop("FinnGen coordinate lookup needs valid chromosome and position in batch ",
               chunk_number, " for ", phenotype)
        }
        lookup[, SNP := NULL]
        mapped <- MungeSumstats::format_sumstats(
          path = lookup, ref_genome = "GRCh38", dbSNP = 155,
          snp_ids_are_rs_ids = FALSE, return_data = TRUE,
          return_format = "data.table", save_path = lookup_path,
          force_new = TRUE, log_folder_ind = FALSE, sort_coordinates = FALSE
        )
        mapped <- as.data.table(mapped)
        required_mapped <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
        if (!all(required_mapped %in% names(mapped))) {
          stop("MungeSumstats returned unexpected fields for FinnGen batch ",
               chunk_number, " of ", phenotype, ": ", paste(names(mapped), collapse = ", "))
        }
        if (!"N" %in% names(mapped)) mapped[, N := batch_n]
        mapped <- mapped[!is.na(SNP) & grepl("^rs[0-9]+$", SNP, ignore.case = TRUE)]
        out <- data.table::rbindlist(list(out[!missing_id], mapped), fill = TRUE, use.names = TRUE)
        rm(lookup, mapped)
      }
      out <- out[!is.na(SNP) & grepl("^rs[0-9]+$", SNP, ignore.case = TRUE)]
      out <- out[, .(SNP, CHR, BP, A1, A2, BETA, SE, P, N)]
      if (nrow(out)) {
        fwrite(out, chunk_path, sep = "\t", na = "NA")
        if (file.exists(empty_path)) unlink(empty_path)
      } else {
        file.create(empty_path)
        if (file.exists(chunk_path)) unlink(chunk_path)
      }
      rm(dt, lines, out)
      invisible(gc())
    }
    if (file.exists(chunk_path)) chunk_paths <- c(chunk_paths, chunk_path)
    rm(lines)
    if (chunk_number %% 10L == 0L) {
      message("Processed FinnGen batch ", chunk_number, " for ", phenotype,
              " (", total_rows, " source rows)")
    }
  }
  if (!length(chunk_paths)) stop("No rsID-mapped FinnGen rows were produced for ", phenotype)

  # A marker can occur more than once across FinnGen batches. Exclude every
  # row for duplicated rsIDs rather than selecting an arbitrary association.
  ids_by_chunk <- lapply(chunk_paths, function(chunk_path) {
    fread(chunk_path, select = "SNP", showProgress = FALSE)$SNP
  })
  ids <- unlist(ids_by_chunk, use.names = FALSE)
  rm(ids_by_chunk)
  duplicate_ids <- unique(ids[duplicated(ids) | duplicated(ids, fromLast = TRUE)])
  if (length(duplicate_ids)) {
    message("Excluding all rows for ", length(duplicate_ids),
            " duplicated FinnGen rsIDs in ", phenotype,
            "; examples: ", paste(head(duplicate_ids, 5L), collapse = ", "))
  }

  output_con <- gzfile(out_path, open = "wt")
  output_open <- TRUE
  on.exit(if (output_open) try(close(output_con), silent = TRUE), add = TRUE)
  header_written <- FALSE
  retained_rows <- 0L
  for (chunk_path in chunk_paths) {
    chunk_dt <- fread(chunk_path, showProgress = FALSE)
    if (length(duplicate_ids)) chunk_dt <- chunk_dt[!SNP %chin% duplicate_ids]
    if (nrow(chunk_dt)) {
      write.table(chunk_dt, file = output_con, sep = "\t", quote = FALSE,
                  row.names = FALSE, col.names = !header_written, na = "NA")
      retained_rows <- retained_rows + nrow(chunk_dt)
      header_written <- TRUE
    }
    rm(chunk_dt)
  }
  close(output_con)
  output_open <- FALSE
  if (!retained_rows) stop("No unique FinnGen rsIDs remain after duplicate removal for ", phenotype)
  list(path = out_path, variants = retained_rows,
       unique_rsids = length(ids) - sum(ids %chin% duplicate_ids), rows_read = total_rows,
       duplicate_rsids_removed = length(duplicate_ids))
}

cleaned <- list()
for (i in seq_len(nrow(ref))) {
  phenotype <- as.character(ref[[ph_col]][[i]])
  path <- trimws(as.character(ref[[path_col]][[i]]))
  if (is.na(path) || !nzchar(path)) stop("No file location for METAL phenotype: ", phenotype)
  path <- resolve_path(path, phenotype)
  if (!file.exists(path) || file.access(path, 4L) != 0L) stop("Missing/unreadable raw sumstats for ", phenotype, ": ", path)

  header <- names(fread(path, nrows = 0L, showProgress = FALSE, check.names = FALSE))
  nms <- norm(header)
  is_panukb <- grepl("pan[ _-]?ukb", paste(ref[[source_col]][[i]], ref[[cohort_col]][[i]]),
                     ignore.case = TRUE)
  is_finngen <- grepl("finngen", paste(phenotype, ref[[source_col]][[i]], ref[[cohort_col]][[i]]),
                      ignore.case = TRUE)
  has_eur_effect_stats <- is_panukb && all(norm(c("beta_EUR", "se_EUR", "pval_EUR")) %in% nms)
  beta_choices <- if (has_eur_effect_stats) c("beta_EUR", "BETA", "beta", "beta_meta", "Effect", "B", "log_odds") else c("BETA", "beta", "beta_meta", "beta_meta_hq", "Effect", "B", "log_odds")
  se_choices <- if (has_eur_effect_stats) c("se_EUR", "SE", "standard_error", "sebeta", "se_meta", "StdErr", "stderr") else c("SE", "standard_error", "sebeta", "se_meta", "se_meta_hq", "StdErr", "stderr")
  p_choices <- if (has_eur_effect_stats) c("pval_EUR", "P", "p_value", "pval", "pval_meta", "P.value", "Pvalue") else c("P", "p_value", "pval", "pval_meta", "pval_metal", "P.value", "Pvalue")
  if (is_panukb) {
    message("Using ", if (has_eur_effect_stats) "EUR-specific" else "meta-analysis",
            " effect, SE, and P columns for ", phenotype)
  }
  pick <- function(choices) {
    ix <- match(norm(choices), nms, nomatch = 0L)
    ix <- ix[ix > 0L]
    if (length(ix)) header[[ix[[1L]]]] else NULL
  }
  select <- unique(na.omit(c(
    pick(c("rs_id", "rsids", "rsid", "dbsnp")),
    pick(c("SNP", "MarkerName", "variant_id", "variant", "ID")),
    pick(c("CHR", "chromosome", "chrom", "#chrom", "chr")),
    pick(c("BP", "POS", "base_pair_location", "pos", "position")),
    pick(c("A1", "effect_allele", "Allele1", "ALT", "inc_allele")),
    pick(c("A2", "other_allele", "Allele2", "REF", "non_effect_allele")),
    pick(beta_choices),
    pick(se_choices),
    pick(p_choices)
  )))
  if (identical(phenotype, "Albumin PanUKB")) {
    message("Albumin PanUKB resolved input path: ", path)
    message("Albumin PanUKB header fields found: ", paste(header, collapse = ", "))
    message("Albumin PanUKB selected fields: ", paste(select, collapse = ", "))
  }
  if (!length(select)) stop("Could not identify summary-stat columns in: ", path)

  if (is_panukb) {
    sample_n <- as_num(ref[[n_col]][[i]])
    if (is.null(sample_n) || !length(sample_n) || is.na(sample_n[[1L]])) {
      sample_n <- NA_real_
    } else {
      sample_n <- sample_n[[1L]]
    }
    out_path <- file.path(input_dir, paste0(safe_name(phenotype), ".tsv.gz"))
    batch_result <- prepare_panukb_in_batches(
      path = path, header = header, select = select, n_value = sample_n,
      phenotype = phenotype, out_path = out_path
    )
    cleaned[[phenotype]] <- list(
      path = normalizePath(out_path), n = sample_n,
      variants = batch_result$variants, unique_rsids = batch_result$unique_rsids,
      cohort = as.character(ref[[cohort_col]][[i]]),
      source = as.character(ref[[source_col]][[i]]),
      build = "GRCh38 (batch-lifted from GRCh37)"
    )
    message("Prepared PanUKB in batches: ", phenotype, " (", batch_result$variants,
            " variants; ", batch_result$rows_read, " input rows; cache: ",
            file.path(out_dir, "batch_cache", safe_name(phenotype)), ")")
    next
  }

  if (is_finngen) {
    sample_n <- as_num(ref[[n_col]][[i]])
    if (is.null(sample_n) || !length(sample_n) || is.na(sample_n[[1L]])) {
      sample_n <- NA_real_
    } else {
      sample_n <- sample_n[[1L]]
    }
    out_path <- file.path(input_dir, paste0(safe_name(phenotype), ".tsv.gz"))
    batch_result <- prepare_finngen_in_batches(
      path = path, header = header, select = select, n_value = sample_n,
      phenotype = phenotype, out_path = out_path
    )
    cleaned[[phenotype]] <- list(
      path = normalizePath(out_path), n = sample_n,
      variants = batch_result$variants, unique_rsids = batch_result$unique_rsids,
      cohort = as.character(ref[[cohort_col]][[i]]),
      source = as.character(ref[[source_col]][[i]]), build = "GRCh38"
    )
    message("Prepared FinnGen in batches: ", phenotype, " (", batch_result$variants,
            " variants; ", batch_result$rows_read, " input rows; cache: ",
            file.path(out_dir, "batch_cache", safe_name(phenotype)),
            "; excluded duplicated rsIDs: ", batch_result$duplicate_rsids_removed, ")")
    next
  }

  dt <- fread(path, select = select, showProgress = TRUE, check.names = FALSE)
  get <- function(choices) {
    c_name <- find_col(names(dt), choices)
    if (!is.null(c_name)) dt[[c_name]] else NULL
  }
  
  # Retrieving necessary columns based on potential naming conventions (from 01_preprocess_sumstats.R checks)
  chr <- get(c("CHR", "chromosome", "chrom", "#chrom", "chr"))
  bp <- get(c("BP", "POS", "base_pair_location", "pos", "position"))
  # Prefer a dedicated rsID column. Other marker fields are accepted only if
  # their values themselves contain a valid rsID (no build-dependent fallback).
  id <- get(c("rs_id", "rsids", "rsid", "dbsnp", "SNP", "MarkerName", "variant_id", "variant", "ID"))
  if (is.null(id)) id <- rep(NA_character_, nrow(dt))
  id <- normalize_rsid(id)
  nn <- norm(names(dt))
  finn_gen_alt <- all(c("alt", "ref") %in% nn) && any(nn %in% c("rsids", "rsid"))
  a1 <- get(if (finn_gen_alt) "alt" else c("A1", "effect_allele", "Allele1", "ALT", "inc_allele"))
  a2 <- get(if (finn_gen_alt) "ref" else c("A2", "other_allele", "Allele2", "REF", "non_effect_allele"))
  beta <- as_num(get(beta_choices))
  se <- as_num(get(se_choices))
  p <- as_num(get(p_choices))
  if (is.null(beta) || is.null(se)) stop("METAL IVW requires effect and SE columns for ", phenotype)
  if (is.null(a1) || is.null(a2)) stop("METAL IVW requires effect and other alleles for ", phenotype)
  if (is.null(p)) p <- 2 * pnorm(-abs(beta / se))
  p[!is.na(p) & p <= 0] <- .Machine$double.xmin
  n_value <- as_num(ref[[n_col]][[i]])
  # A missing sheet N can become numeric(0), which makes data.table() create
  # a zero-row table. Expand the scalar (or missing value) across all variants.
  if (is.null(n_value) || !length(n_value) || is.na(n_value[[1L]])) {
    n_value <- NA_real_
  } else {
    n_value <- n_value[[1L]]
  }
  n_vector <- rep(n_value, nrow(dt))
  field_lengths <- c(
    SNP = length(id),
    CHR = if (is.null(chr)) nrow(dt) else length(chr),
    BP = if (is.null(bp)) nrow(dt) else length(bp),
    A1 = length(a1), A2 = length(a2), BETA = length(beta),
    SE = length(se), P = length(p), N = length(n_vector)
  )
  if (identical(phenotype, "Albumin PanUKB")) {
    message("Albumin PanUKB rows/field lengths: expected ", nrow(dt), "; ",
            paste(names(field_lengths), field_lengths, sep = "=", collapse = ", "))
  }
  if (any(field_lengths != nrow(dt))) {
    stop("Input-vector length mismatch for ", phenotype, ". Expected ", nrow(dt),
         " values per field; got ",
         paste(names(field_lengths), field_lengths, sep = "=", collapse = ", "))
  }
  out <- data.table(SNP = id,
                    CHR = if (is.null(chr)) rep(NA_character_, nrow(dt)) else as.character(chr),
                    BP = if (is.null(bp)) rep(NA_real_, nrow(dt)) else as_num(bp),
                    A1 = as.character(a1), A2 = as.character(a2),
                    BETA = beta, SE = se, P = p,
                    N = n_vector)
  if (nrow(out) != nrow(dt)) {
    stop("Internal row-count mismatch while preparing ", phenotype,
         ": fread read ", nrow(dt), " rows but constructed ", nrow(out), ".")
  }
  valid_alleles <- !is.na(out$A1) & nzchar(out$A1) & !is.na(out$A2) & nzchar(out$A2)
  valid_beta <- is.finite(out$BETA)
  valid_se <- is.finite(out$SE) & out$SE > 0
  valid_p <- is.finite(out$P) & out$P > 0 & out$P <= 1
  usable_stats <- valid_alleles & valid_beta & valid_se & valid_p
  if (identical(phenotype, "Albumin PanUKB")) {
    message("Albumin PanUKB first rows as read by fread:")
    print(dt[seq_len(min(nrow(dt), 5L))])
    message("Albumin PanUKB parsed first rows (ALT, REF, beta, SE, P):")
    print(data.table(
      ALT = as.character(a1)[seq_len(min(length(a1), 5L))],
      REF = as.character(a2)[seq_len(min(length(a2), 5L))],
      beta = beta[seq_len(min(length(beta), 5L))],
      SE = se[seq_len(min(length(se), 5L))],
      P = p[seq_len(min(length(p), 5L))]
    ))
    message("Albumin PanUKB rows read: ", nrow(dt),
            "; valid alleles: ", sum(valid_alleles, na.rm = TRUE),
            "; finite beta: ", sum(valid_beta, na.rm = TRUE),
            "; positive finite SE: ", sum(valid_se, na.rm = TRUE),
            "; valid P: ", sum(valid_p, na.rm = TRUE),
            "; all required fields valid: ", sum(usable_stats, na.rm = TRUE))
  }

  # PanUKB coordinates are converted from GRCh37 to GRCh38 before METAL.
  # MungeSumstats also checks allele orientation against the target reference.
  if (is_panukb) {
    if (is.null(chr) || is.null(bp)) {
      stop("PanUKB METAL input needs chromosome and position columns for liftover: ", phenotype)
    }
  }

  out <- out[usable_stats]
  if (!nrow(out)) {
    stop("No usable variants remain for ", phenotype,
         ". Inspect the preceding counts for missing/invalid alleles, beta, SE, or P values.")
  }

  # FinnGen lab-value files may have only GRCh38 coordinates. PanUKB inputs
  # may also need coordinates mapped to rsIDs before their build conversion.
  missing_rsid <- is.na(out$SNP) | !grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE)
  if (any(missing_rsid) && (is_panukb || is_finngen)) {
    lookup_input <- data.table::copy(out[missing_rsid])
    if (anyNA(lookup_input$CHR) || anyNA(lookup_input$BP) || any(!is.finite(lookup_input$BP))) {
      stop("Coordinate-based rsID lookup needs valid chromosome and position columns for ", phenotype)
    }
    lookup_build <- if (is_panukb) "GRCh37" else "GRCh38"
    message("Resolving rsIDs from ", lookup_build, " coordinates for ", phenotype)
    lookup_path <- tempfile(tmpdir = out_dir, fileext = ".tsv.gz")
    lookup_input[, SNP := NULL]
    formatted <- MungeSumstats::format_sumstats(
      path = lookup_input,
      ref_genome = lookup_build,
      dbSNP = 155,
      snp_ids_are_rs_ids = FALSE,
      return_data = TRUE,
      return_format = "data.table",
      save_path = lookup_path,
      force_new = TRUE,
      log_folder_ind = FALSE,
      sort_coordinates = FALSE
    )
    unlink(lookup_path)
    if (!data.table::is.data.table(formatted) && !is.data.frame(formatted)) {
      stop("MungeSumstats did not return a summary-statistics table for ", phenotype)
    }
    mapped <- as.data.table(formatted)
    required_lookup <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
    if (!all(required_lookup %in% names(mapped))) {
      stop("Unexpected MungeSumstats output for ", phenotype, ": ", paste(names(mapped), collapse = ", "))
    }
    valid_rsid <- !is.na(mapped$SNP) & grepl("^rs[0-9]+$", mapped$SNP, ignore.case = TRUE)
    message("Mapped ", sum(valid_rsid), " of ", nrow(lookup_input), " coordinate-only rows to rsIDs for ", phenotype)
    out <- data.table::rbindlist(list(out[!missing_rsid], mapped[valid_rsid]), fill = TRUE, use.names = TRUE)
  }
  if (!is_panukb) out <- out[!is.na(SNP) & nzchar(SNP)]
  if (!nrow(out)) stop("No variants with usable rsIDs remain for ", phenotype)

  if (is_panukb) {
    if (anyNA(out$CHR) || anyNA(out$BP)) {
      stop("PanUKB METAL input has missing chromosome or position values: ", phenotype)
    }
    required_lifted <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N")
    if (!all(required_lifted %in% names(out))) {
      stop("Unexpected MungeSumstats output for ", phenotype, ": ", paste(names(out), collapse = ", "))
    }
    if (anyNA(out$SNP) || any(!grepl("^rs[0-9]+$", out$SNP, ignore.case = TRUE))) {
      stop("GRCh37 coordinate lookup did not produce valid rsIDs for all retained rows in ", phenotype)
    }
    data.table::setcolorder(out, c("SNP", "CHR", "BP", setdiff(names(out), c("SNP", "CHR", "BP"))))
    # Lift the validated, non-empty table separately so all statistics columns
    # remain attached to their coordinates.
    out <- MungeSumstats::liftover(
      sumstats_dt = out,
      ref_genome = "GRCh37",
      convert_ref_genome = "GRCh38",
      chain_source = "ensembl",
      imputation_ind = FALSE
    )
    out <- as.data.table(out)
    if (!nrow(out) || !all(c("SNP", "CHR", "BP") %in% names(out))) {
      stop("Liftover returned no usable variants for ", phenotype)
    }
    out <- out[!is.na(SNP) & !is.na(A1) & !is.na(A2) & is.finite(BETA) &
                 is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
    message("Lifted PanUKB input to GRCh38: ", phenotype, " (", nrow(out), " variants retained)")
  }

  out <- out[, .(SNP, A1, A2, BETA, SE, P, N)]
  if (!nrow(out)) stop("No valid variant rows remain after cleaning ", phenotype,
                       ". Check that the source includes usable rsIDs, effect/other alleles, beta, and SE.")
  # Exclude every row for duplicated rsIDs. Keeping one row arbitrarily could
  # retain the wrong allele/effect when a marker is multi-allelic or repeated.
  duplicate_ids <- out[, .N, by = SNP][N > 1L, SNP]
  if (length(duplicate_ids)) {
    duplicate_rows <- out[SNP %chin% duplicate_ids, .N]
    message(
      "Removing ", duplicate_rows, " rows across ", length(duplicate_ids),
      " duplicated rsIDs for ", phenotype,
      "; examples: ", paste(head(duplicate_ids, 5L), collapse = ", ")
    )
    out <- out[!SNP %chin% duplicate_ids]
  }

  if (nrow(out) == 0L) {
    stop("No variants remain for ", phenotype,
         " after removing duplicated rsIDs and invalid summary-statistic rows.")
  }
  out_path <- file.path(input_dir, paste0(safe_name(phenotype), ".tsv.gz"))
  fwrite(out, out_path, sep = "\t", na = "NA")
  cleaned[[phenotype]] <- list(path = normalizePath(out_path), n = n_value,
                               variants = nrow(out), unique_rsids = uniqueN(out$SNP),
                               cohort = as.character(ref[[cohort_col]][[i]]), source = as.character(ref[[source_col]][[i]]),
                               build = if (is_panukb) "GRCh38 (lifted from GRCh37)" else "unchanged")
  message("Prepared ", phenotype, ": ", nrow(out), " variants with rsIDs; coordinate-only rows were excluded")
  rm(dt, out)
  invisible(gc(FALSE))
}

metal_exe <- if (file.exists(metal_bin)) normalizePath(metal_bin, mustWork = TRUE) else Sys.which(metal_bin)
if (!nzchar(metal_exe)) stop("METAL executable not found: ", metal_bin, ". Set METAL_BIN to its full path.")
combined_manifest <- list()
overlap_manifest <- list()
for (combined_name in names(pairs)) {
  members <- pairs[[combined_name]]
  prefix <- file.path(out_dir, safe_name(combined_name))
  # Calculate the exact shared rsID set before METAL and record overlap QC.
  rsids_1 <- unique(fread(cleaned[[members[[1L]]]]$path, select = "SNP", showProgress = FALSE)$SNP)
  rsids_2 <- unique(fread(cleaned[[members[[2L]]]]$path, select = "SNP", showProgress = FALSE)$SNP)
  shared_rsids <- intersect(rsids_1, rsids_2)
  overlap_manifest[[length(overlap_manifest) + 1L]] <- data.table(
    combined_phenotype = combined_name, phenotype_1 = members[[1L]], phenotype_2 = members[[2L]],
    rsids_1 = length(rsids_1), rsids_2 = length(rsids_2), shared_rsids = length(shared_rsids),
    overlap_fraction_1 = length(shared_rsids) / length(rsids_1),
    overlap_fraction_2 = length(shared_rsids) / length(rsids_2))
  message("Shared rsIDs for ", combined_name, ": ", length(shared_rsids),
          " (", length(rsids_1), " / ", length(rsids_2), " available)")
  if (!length(shared_rsids)) stop("No shared rsIDs for ", combined_name,
                                  "; check rsID availability or obtain a build-aware rsID mapping.")

  # Write both inputs restricted to the shared rsID set before running METAL.
  # Each large table is loaded and released separately to limit peak memory.
  shared_dir <- file.path(out_dir, "shared_rsid_inputs", safe_name(combined_name))
  dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
  shared_paths <- character(length(members))
  for (j in seq_along(members)) {
    dat <- fread(cleaned[[members[[j]]]]$path, showProgress = TRUE)
    dat <- dat[SNP %chin% shared_rsids]
    shared_paths[[j]] <- file.path(shared_dir, paste0(safe_name(members[[j]]), "_shared.tsv.gz"))
    fwrite(dat, shared_paths[[j]], sep = "\t", na = "NA")
    rm(dat)
    invisible(gc(FALSE))
  }

  metal_script <- file.path(out_dir, paste0(safe_name(combined_name), ".metal.txt"))
  lines <- c("SCHEME STDERR", "VERBOSE OFF", "GENOMICCONTROL OFF", "MARKER SNP",
             "ALLELE A1 A2", "EFFECT BETA", "STDERR SE")
  for (input_path in shared_paths) lines <- c(lines, paste("PROCESS", input_path))
  lines <- c(lines, paste("OUTFILE", prefix, ".tbl"), "ANALYZE", "QUIT")
  writeLines(lines, metal_script)
  log_path <- paste0(prefix, ".metal.log")
  metal_stderr <- system2(metal_exe, stdin = metal_script, stdout = log_path, stderr = TRUE)
  writeLines(metal_stderr, paste0(prefix, ".metal.stderr.log"))
  status <- attr(metal_stderr, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) stop("METAL failed for ", combined_name, "; see ", log_path, " and ", paste0(prefix, ".metal.stderr.log"))
  table_path <- paste0(prefix, ".tbl")
  if (!file.exists(table_path)) stop("METAL completed without its output table: ", table_path)
  tab <- fread(table_path, check.names = FALSE, showProgress = FALSE)
  col_snp <- find_col(names(tab), c("MarkerName", "SNP"))
  col_a1 <- find_col(names(tab), c("Allele1", "A1"))
  col_a2 <- find_col(names(tab), c("Allele2", "A2"))
  col_beta <- find_col(names(tab), c("Effect", "BETA"))
  col_se <- find_col(names(tab), c("StdErr", "SE"))
  col_p <- find_col(names(tab), c("P-value", "Pvalue", "P"))
  if (any(vapply(list(col_snp, col_a1, col_a2, col_beta, col_se, col_p), is.null, logical(1)))) {
    stop("Unexpected METAL output columns in ", table_path, ": ", paste(names(tab), collapse = ", "))
  }
  n_total <- sum(vapply(members, function(x) cleaned[[x]]$n, numeric(1)), na.rm = TRUE)
  final <- data.table(SNP = as.character(tab[[col_snp]]), A1 = as.character(tab[[col_a1]]),
                      A2 = as.character(tab[[col_a2]]), BETA = as_num(tab[[col_beta]]),
                      SE = as_num(tab[[col_se]]), P = as_num(tab[[col_p]]), N = n_total)
  final <- final[!is.na(SNP) & nzchar(SNP) & is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  if (!nrow(final)) stop("METAL returned no usable variants for ", combined_name)
  final_path <- paste0(prefix, "_combined.tsv.gz")
  fwrite(final, final_path, sep = "\t", na = "NA")
  combined_manifest[[length(combined_manifest) + 1L]] <- data.table(
    phenotype = combined_name, source = "METAL_IVW", members = paste(members, collapse = " + "),
    N = n_total, build = "GRCh38", path = normalizePath(final_path, mustWork = TRUE), variants = nrow(final))
  message("METAL IVW saved: ", combined_name, " (", nrow(final), " variants)")
}
# Writing manifest and METAL-ready cleaned data to scratch (/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw)
fwrite(rbindlist(combined_manifest), file.path(out_dir, "combined_manifest.tsv"), sep = "\t", na = "")
fwrite(rbindlist(overlap_manifest), file.path(out_dir, "rsid_overlap_qc.tsv"), sep = "\t", na = "")
fwrite(data.table(phenotype = names(cleaned), source = vapply(cleaned, `[[`, "", "source"),
                  cohort = vapply(cleaned, `[[`, "", "cohort"), variants = vapply(cleaned, `[[`, numeric(1), "variants"),
                  unique_rsids = vapply(cleaned, `[[`, numeric(1), "unique_rsids"),
                  build = vapply(cleaned, `[[`, "", "build"),
                  path = vapply(cleaned, `[[`, "", "path")),
       file.path(out_dir, "metal_inputs_manifest.tsv"), sep = "\t", na = "")
message("All METAL combinations complete: ", file.path(out_dir, "combined_manifest.tsv"))
