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
message("Primary cleaner version: explicit analysis-row filter")

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
  CHR = c("#chrom", "CHR", "chromosome", "chrom", "chr"),
  BP = c("POS", "BP", "base_pair_location", "position", "pos"),
  A1 = c("effect_allele", "A1", "Allele1", "ALT", "inc_allele", "EA"),
  A2 = c("other_allele", "A2", "Allele2", "REF", "non_effect_allele", "OA"),
  BETA = c("beta", "BETA", "Effect", "B", "log_odds"),
  OR = c("OR", "odds_ratio"),
  SE = c("standard_error", "SE", "sebeta", "StdErr", "stderr"),
  P = c("p_value", "P", "pval", "P.value", "Pvalue"),
  LOGP = c("neg_log_10_p_value", "mlogp", "minus_log10_p", "neg_log10_p", "LOG10P"),
  Z = c("Z", "Zscore", "Z_score", "z_score", "tstat"),
  FRQ = c("effect_allele_frequency", "EAF", "af_alt", "AF_Allele2", "FRQ", "af"),
  INFO = c("INFO", "imputation_info"),
  N = c("N", "n", "sample_size", "total_n", "nCompleteSamples", "Neff"),
  N_CASE = c("N_CASE", "N_CASES", "N_CAS", "Nca", "cases_n"),
  N_CONTROL = c("N_CONTROL", "N_CONTROLS", "N_CON", "Nco", "controls_n")
)

cleaned <- list()
headers_report <- list()
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$phenotype[[i]]
  out_path <- file.path(out_dir, sprintf("%03d_%s_clean.tsv.gz", i, safe_name(trait)))

  # Reuse a complete-looking cleaned file from an earlier run. This lets a
  # rerun resume after a failed input without reprocessing finished traits.
  if (file.exists(out_path) && !is.na(file.info(out_path)$size) && file.info(out_path)$size > 0) {
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
  if (is.null(p) && !is.null(logp)) p <- 10^(-logp)
  if (is.null(z) && !is.null(beta) && !is.null(se)) z <- beta / se
  if (is.null(z) && !is.null(beta) && !is.null(p)) z <- sign(beta) * qnorm(pmax(p / 2, .Machine$double.xmin), lower.tail = FALSE)
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

  if (all(is.na(out$SNP)) && all(!is.na(out$CHR)) && all(!is.na(out$BP))) {
    is_finngen <- grepl(
      "finngen",
      paste(trait, manifest$source[[i]], manifest$cohort[[i]]),
      ignore.case = TRUE
    )
    is_lymphocyte_gwas_catalog <- identical(trait, "Lymphocyte Count") &&
      grepl("GWAS Catalog", manifest$source[[i]], ignore.case = TRUE)
    coordinate_build <- if (is_finngen) {
      "GRCh38"
    } else if (is_lymphocyte_gwas_catalog) {
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
  if (all(is.na(out$SNP))) stop("No usable variant IDs or coordinate-based rsID mapping for ", trait)
  if (all(is.na(out$P)) || all(is.na(out$Z))) stop("Need P and signed effect information for ", trait)
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
