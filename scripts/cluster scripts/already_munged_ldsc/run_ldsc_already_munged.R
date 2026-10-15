#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Usage: run_ldsc_already_munged.R <already_munged.csv> <results_dir> <python> <ldsc.py> <ref_ld_chr_prefix> <w_ld_chr_prefix>")
}

csv_path <- args[[1L]]
results_dir <- args[[2L]]
python <- args[[3L]]
ldsc_py <- args[[4L]]
ref_ld_chr <- args[[5L]]
w_ld_chr <- args[[6L]]
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
find_col <- function(nms, candidates) {
  idx <- match(norm(candidates), norm(nms), nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx)) nms[[idx[[1L]]]] else NULL
}
safe_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "trait"
  x
}

if (!file.exists(csv_path)) stop("already_munged.csv not found: ", csv_path)
for (p in c(python, ldsc_py, paste0(ref_ld_chr, "1.l2.ldscore.gz"), paste0(w_ld_chr, "1.l2.ldscore.gz"))) {
  if (!file.exists(p)) stop("Required file not found: ", p)
}
ref <- fread(csv_path, na.strings = c("", "NA", "N/A"), showProgress = FALSE)
phenotype_col <- find_col(names(ref), "phenotype")
source_col <- find_col(names(ref), "source")
path_col <- find_col(names(ref), "path")
if (is.null(phenotype_col) || is.null(source_col) || is.null(path_col)) {
  stop("CSV must contain phenotype, source, and path columns. Found: ", paste(names(ref), collapse = ", "))
}

traits <- data.table(
  phenotype = trimws(as.character(ref[[phenotype_col]])),
  source = trimws(as.character(ref[[source_col]])),
  path = trimws(as.character(ref[[path_col]]))
)
if (anyNA(traits$phenotype) || any(!nzchar(traits$phenotype)) || anyNA(traits$source) ||
    any(!nzchar(traits$source)) || anyNA(traits$path) || any(!nzchar(traits$path))) {
  stop("already_munged.csv contains a blank phenotype, source, or path.")
}
if (anyDuplicated(traits[, .(phenotype, source)])) {
  stop("Duplicate phenotype+source labels found; make these pairs unique before running.")
}
if (any(grepl("^https?://", traits$path, ignore.case = TRUE))) {
  stop("This sanity-check job expects local LDSC-munged files. Download remote files first; do not pass .premunged.gz files as final LDSC inputs.")
}
missing <- traits$path[!file.exists(traits$path) | file.access(traits$path, 4L) != 0L]
if (length(missing)) stop("Missing or unreadable LDSC input(s):\n", paste(unique(missing), collapse = "\n"))
if (nrow(traits) < 2L) stop("At least two LDSC-munged datasets are required for pairwise rg.")

# Check each direct input before making links: LDSC-munged files need these five fields.
for (i in seq_len(nrow(traits))) {
  hdr <- names(fread(traits$path[[i]], nrows = 0L, showProgress = FALSE))
  missing_cols <- setdiff(c("SNP", "A1", "A2", "N", "Z"), hdr)
  if (length(missing_cols)) {
    stop("Not an LDSC-ready file for ", traits$phenotype[[i]], " [", traits$source[[i]], "]: ",
         traits$path[[i]], "\nMissing columns: ", paste(missing_cols, collapse = ", "),
         "\nFinnGen .premunged.gz files are intermediate inputs and must first go through munge_sumstats.py.")
  }
}

# Give LDSC readable, collision-proof basenames while reading the original files through symlinks.
link_dir <- file.path(results_dir, "ldsc_inputs")
dir.create(link_dir, recursive = TRUE, showWarnings = FALSE)
traits[, ldsc_name := sprintf("trait_%03d_%s__%s", .I, safe_label(phenotype), safe_label(source))]
traits[, ldsc_path := NA_character_]
for (i in seq_len(nrow(traits))) {
  original <- traits$path[[i]]
  suffix <- if (grepl("\\.sumstats\\.gz$", original, ignore.case = TRUE)) ".sumstats.gz" else if (grepl("\\.sumstats$", original, ignore.case = TRUE)) ".sumstats" else ".sumstats.gz"
  link_path <- file.path(link_dir, paste0(traits$ldsc_name[[i]], suffix))
  if (file.exists(link_path) || nzchar(Sys.readlink(link_path))) unlink(link_path)
  if (!file.symlink(normalizePath(original, mustWork = TRUE), link_path)) stop("Could not create LDSC input symlink: ", link_path)
  traits$ldsc_path[[i]] <- link_path
}
fwrite(traits, file.path(results_dir, "input_manifest.tsv"), sep = "\t", na = "")

extract_rg_log_table <- function(log_path) {
  if (!file.exists(log_path)) stop("LDSC log not found: ", log_path)
  lines <- readLines(log_path, warn = FALSE)
  lines <- sub("^.*?INFO\\s+-\\s+", "", lines, perl = TRUE)
  marker <- grep("Summary of Genetic Correlation Results", lines, fixed = TRUE)
  if (!length(marker)) stop("No genetic-correlation summary table found in ", log_path)
  block <- lines[seq.int(marker[[1L]] + 1L, length(lines))]
  words <- strsplit(trimws(block), "\\s+")
  header_idx <- which(vapply(words, function(x) {
    all(c("p1", "p2", "rg", "se", "z", "p") %in% tolower(x))
  }, logical(1)))
  if (!length(header_idx)) stop("Could not identify the rg table header in ", log_path)
  header_i <- header_idx[[1L]]
  rows <- block[seq.int(header_i + 1L, length(block))]
  row_words <- strsplit(trimws(rows), "\\s+")
  keep <- vapply(row_words, function(x) {
    length(x) >= 6L && grepl("^[-+0-9.eE]+$|^NA$|^nan$", x[[3L]], ignore.case = TRUE)
  }, logical(1))
  rows <- rows[keep]
  if (!length(rows)) stop("No rg result rows could be parsed from ", log_path)
  read.table(text = paste(c(block[[header_i]], rows), collapse = "\n"),
             header = TRUE, sep = "", check.names = FALSE,
             stringsAsFactors = FALSE, na.strings = c("NA", "nan"))
}

# LDSC --rg compares the first trait with each later trait, not all possible pairs.
# Anchor each trait in turn so every unique pair is estimated exactly once.
rg_parts <- vector("list", nrow(traits) - 1L)
for (i in seq_len(nrow(traits) - 1L)) {
  other_idx <- seq.int(i + 1L, nrow(traits))
  out_prefix <- file.path(results_dir, sprintf("anchor_%03d", i))
  cmd <- c(ldsc_py, "--rg", paste(traits$ldsc_path[c(i, other_idx)], collapse = ","),
           "--ref-ld-chr", ref_ld_chr, "--w-ld-chr", w_ld_chr, "--out", out_prefix)
  console <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
  status <- attr(console, "status")
  if (is.null(status)) status <- 0L
  writeLines(console, paste0(out_prefix, ".console.log"))
  if (status != 0L) {
    stop("LDSC anchor run ", i, " failed with exit status ", status,
         ". See ", out_prefix, ".console.log")
  }
  rg_parts[[i]] <- extract_rg_log_table(paste0(out_prefix, ".log"))
  message("Completed anchor ", i, "/", nrow(traits) - 1L,
          " (", nrow(rg_parts[[i]]), " correlations).")
}
rg <- rbindlist(rg_parts, fill = TRUE)
expected_pairs <- choose(nrow(traits), 2L)
if (nrow(rg) != expected_pairs) {
  stop("Expected ", expected_pairs, " unique trait pairs but parsed ", nrow(rg),
       " from LDSC logs. Check the anchor logs before using the heatmap outputs.")
}
p1_col <- find_col(names(rg), "p1")
p2_col <- find_col(names(rg), "p2")
rg_col <- find_col(names(rg), "rg")
p_col <- find_col(names(rg), "p")
if (is.null(p1_col) || is.null(p2_col) || is.null(rg_col) || is.null(p_col)) {
  stop("Unexpected LDSC rg table columns: ", paste(names(rg), collapse = ", "))
}
lookup <- traits[, .(ldsc_name, phenotype, source)]
name_from_result <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.sumstats\\.gz$", "", x, ignore.case = TRUE)
  x <- sub("\\.sumstats$", "", x, ignore.case = TRUE)
  x
}
rg[, ldsc_p1 := name_from_result(get(p1_col))]
rg[, ldsc_p2 := name_from_result(get(p2_col))]
rg <- merge(rg, lookup, by.x = "ldsc_p1", by.y = "ldsc_name", all.x = TRUE, sort = FALSE)
setnames(rg, c("phenotype", "source"), c("phenotype_1", "source_1"))
rg <- merge(rg, lookup, by.x = "ldsc_p2", by.y = "ldsc_name", all.x = TRUE, sort = FALSE)
setnames(rg, c("phenotype", "source"), c("phenotype_2", "source_2"))
if (anyNA(rg$phenotype_1) || anyNA(rg$phenotype_2)) {
  bad <- unique(c(rg[is.na(phenotype_1), ldsc_p1], rg[is.na(phenotype_2), ldsc_p2]))
  stop("Could not map LDSC result labels back to the CSV: ", paste(bad, collapse = ", "))
}

rg[, p_fdr := NA_real_]
non_diagonal <- rg$ldsc_p1 != rg$ldsc_p2 & is.finite(as.numeric(rg[[p_col]]))
rg[non_diagonal, p_fdr := p.adjust(as.numeric(get(p_col)), method = "BH")]
fwrite(rg, file.path(results_dir, "ldsc_all_pairs_labeled_fdr.tsv"), sep = "\t", na = "")

trait_names <- traits$phenotype
names(trait_names) <- traits$ldsc_name
mat <- matrix(NA_real_, nrow = nrow(traits), ncol = nrow(traits),
              dimnames = list(traits$phenotype, traits$phenotype))
diag(mat) <- 1
for (i in seq_len(nrow(rg))) {
  r <- rg[i]
  a <- match(r$ldsc_p1, traits$ldsc_name)
  b <- match(r$ldsc_p2, traits$ldsc_name)
  if (!is.na(a) && !is.na(b)) {
    value <- suppressWarnings(as.numeric(r[[rg_col]]))
    mat[a, b] <- value
    mat[b, a] <- value
  }
}
matrix_dt <- as.data.table(mat, keep.rownames = "phenotype")
fwrite(matrix_dt, file.path(results_dir, "rg_matrix.tsv"), sep = "\t", na = "NA")

lower <- mat
lower[upper.tri(lower, diag = TRUE)] <- NA_real_
lower_dt <- as.data.table(lower, keep.rownames = "phenotype")
fwrite(lower_dt, file.path(results_dir, "rg_lower_triangle.tsv"), sep = "\t", na = "NA")

sig <- matrix(NA_real_, nrow = nrow(traits), ncol = nrow(traits),
              dimnames = list(traits$phenotype, traits$phenotype))
for (i in seq_len(nrow(rg))) {
  r <- rg[i]
  a <- match(r$ldsc_p1, traits$ldsc_name)
  b <- match(r$ldsc_p2, traits$ldsc_name)
  if (!is.na(a) && !is.na(b) && a != b && is.finite(r$p_fdr) && r$p_fdr < 0.05) {
    sig[max(a, b), min(a, b)] <- suppressWarnings(as.numeric(r[[rg_col]]))
  }
}
sig_dt <- as.data.table(sig, keep.rownames = "phenotype")
fwrite(sig_dt, file.path(results_dir, "rg_lower_triangle_fdr05.tsv"), sep = "\t", na = "NA")
message("Finished ", nrow(traits), " traits. Pair table and heatmap matrices saved under: ", results_dir)
