#!/usr/bin/env Rscript

# Re-run only the METAL stage from the existing lifted/shared input files.
# This is useful when input formatting needs correction; it avoids repeating
# PanUKB liftover and coordinate matching.

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: rerun_loneliness_metal_only.R <output directory> <METAL executable>")
}
out_dir <- args[[1L]]
metal_bin <- args[[2L]]
if (!dir.exists(out_dir)) stop("Output directory does not exist: ", out_dir)
if (!file.exists(metal_bin)) stop("METAL executable not found: ", metal_bin)

pan_file <- file.path(out_dir, "PanUKB_Loneliness_shared.tsv")
aou_file <- file.path(out_dir, "AoU_Loneliness_shared.tsv")
for (path in c(pan_file, aou_file)) {
  if (!file.exists(path)) stop("Prepared shared input not found: ", path)
}

as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
find_col <- function(nms, candidates) {
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  i <- match(norm(candidates), norm(nms), nomatch = 0L)
  i <- i[i > 0L]
  if (length(i)) nms[[i[[1L]]]] else NULL
}

# Read and rewrite each source separately with unquoted tab-delimited fields.
unquoted_paths <- character(2L)
for (i in seq_along(c(pan_file, aou_file))) {
  input <- c(pan_file, aou_file)[[i]]
  dat <- fread(input, showProgress = TRUE, check.names = FALSE)
  needed <- c("SNP", "A1", "A2", "BETA", "SE")
  if (!all(needed %in% names(dat))) {
    stop("Prepared input has unexpected columns: ", input, " (",
         paste(names(dat), collapse = ", "), ")")
  }
  unquoted_paths[[i]] <- sub("\\.tsv$", ".unquoted.tsv", input)
  fwrite(dat, unquoted_paths[[i]], sep = "\t", na = "NA", quote = FALSE)
  rm(dat)
  gc(verbose = FALSE)
}

prefix <- file.path(out_dir, "Loneliness_PanUKB_AoU_METAL_retry")
metal_script <- paste0(prefix, ".metal.txt")
writeLines(c(
  "SCHEME STDERR", "VERBOSE OFF", "GENOMICCONTROL OFF", "SEPARATOR TAB",
  "MARKER SNP", "ALLELE A1 A2", "EFFECT BETA", "STDERRLABEL SE",
  paste("PROCESS", unquoted_paths[[1L]]),
  paste("PROCESS", unquoted_paths[[2L]]),
  paste("OUTFILE", prefix, ".tbl"), "ANALYZE", "QUIT"
), metal_script)

log_file <- paste0(prefix, ".metal.log")
stderr_file <- paste0(prefix, ".metal.stderr.log")
status <- system2(metal_bin, stdin = metal_script, stdout = log_file, stderr = stderr_file)
if (!identical(as.integer(status), 0L)) {
  stop("METAL failed with exit status ", status, ". Inspect ", log_file, " and ", stderr_file)
}

tables <- list.files(out_dir,
  pattern = "^Loneliness_PanUKB_AoU_METAL_retry[0-9]*\\.tbl$", full.names = TRUE)
if (!length(tables)) stop("METAL did not create a results table. Inspect ", log_file)
table_path <- tables[[which.max(file.info(tables)$mtime)]]
metal <- fread(table_path, check.names = FALSE, showProgress = FALSE)
find_output <- function(candidates) find_col(names(metal), candidates)
marker_col <- find_output(c("MarkerName", "SNP"))
a1_col <- find_output(c("Allele1", "A1"))
a2_col <- find_output(c("Allele2", "A2"))
beta_col <- find_output(c("Effect", "BETA"))
se_col <- find_output(c("StdErr", "SE"))
p_col <- find_output(c("P-value", "Pvalue", "P"))
columns <- list(marker_col, a1_col, a2_col, beta_col, se_col, p_col)
if (any(vapply(columns, is.null, logical(1)))) {
  stop("Unexpected METAL output fields: ", paste(names(metal), collapse = ", "))
}

combined <- data.table(
  SNP = as.character(metal[[marker_col]]),
  A1 = as.character(metal[[a1_col]]), A2 = as.character(metal[[a2_col]]),
  BETA = as_num(metal[[beta_col]]), SE = as_num(metal[[se_col]]),
  P = as_num(metal[[p_col]])
)
message("METAL output rows: ", nrow(combined),
        "; finite BETA: ", sum(is.finite(combined$BETA)),
        "; positive finite SE: ", sum(is.finite(combined$SE) & combined$SE > 0),
        "; valid P: ", sum(is.finite(combined$P) & combined$P > 0 & combined$P <= 1))
combined <- combined[!is.na(SNP) & is.finite(BETA) & is.finite(SE) & SE > 0 &
                       is.finite(P) & P > 0 & P <= 1]
if (!nrow(combined)) stop("METAL returned no usable variants. Inspect its table and logs.")

combined_path <- file.path(out_dir, "Loneliness_PanUKB_AoU_METAL_combined.tsv.gz")
fwrite(combined, combined_path, sep = "\t", na = "NA")
message("Combined sumstats saved: ", combined_path)
