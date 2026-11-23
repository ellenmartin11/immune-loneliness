#!/usr/bin/env Rscript

# Combine PanUKB and AoU loneliness summary statistics with METAL's STDERR scheme.
# PanUKB coordinates are lifted to GRCh38; both inputs are matched by GRCh38
# chromosome and position, with alleles supplied separately for METAL alignment.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: prepare_loneliness_metal.R <AoU file> <PanUKB file> <output directory>")
}

aou_path <- args[[1L]]
pan_path <- args[[2L]]
out_dir <- args[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PLINK_BUILD_FROM <- "GRCh37"
PLINK_BUILD_TO <- "GRCh38"

metal_bin <- Sys.getenv("METAL_BIN", "")
if (!nzchar(metal_bin) || !file.exists(metal_bin)) {
  metal_bin <- Sys.which("metal")
}
if (!nzchar(metal_bin) || !file.exists(metal_bin)) {
  metal_bin <- Sys.which("METAL")
}
if (!nzchar(metal_bin)) stop("Cannot find METAL executable. Load its module or set METAL_BIN.")

required_packages <- c(
  "MungeSumstats", "SNPlocs.Hsapiens.dbSNP155.GRCh37",
  "BSgenome.Hsapiens.1000genomes.hs37d5",
  "SNPlocs.Hsapiens.dbSNP155.GRCh38", "BSgenome.Hsapiens.NCBI.GRCh38"
)
# The package installation may be in the project renv library or the shared
# R 4.4 user library. Add both known locations before checking availability.
project_r_libraries <- c(
  "/vast/palmer/pi/polimanti/ellen_m/immune/renv/library/R-4.4/x86_64-pc-linux-gnu",
  "/vast/palmer/pi/polimanti/ellen_m/immune/r_libs/R-4.4.1"
)
project_r_libraries <- project_r_libraries[dir.exists(project_r_libraries)]
.libPaths(unique(c(project_r_libraries, .libPaths())))
message("R library paths: ", paste(.libPaths(), collapse = " | "))
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install these packages in the active R library first: ",
       paste(missing_packages, collapse = ", "),
       ". Package search paths were: ", paste(.libPaths(), collapse = " | "))
}
if (!file.exists(aou_path)) stop("AoU input not found: ", aou_path)
if (!file.exists(pan_path)) stop("PanUKB input not found: ", pan_path)

message("METAL executable: ", normalizePath(metal_bin))
message("AoU input: ", aou_path)
message("PanUKB input: ", pan_path)

as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
find_col <- function(nms, candidates) {
  i <- match(norm(candidates), norm(nms), nomatch = 0L)
  i <- i[i > 0L]
  if (length(i)) nms[[i[[1L]]]] else NULL
}
chr_clean <- function(x) {
  x <- sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)
  x
}
coordinate_id <- function(chr, bp) paste0("chr", chr_clean(chr), ":", as.integer(bp))

valid_alleles <- function(x) {
  x[, `:=`(A1 = toupper(as.character(A1)), A2 = toupper(as.character(A2)))]
  x[A1 %chin% c("A", "C", "G", "T") & A2 %chin% c("A", "C", "G", "T") & A1 != A2]
}

read_aou <- function(path) {
  x <- fread(path, showProgress = TRUE, check.names = FALSE)
  required <- c("CHROM", "GENPOS", "ALLELE0", "ALLELE1", "TEST", "BETA", "SE", "LOG10P")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("AoU file missing columns: ", paste(missing, collapse = ", "))

  x <- x[TEST == "ADD", .(
    CHR = chr_clean(CHROM), BP = as_num(GENPOS),
    A1 = toupper(as.character(ALLELE1)), A2 = toupper(as.character(ALLELE0)),
    BETA = as_num(BETA), SE = as_num(SE),
    P = pmax(10^(-as_num(LOG10P)), .Machine$double.xmin),
    N = if ("N" %in% names(x)) as_num(N) else NA_real_
  )]
  x <- x[!is.na(CHR) & nzchar(CHR) & is.finite(BP) & BP > 0 &
           is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  x <- valid_alleles(x)
  if (!nrow(x)) stop("No usable additive AoU variants remain.")
  x[, SNP := coordinate_id(CHR, BP)]
  x
}

read_panukb <- function(path) {
  x <- fread(path, showProgress = TRUE, check.names = FALSE)
  chr_col <- find_col(names(x), c("CHR", "chromosome", "chrom", "chr", "#chrom"))
  bp_col <- find_col(names(x), c("BP", "POS", "base_pair_location", "pos", "position", "GENPOS"))
  a1_col <- find_col(names(x), c("A1", "effect_allele", "Allele1", "ALT", "alt"))
  a2_col <- find_col(names(x), c("A2", "other_allele", "Allele2", "REF", "ref"))
  beta_col <- find_col(names(x), c("BETA", "beta", "beta_meta", "Effect"))
  se_col <- find_col(names(x), c("SE", "standard_error", "sebeta", "se_meta", "StdErr"))
  p_col <- find_col(names(x), c("P", "p_value", "pval", "pval_meta", "P.value"))
  needed <- c(chr_col, bp_col, a1_col, a2_col, beta_col, se_col, p_col)
  if (any(vapply(needed, is.null, logical(1)))) {
    stop("PanUKB file needs chromosome, position, alleles, BETA, SE, and P columns. Found: ",
         paste(names(x), collapse = ", "))
  }

  out <- x[, .(
    CHR = chr_clean(get(chr_col)), BP = as_num(get(bp_col)),
    A1 = toupper(as.character(get(a1_col))), A2 = toupper(as.character(get(a2_col))),
    BETA = as_num(get(beta_col)), SE = as_num(get(se_col)), P = as_num(get(p_col)),
    N = if ("N" %in% names(x)) as_num(N) else NA_real_
  )]
  out <- out[!is.na(CHR) & nzchar(CHR) & is.finite(BP) & BP > 0 &
               is.finite(BETA) & is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
  out <- valid_alleles(out)
  if (!nrow(out)) stop("No usable PanUKB variants remain.")
  out[, SNP := coordinate_id(CHR, BP)]
  out
}

pan <- read_panukb(pan_path)
aou <- read_aou(aou_path)
panukb_input_rows <- nrow(pan)
aou_input_rows <- nrow(aou)
message("Usable PanUKB rows before liftover: ", format(nrow(pan), big.mark = ","))
message("Usable AoU ADD rows: ", format(nrow(aou), big.mark = ","))

# Use the cached chain file so reruns do not download it again.
chain_dir <- file.path(out_dir, "reference")
dir.create(chain_dir, recursive = TRUE, showWarnings = FALSE)
chain_gz <- file.path(chain_dir, "GRCh37_to_GRCh38.chain.gz")
chain_file <- sub("\\.gz$", "", chain_gz)
if (!file.exists(chain_gz)) {
  url <- "https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz"
  status <- system2("curl", c("--location", "--fail", "--silent", "--show-error",
                               "--retry", "3", "--output", chain_gz, url))
  if (!identical(as.integer(status), 0L) || !file.exists(chain_gz)) {
    stop("Could not download the GRCh37-to-GRCh38 chain file: ", url)
  }
}
if (!file.exists(chain_file)) {
  status <- system2("gzip", c("-dc", chain_gz), stdout = chain_file)
  if (!identical(as.integer(status), 0L) || !file.exists(chain_file)) {
    stop("Could not decompress chain file: ", chain_gz)
  }
}

# MungeSumstats expects SNP, CHR, BP as the first three fields.
setcolorder(pan, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"))
message("Lifting PanUKB from ", PLINK_BUILD_FROM, " to ", PLINK_BUILD_TO, "...")
pan38 <- MungeSumstats::liftover(
  sumstats_dt = pan,
  ref_genome = PLINK_BUILD_FROM,
  convert_ref_genome = PLINK_BUILD_TO,
  chain_source = "ensembl",
  local_chain = chain_file,
  imputation_ind = FALSE
)
pan38 <- as.data.table(pan38)
if (!all(c("CHR", "BP", "A1", "A2", "BETA", "SE", "P") %in% names(pan38))) {
  stop("Unexpected PanUKB liftover output columns: ", paste(names(pan38), collapse = ", "))
}
pan38[, `:=`(CHR = chr_clean(CHR), BP = as_num(BP), A1 = toupper(A1), A2 = toupper(A2))]
pan38 <- pan38[!is.na(CHR) & is.finite(BP) & BP > 0 & A1 %chin% c("A", "C", "G", "T") &
                 A2 %chin% c("A", "C", "G", "T") & A1 != A2 & is.finite(BETA) &
                 is.finite(SE) & SE > 0 & is.finite(P) & P > 0 & P <= 1]
pan38[, SNP := coordinate_id(CHR, BP)]

# Drop all records at duplicated coordinates; choosing one could select the wrong allele.
drop_duplicate_markers <- function(x, label) {
  duplicates <- x[, .N, by = SNP][N > 1L, SNP]
  if (length(duplicates)) {
    message(label, ": excluding ", sum(x$SNP %chin% duplicates),
            " rows at ", length(duplicates), " duplicated coordinates.")
    x <- x[!SNP %chin% duplicates]
  }
  x
}
pan38 <- drop_duplicate_markers(pan38, "PanUKB")
aou <- drop_duplicate_markers(aou, "AoU")

shared <- intersect(pan38$SNP, aou$SNP)
if (!length(shared)) {
  stop("No shared GRCh38 coordinates remain. Check PanUKB build and coordinate/allele fields.")
}
message("Shared GRCh38 coordinates: ", format(length(shared), big.mark = ","))

pan38 <- pan38[SNP %chin% shared, .(SNP, A1, A2, BETA, SE, P, N)]
aou <- aou[SNP %chin% shared, .(SNP, A1, A2, BETA, SE, P, N)]
setorder(pan38, SNP)
setorder(aou, SNP)

pan_file <- file.path(out_dir, "PanUKB_Loneliness_shared.tsv")
aou_file <- file.path(out_dir, "AoU_Loneliness_shared.tsv")
fwrite(pan38, pan_file, sep = "\t", na = "NA")
fwrite(aou, aou_file, sep = "\t", na = "NA")

prefix <- file.path(out_dir, "Loneliness_PanUKB_AoU_METAL")
metal_script <- paste0(prefix, ".metal.txt")
metal_lines <- c(
  "SCHEME STDERR", "VERBOSE OFF", "GENOMICCONTROL OFF", "SEPARATOR TAB",
  "MARKER SNP", "ALLELE A1 A2", "EFFECT BETA", "STDERRLABEL SE",
  paste("PROCESS", pan_file), paste("PROCESS", aou_file),
  paste("OUTFILE", prefix, ".tbl"), "ANALYZE", "QUIT"
)
writeLines(metal_lines, metal_script)

log_file <- paste0(prefix, ".metal.log")
stderr_file <- paste0(prefix, ".metal.stderr.log")
status <- system2(metal_bin, stdin = metal_script, stdout = log_file, stderr = stderr_file)
if (!identical(as.integer(status), 0L)) {
  stop("METAL failed with exit status ", status, ". Inspect ", log_file, " and ", stderr_file)
}

possible_tables <- list.files(out_dir,
  pattern = "^Loneliness_PanUKB_AoU_METAL[0-9]*\\.tbl$", full.names = TRUE)
if (!length(possible_tables)) {
  stop("METAL did not produce a results table. See ", log_file, " and ", stderr_file)
}
table_path <- possible_tables[[which.max(file.info(possible_tables)$mtime)]]
metal <- fread(table_path, check.names = FALSE, showProgress = FALSE)
find_output <- function(candidates) find_col(names(metal), candidates)
marker_col <- find_output(c("MarkerName", "SNP"))
a1_col <- find_output(c("Allele1", "A1"))
a2_col <- find_output(c("Allele2", "A2"))
beta_col <- find_output(c("Effect", "BETA"))
se_col <- find_output(c("StdErr", "SE"))
p_col <- find_output(c("P-value", "Pvalue", "P"))
if (any(vapply(list(marker_col, a1_col, a2_col, beta_col, se_col, p_col), is.null, logical(1)))) {
  stop("Unexpected METAL output fields: ", paste(names(metal), collapse = ", "))
}

combined <- data.table(
  SNP = as.character(metal[[marker_col]]),
  A1 = as.character(metal[[a1_col]]), A2 = as.character(metal[[a2_col]]),
  BETA = as_num(metal[[beta_col]]), SE = as_num(metal[[se_col]]),
  P = as_num(metal[[p_col]])
)
combined <- combined[!is.na(SNP) & is.finite(BETA) & is.finite(SE) & SE > 0 &
                       is.finite(P) & P > 0 & P <= 1]
if (!nrow(combined)) stop("METAL completed but returned no usable combined rows.")

combined_path <- paste0(prefix, "_combined.tsv.gz")
fwrite(combined, combined_path, sep = "\t", na = "NA")
fwrite(data.table(
  trait = "Loneliness PanUKB + AoU",
  method = "METAL SCHEME STDERR",
  panukb_build_before_liftover = PLINK_BUILD_FROM,
  common_build = PLINK_BUILD_TO,
  panukb_input_rows = panukb_input_rows, aou_add_input_rows = aou_input_rows,
  shared_coordinates = length(shared), combined_variants = nrow(combined),
  output = normalizePath(combined_path)
), file.path(out_dir, "Loneliness_PanUKB_AoU_METAL_summary.tsv"), sep = "\t")

message("Combined loneliness sumstats saved: ", combined_path)
message("METAL output table: ", table_path)
