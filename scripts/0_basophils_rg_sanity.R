############################################################
# Basophils PanUKB vs FinnGen: LDSC sanity check
############################################################

# Ellen Martin
# 12/10/2026
# This script is used to compute the correlation between PanUKB and FinnGen basophil count, measured differently.
# This script also uses custom functions in cleaning_functions.R to prepare PanUKB and lift build.
# Correlations are shown in the LDSC log (and sample sizes after munging).
# NOTE: for some reason, the effect was flipped for FinnGen leading to negative correlations, so I flipped it back and remunged, then ran LDSC
# High correlation indicates suitability for meta-analysis using METAL with other biomarkers. 
# In the spreadsheet, can change basophils to METAL analysis if correlation his high enough. If not, just use PanUKB.
# METAL scripts are in `cluster scripts`/metal_ldsc_pipeline

suppressPackageStartupMessages({
  library(data.table)
  library(MungeSumstats)
})

# cleaning functions
source("/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/functions/cleaning_functions.R")

panukb_raw <- "/gpfs/gibbs/pi/polimanti/Pan_UKBiobank_sumstats/continuous-30160-both_sexes.tsv.gz"
finngen_raw <- "/vast/palmer/scratch/polimanti/ejm96/raw_sumstats/lab_values_summary_stats_finngen_R12_3006315.gz"

# from spreadsheet
panukb_N <- 407277
finngen_N <- 175163

# temporary and final output locations
work_dir <- "/vast/palmer/scratch/polimanti/ejm96/basophils_rg_sanity"
results_dir <- "/gpfs/gibbs/project/polimanti/ejm96/immune/results/basophils_rg_sanity"
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # Existing LDSC software and EUR reference files on the cluster
  ldsc_python <- "/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/python"
  ldsc_dir <- "/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-CBIIT"
  munge_py <- file.path(ldsc_dir, "munge_sumstats.py")
  ldsc_py <- file.path(ldsc_dir, "ldsc.py")
  hm3_snplist <- "/gpfs/gibbs/pi/polimanti/LDSC_files/LDSCORE_w_hm3.snplist"
  ref_ld_chr <- "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
  w_ld_chr <- "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"

# Python 3.10 shared libraries required by the existing LDSC environment
python310_lib <- "/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
libffi_lib <- "/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"
Sys.unsetenv(c("PYTHONPATH", "PYTHONHOME"))
old_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
Sys.setenv(LD_LIBRARY_PATH = paste(
  c(python310_lib, libffi_lib, old_ld[nzchar(old_ld)]), collapse = ":"
))


# Detect whether input has rsIDs otherwise MungeSumstats can get from CHR and POS
has_rsid_column <- function(path) {
  hdr <- names(fread(path, nrows = 0L, showProgress = FALSE, check.names = FALSE))
  if (any(tolower(gsub("[^a-z0-9]", "", hdr)) %in% c("rsid", "rsids", "dbsnp"))) {
    return(TRUE)
  }
  sample <- fread(path, nrows = 100L, showProgress = FALSE, check.names = FALSE)
  id_col <- intersect(hdr, c("SNP", "snp", "variant_id", "MarkerName", "rs_id"))
  length(id_col) > 0L && any(grepl("^rs[0-9]+$", sample[[id_col[[1L]]]], ignore.case = TRUE))
}

# Preprocessing PanUKB
pan_clean <- file.path(
  work_dir,
  "basophils_panukb_EUR_GRCh37.clean.tsv.gz"
)

prepare_panukb(
  input_path = panukb_raw,
  output_path = pan_clean
)

# Formatting for LDSC (uses format_for_ldsc function in cleaning_functions.R)
pan_formatted <- file.path(
  work_dir,
  "basophils_panukb_EUR_GRCh38.standardized.tsv.gz"
)

finngen_formatted <- file.path(
  work_dir,
  "basophils_finngen_GRCh38.standardized.tsv.gz"
)

# PanUKB: GRCh37 -> GRCh38 DONE
format_for_ldsc(
  pan_clean,
  "PanUKB EUR basophils",
  "GRCh37",
  pan_formatted
)

# FinnGen: already GRCh38 - DONE
format_for_ldsc(
  finngen_raw,
  "FinnGen R12 basophils",
  "GRCh38",
  finngen_formatted
)

# Munging (uses run_ldsc_munge function in cleaning_functions.R)
# PanUKB
pan_munged <- run_ldsc_munge(
  pan_formatted, file.path(work_dir, "basophils_panukb_GRCh38"), panukb_N
)

# FinnGen
finngen_munged <- run_ldsc_munge(
  finngen_formatted, file.path(work_dir, "basophils_finngen_GRCh38"), finngen_N
)

# LDSC
# Estimate pairwise genetic correlation using the existing EUR LD scores.
rg_prefix <- file.path(results_dir, "basophils_finngen_vs_panukb")
rg_args <- c(
  ldsc_py,
  "--rg", paste(pan_munged, finngen_munged, sep = ","),
  "--ref-ld-chr", ref_ld_chr,
  "--w-ld-chr", w_ld_chr,
  "--out", rg_prefix
)
rg_log <- system2(ldsc_python, args = rg_args, stdout = TRUE, stderr = TRUE)
rg_status <- attr(rg_log, "status")
if (is.null(rg_status)) rg_status <- 0L
if (rg_status != 0L) {
  stop("LDSC rg failed:\n", paste(rg_log, collapse = "\n"))
}

message("Complete. LDSC output prefix: ", rg_prefix)
message("Read the .log and .results files there for rG, standard error, and P value.")


### CORRECTING SIGN REVERSAL IN FINNGEN
library(data.table)

original_path <- file.path(
  work_dir,
  "basophils_finngen_GRCh38.standardized.tsv.gz"
)

corrected_path <- file.path(
  work_dir,
  "basophils_finngen_GRCh38.corrected.tsv.gz"
)

fg_original <- fread(original_path)

# Create an independent copy
fg_corrected <- copy(fg_original)

# Reverse signed statistics
fg_corrected[, `:=`(
  BETA = -BETA,
  Z = -Z
)]

# Save corrected file
fwrite(fg_corrected, corrected_path, sep = "\t")

# check reversals
original <- fread(original_path)
corrected <- fread(corrected_path)

# Compare row-by-row (same order, same SNPs)
stopifnot(
  nrow(original) == nrow(corrected),
  identical(original$SNP, corrected$SNP)
)

# Verify exact sign reversal
stopifnot(
  isTRUE(all.equal(
    corrected$Z,
    -original$Z
  ))
)

cat("Success: all FinnGen Z-scores reversed.\n")

# re-munge FinnGen for LDSC
finngen_munged_corrected <- run_ldsc_munge(
  corrected_path,
  file.path(work_dir, "basophils_finngen_GRCh38_corrected"),
  finngen_N
)

# then re-run LDSC 
rg_prefix <- file.path(results_dir, "basophils_finngen_vs_panukb")
rg_args <- c(
  ldsc_py,
  "--rg", paste(pan_munged, finngen_munged_corrected, sep = ","), # using corrected FinnGen
  "--ref-ld-chr", ref_ld_chr,
  "--w-ld-chr", w_ld_chr,
  "--out", rg_prefix
)
rg_log <- system2(ldsc_python, args = rg_args, stdout = TRUE, stderr = TRUE)
rg_status <- attr(rg_log, "status")
if (is.null(rg_status)) rg_status <- 0L
if (rg_status != 0L) {
  stop("LDSC rg failed:\n", paste(rg_log, collapse = "\n"))
}

message("Complete. LDSC output prefix: ", rg_prefix)
message("Read the .log and .results files there for rG, standard error, and P value.")
