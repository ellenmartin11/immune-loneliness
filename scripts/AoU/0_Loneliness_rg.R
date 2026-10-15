############################################################
# Loneliness AoU vs PanUKB: LDSC sanity check
############################################################

# Ellen Martin
# 13/10/2026
# This script is used to compute the correlation between PanUKB and AoU Loneliness
# This script also uses custom functions in cleaning_functions.R to prepare PanUKB and lift build.
# Correlations are shown in the LDSC log (and sample sizes after munging).

suppressPackageStartupMessages({
  library(data.table)
  library(MungeSumstats)
})

# cleaning functions
source("/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/functions/cleaning_functions.R")

panukb_munged <- "/gpfs/gibbs/pi/polimanti/Pan_UKBiobank_sumstats/LDSC.munged/EUR.munged/EUR_categorical-2020-both_sexes-2020.tsv.sumstats.gz"
aou_munged <- "/vast/palmer/pi/polimanti/jun_h/programs/24_sdoh/data/gwas/aou_sdoh_gwas_munged/UCLA_loneliness_scale_eur.sumstats.gz"


# temporary and final output locations
work_dir <- "/vast/palmer/scratch/polimanti/ejm96/loneliness_rg_sanity"
results_dir <- "/gpfs/gibbs/project/polimanti/ejm96/immune/results/loneliness_rg_sanity"
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



# LDSC
# Estimate pairwise genetic correlation using the existing EUR LD scores.
rg_prefix <- file.path(results_dir, "loneliness_aou_vs_panukb")
rg_args <- c(
  ldsc_py,
  "--rg", paste(aou_munged, panukb_munged, sep = ","),
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

