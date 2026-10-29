# One-off repair: reconstruct SE = BETA / Z in a cleaned sumstats file whose
# SE column is empty while BETA and Z are intact.
#
# Interactive use: set the path below, then select and run the whole script
# (or source it). Also still works headless: Rscript repair_se.R <file.tsv.gz>
#
# Writes <clean_file>.repair.tsv (plain), gzips it, then swaps it in.
# The original is kept as <clean_file>.bak.gz until you delete it manually.

path <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/022_Obesity_clean.tsv.gz"

# ---- nothing below this line needs editing ---------------------------------
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  path <- commandArgs(trailingOnly = TRUE)[[1L]]
}

suppressPackageStartupMessages(library(data.table))
if (!file.exists(path)) stop("File not found: ", path)

cat("Reading", path, "...\n")
dt <- fread(path, showProgress = FALSE, na.strings = c("NA", ""))
cat("Rows:", nrow(dt), "\n")

if (!all(c("BETA", "Z", "SE") %in% names(dt)))
  stop("File must contain BETA, Z and SE columns. Found: ", paste(names(dt), collapse = "|"))

n_z_na  <- sum(is.na(dt$Z))
n_z_zero <- sum(!is.na(dt$Z) & dt$Z == 0)
n_se_na_before <- sum(is.na(dt$SE))
cat("SE NA before repair: ", n_se_na_before, "\n",
    "Z NA: ", n_z_na, "; Z == 0: ", n_z_zero,
    " (these rows get SE = NA and are dropped at munging)\n", sep = "")

# SE = BETA / Z; 0/0 and x/0 become NaN -> set to NA explicitly
se_new <- dt$BETA / dt$Z
se_new[!is.finite(se_new)] <- NA_real_
dt[, SE := se_new]

n_se_after <- sum(is.na(dt$SE))
cat("SE NA after repair: ", n_se_after,
    " (", round(100 * n_se_after / nrow(dt), 4), "%)\n", sep = "")

# Sign audit: in a consistent file sign(Z) == sign(BETA), so BETA/Z > 0 and
# every repaired SE is positive. Negative SEs expose Z/BETA sign disagreement
# (e.g., allele flipping in the source) that the QC could not check before.
n_se_neg <- sum(!is.na(dt$SE) & dt$SE < 0)
cat("Negative SE after repair: ", n_se_neg,
    " (", round(100 * n_se_neg / nrow(dt), 4), "%)",
    if (n_se_neg / nrow(dt) > 0.001)
      " -- WARNING: Z/BETA sign disagreement, inspect before munging" else "\n",
    sep = "")

tmp <- paste0(path, ".repair.tsv")
cat("Writing", tmp, "...\n")
fwrite(dt, tmp, sep = "\t", na = "NA", quote = FALSE)
rm(dt); gc(FALSE)

cat("Gzipping...\n")
system2("gzip", c("-f", tmp))          # -> .repair.tsv.gz
gz <- paste0(tmp, ".gz")

# Verify before swapping: column count and spot-check SE == BETA/Z
chk <- fread(gz, nrows = 1000)
stopifnot(ncol(chk) == ncol(fread(path, nrows = 1)))
cat("Spot check (first 3 rows):\n")
print(chk[1:3, .(SNP, BETA, SE, Z, P)])

backup <- paste0(path, ".bak.gz")
file.rename(path, backup)
file.rename(gz, path)
cat("Done.\n",
    "  repaired file: ", path, "\n",
    "  original kept: ", backup, " (delete after verifying)\n", sep = "")
cat("Verify with the deep QC single-trait run before deleting the backup.\n")


