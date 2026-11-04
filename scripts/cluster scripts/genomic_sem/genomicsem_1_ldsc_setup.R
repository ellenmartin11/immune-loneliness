############################################################################
### genomicsem_1_ldsc_setup.R
###
### One-time setup for the loneliness GenomicSEM adjustment models.
### Runs GenomicSEM::ldsc() (multivariate LDSC, pure R - no python) on the
### traits listed in the munged-file manifest and saves the covariance
### structure to LDSCoutput.rds.
###
### Run this ONCE via SLURM (see submit_genomicsem_ldsc.sh) or interactively:
###   Rscript genomicsem_1_ldsc_setup.R
###
### Env vars (all optional):
###   MANIFEST  - path to the munged-file manifest
###               (default: <OUT_DIR>/genomicsem_munged_manifest.tsv)
###   LD_DIR    - LD score directory (default: eur_w_ld_chr used for LDSR)
###   OUT_DIR   - output directory (default: <ROOT>/genomicsem)
###   GSEM_CHR  - number of chromosomes in the LD scores (default 22)
############################################################################

ROOT <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline"

LD_DIR  <- Sys.getenv("LD_DIR",  "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr")
OUT_DIR <- Sys.getenv("OUT_DIR", file.path(ROOT, "genomicsem"))
CHR     <- as.integer(Sys.getenv("GSEM_CHR", "22"))

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)   # ldsc() writes its log to the working directory

suppressMessages(library(GenomicSEM))

## ---- manifest: short_name, trait, path, trait_type -------------------------
## short_name  : lavaan-safe trait label (letters/digits/underscore only)
## trait       : display name (for the record)
## path        : full path to the munged .sumstats.gz
## trait_type  : case_control | quantitative
MANIFEST <- Sys.getenv("MANIFEST", file.path(OUT_DIR, "genomicsem_munged_manifest.tsv"))
stopifnot(file.exists(MANIFEST))
man <- read.delim(MANIFEST, stringsAsFactors = FALSE, colClasses = "character")

stopifnot(all(c("short_name", "trait", "path", "trait_type") %in% names(man)))
stopifnot(!anyDuplicated(man$short_name))
if (any(grepl("[^A-Za-z0-9_]", man$short_name))) {
  stop("short_name values must be lavaan-safe (letters, digits, underscore only): ",
       paste(man$short_name[grepl("[^A-Za-z0-9_]", man$short_name)], collapse = ", "))
}
if (!all(file.exists(man$path))) {
  stop("Missing munged files:\n",
       paste(man$path[!file.exists(man$path)], collapse = "\n"))
}
stopifnot(all(man$trait_type %in% c("case_control", "quantitative")))

traits <- setNames(man$path, man$short_name)

## ---- sample/population prevalence -----------------------------------------
## Case-control traits: identity-liability convention (Neff-based munging,
## Grotzinger et al. 2023) -> 0.5 / 0.5. Quantitative -> NA.
sample.prev     <- ifelse(man$trait_type == "case_control", 0.5, NA_real_)
population.prev <- sample.prev

## ---- multivariate LDSC -----------------------------------------------------
## ld/wld: same EUR LD scores used for the univariate LDSR runs.
## ldsc() reads each munged file's N column (the constant --N passed at
## munging) - EA is not in this set, so no re-munge dependency here.
LDSCoutput <- ldsc(
  traits          = unname(traits),
  sample.prev     = sample.prev,
  population.prev = population.prev,
  ld              = LD_DIR,
  wld             = LD_DIR,
  trait.names     = names(traits),
  chr             = CHR
)

saveRDS(LDSCoutput, file.path(OUT_DIR, "LDSCoutput.rds"))

## Quick-look outputs
S <- LDSCoutput$S
colnames(S) <- rownames(S) <- names(traits)
write.csv(S, file.path(OUT_DIR, "S_genetic_covariance.csv"))

I <- LDSCoutput$I
colnames(I) <- rownames(I) <- names(traits)
write.csv(I, file.path(OUT_DIR, "I_intercepts.csv"))

## Sanity check: genetic correlation between the two loneliness GWASs
if (all(c("LONELY", "LONELY_AOU") %in% colnames(S))) {
  rg_lonely <- S["LONELY", "LONELY_AOU"] /
    sqrt(S["LONELY", "LONELY"] * S["LONELY_AOU", "LONELY_AOU"])
  cat(sprintf("\nSanity check: rg(LONELY, LONELY_AOU) = %.3f\n", rg_lonely))
}

cat("\nDone. Saved:\n",
    file.path(OUT_DIR, "LDSCoutput.rds"), "\n",
    file.path(OUT_DIR, "S_genetic_covariance.csv"), "\n",
    file.path(OUT_DIR, "I_intercepts.csv"), "\n")
cat("N vector used:\n")
print(setNames(LDSCoutput$N, names(traits)))
