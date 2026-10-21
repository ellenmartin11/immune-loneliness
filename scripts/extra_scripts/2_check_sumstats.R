#2_check_sumstats.R

# this script is used to check how well sumstats were prepared and preprocessed for LDSC
# the cleaned data come from the run_primary_clean.sh job, which prepares all sumstats marked as 'Primary' in the sumstats sheet for LDSC
# this is a sanity check prior to running LDSC via slurm. If this all looks good, proceed to cluster scripts and run submit_ldsc_exposures.sh in generic_ldsc_cluster
# we can also separately check this worked for METAL

library(data.table)

manifest_path <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/clean_manifest.tsv"
manifest <- fread(manifest_path)

stopifnot(all(c("phenotype", "source", "input_path", "clean_path") %in% names(manifest)))
cat("Manifest columns:", paste(names(manifest), collapse = ", "), "\n")
cat("N phenotypes:", nrow(manifest), "\n\n")

for (i in seq_len(nrow(manifest))) {
  pheno <- manifest$phenotype[i]
  path  <- manifest$clean_path[i]   # switch to input_path to QC the raw files
  
  cat("==============================================================\n")
  cat("Phenotype:", pheno, "\n")
  cat("Source:   ", manifest$source[i], "\n")
  cat("Path:     ", path, "\n")
  
  if (!file.exists(path)) {
    cat("!! FILE NOT FOUND -- skipping\n\n")
    next
  }
  
  dt <- fread(path, nrows = 5)  # header + first rows only for quick QC
  cat("Columns (", ncol(dt), "):", paste(names(dt), collapse = ", "), "\n")
  cat("First 5 rows:\n")
  print(dt)
  cat("\n")
}

# NOTES/ISSUES
# BMI GIANT completely missing A1
# Crohn disease incorrect SNP formatting
# Depression PGC is missing A2
# Glycoprotein acetylation is not the right SNP format
# Hypothyroidism not the right SNP format
# Tyrosin not the right SNP format
# GSCAN shows chr before the chromosome number - needs to be stripped 
# some lowercase actg
