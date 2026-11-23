#!/bin/bash
# One-time installation of TwoSampleMR and its required R packages.
# Submit this before submit_bidirectional_mr.sh.
# Packages are installed into the same personal R library used by the MR job.

#SBATCH --partition=day
#SBATCH --job-name=install_tsmr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=install_tsmr_%j.out
#SBATCH --error=install_tsmr_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
module load R/4.4.1-foss-2022b
R_INSTALL_LIB="/vast/palmer/pi/polimanti/ellen_m/immune/r_libs/R-4.4.1"
mkdir -p "${R_INSTALL_LIB}"
export R_LIBS_USER="${R_INSTALL_LIB}"

module list 2>&1
echo "PATH=${PATH}"
command -v Rscript
Rscript --version

Rscript -e '
  lib <- "/vast/palmer/pi/polimanti/ellen_m/immune/r_libs/R-4.4.1"
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (file.access(lib, 2L) != 0L) stop("R library is not writable: ", lib)
  .libPaths(c(lib, .libPaths()))
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages("remotes", lib = lib)
  install.packages("ieugwasr", lib = lib)
  remotes::install_github(
    "MRCIEU/TwoSampleMR",
    lib = lib,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never"
  )
  library(TwoSampleMR, lib.loc = lib)
  library(ieugwasr, lib.loc = lib)
  cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
  cat("ieugwasr:", as.character(packageVersion("ieugwasr")), "\n")
'
