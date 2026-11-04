#!/bin/bash
#SBATCH --job-name=gsem_ldsc
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=gsem_ldsc_%j.out
#SBATCH --error=gsem_ldsc_%j.err

# GenomicSEM setup: multivariate LDSC over the manifest traits.
# Outputs land in $OUT_DIR: LDSCoutput.rds, S/I matrices, ldsc log.

module load StdEnv
export SLURM_EXPORT_ENV=ALL
set -euo pipefail

SCRIPTS_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/genomic_sem"
OUT_DIR="/vast/palmer/scratch/polimanti/ejm96/genomic_sem"
# ---------------------------------------------------------------

export OUT_DIR
export MANIFEST="$SCRIPTS_DIR/genomicsem_munged_manifest.tsv"

# Personal R library where GenomicSEM was installed (login-node install)
export R_LIBS_USER="$HOME/Rlibs/4.4.1"
mkdir -p "$R_LIBS_USER"

# Preflight: fail early with a clear message if GenomicSEM is missing
module load R/4.4.1-foss-2022b
Rscript -e 'if (!requireNamespace("GenomicSEM", quietly = TRUE)) stop("GenomicSEM not found in .libPaths: ", paste(.libPaths(), collapse = ": "), " - install it on a login node with: remotes::install_github(\"GenomicSEM/GenomicSEM\", lib=Sys.getenv(\"R_LIBS_USER\"))")'

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

Rscript "$SCRIPTS_DIR/genomicsem_1_ldsc_setup.R"