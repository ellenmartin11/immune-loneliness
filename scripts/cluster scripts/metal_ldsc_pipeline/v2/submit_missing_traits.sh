#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=clean_missing_traits
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --output=clean_missing_%j.out
#SBATCH --error=clean_missing_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
module load R/4.4.1-foss-2022b
R_LIBS_USER="/gpfs/gibbs/project/polimanti/ejm96/immune/r_libs/R-4.4.1"
export R_LIBS_USER
module list 2>&1
command -v Rscript >/dev/null 2>&1 || { echo "Rscript unavailable" >&2; exit 1; }

PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline/v2"
PIPELINE_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean"
MAPPING_CSV="${PIPELINE_DIR}/missing_traits.csv"
ALREADY_MUNGED_CSV="${SCRIPT_DIR}/already_munged.csv"
OUT_DIR="${PIPELINE_DIR}/primary_clean"
QC_PASSED_TSV="${PIPELINE_DIR}/qc_deep_passed.tsv"
CLEAN_CHUNK_ROWS="${CLEAN_CHUNK_ROWS:-250000}"
export CLEAN_CHUNK_ROWS

# Route the ~1.9 GB PanUKB flat files (no rsID column, ~28M rows) through the
# batched streaming path: in-memory coordinate mapping of files this size
# needs >120 GB. Default threshold is 2 GB, just above these files.
export CLEAN_BATCH_GB="${CLEAN_BATCH_GB:-1.5}"

# Serial processing: the remaining traits are cleaned one at a time. The
# batched streaming path keeps memory low; 60 GB leaves ample headroom.
CLEAN_CORES=1
export CLEAN_CORES
echo "Parallel workers (CLEAN_CORES): ${CLEAN_CORES} (serial by design)"

mkdir -p "${OUT_DIR}/tmp"
export TMPDIR="${OUT_DIR}/tmp"

[[ -r "${SCRIPT_DIR}/clean_sumstats_v2.R" ]] || { echo "Cleaner not found: ${SCRIPT_DIR}/clean_sumstats_v2.R" >&2; exit 1; }
[[ -r "${MAPPING_CSV}" ]] || { echo "Missing-traits mapping not found: ${MAPPING_CSV}" >&2; exit 1; }
[[ -r "${ALREADY_MUNGED_CSV}" ]] || { echo "Cannot read ${ALREADY_MUNGED_CSV}" >&2; exit 1; }

# Preserve the full 26-trait manifest: the cleaner overwrites clean_manifest.tsv
# and clean_status.tsv with only what IT processed, so back both up and merge
# the subset rows back in afterwards.
cp "${OUT_DIR}/clean_manifest.tsv" "${OUT_DIR}/clean_manifest_full.tsv"
cp "${OUT_DIR}/clean_status.tsv"   "${OUT_DIR}/clean_status_full.tsv"   2>/dev/null || true

QC_PASSED_ARG=""
if [[ -r "${QC_PASSED_TSV}" ]]; then
  QC_PASSED_ARG="${QC_PASSED_TSV}"
  echo "QC-passed registry: ${QC_PASSED_TSV}"
fi

Rscript "${SCRIPT_DIR}/clean_sumstats_v2.R" \
  "${MAPPING_CSV}" "${OUT_DIR}" "${ALREADY_MUNGED_CSV}" "Primary" ${QC_PASSED_ARG}

# Merge: full manifest minus the 7 (in case any were re-cleaned) + subset results
Rscript - <<'MERGE'
suppressPackageStartupMessages(library(data.table))
out_dir <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean"
a <- fread(file.path(out_dir, "clean_manifest_full.tsv"))
b <- fread(file.path(out_dir, "clean_manifest.tsv"))
fwrite(rbindlist(list(a[!phenotype %in% b$phenotype], b), fill = TRUE),
       file.path(out_dir, "clean_manifest.tsv"))
cat("Merged manifest:", nrow(a[!phenotype %in% b$phenotype]) + nrow(b), "traits\n")
MERGE

echo "Done. Merged clean_manifest.tsv has all traits; verify the trait count above."
