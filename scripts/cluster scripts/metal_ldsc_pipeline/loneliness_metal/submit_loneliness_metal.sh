#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=loneliness_metal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=loneliness_metal_%j.out
#SBATCH --error=loneliness_metal_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
module load R/4.4.1-foss-2022b
module load Metal/2020-05-05-GCC-10.2.0
# Search the project renv library first, then the shared R 4.4 user library.
export R_LIBS_USER="/vast/palmer/pi/polimanti/ellen_m/immune/renv/library/R-4.4/x86_64-pc-linux-gnu:/vast/palmer/pi/polimanti/ellen_m/immune/r_libs/R-4.4.1"

module list 2>&1
echo "PATH=${PATH}"
command -v R || true
command -v Rscript || true
command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript unavailable after loading R/4.4.1-foss-2022b" >&2
  exit 1
}
Rscript --version
Rscript -e 'cat("R library paths:\n", paste(.libPaths(), collapse = "\n"), "\n")'

PROJECT_DIR="/vast/palmer/pi/polimanti/ellen_m/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline/loneliness_metal"
OUT_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/loneliness_metal"
AOULONELY="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/UCLA_loneliness_scale_eur.tsv"
PANLONELY="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/005_Loneliness_isolation_PanUKB_clean.tsv.gz"
METAL_BIN="${METAL_BIN:-$(command -v metal || command -v METAL || true)}"

mkdir -p "${OUT_DIR}"
[[ -r "${AOULONELY}" ]] || { echo "Cannot read AoU input: ${AOULONELY}" >&2; exit 1; }
[[ -r "${PANLONELY}" ]] || { echo "Cannot read PanUKB input: ${PANLONELY}" >&2; exit 1; }
[[ -x "${METAL_BIN}" ]] || { echo "Cannot find METAL executable; set METAL_BIN." >&2; exit 1; }

export METAL_BIN
echo "AoU input: ${AOULONELY}"
echo "PanUKB input: ${PANLONELY}"
echo "METAL executable: ${METAL_BIN}"
echo "Output directory: ${OUT_DIR}"

Rscript "${SCRIPT_DIR}/prepare_loneliness_metal.R" \
  "${AOULONELY}" "${PANLONELY}" "${OUT_DIR}"
