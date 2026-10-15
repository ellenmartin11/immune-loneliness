#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=metal_primary_clean
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=metal_primary_clean_%j.out
#SBATCH --error=metal_primary_clean_%j.err

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
mkdir -p "${R_LIBS_USER}"
module list 2>&1
echo "PATH=${PATH}"
command -v R || true
command -v Rscript || true
command -v Rscript >/dev/null 2>&1 || { echo "Rscript unavailable after loading R/4.4.1-foss-2022b" >&2; exit 1; }
Rscript --version

PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline"
MAPPING_CSV="https://docs.google.com/spreadsheets/d/1tJhUVb1uH2Eku4P_pMAYwJCagoFiO_nYxyqmdHUgAFk/export?format=csv"
ALREADY_MUNGED_CSV="${SCRIPT_DIR}/already_munged.csv"
OUT_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean"
mkdir -p "${OUT_DIR}/tmp"
export TMPDIR="${OUT_DIR}/tmp"

if [[ ! -r "${SCRIPT_DIR}/clean_sumstats.R" ]]; then
  echo "Pipeline cleaner not found: ${SCRIPT_DIR}/clean_sumstats.R" >&2
  exit 1
fi
[[ -r "${ALREADY_MUNGED_CSV}" ]] || { echo "Cannot read ${ALREADY_MUNGED_CSV}" >&2; exit 1; }
echo "Primary cleaner script: ${SCRIPT_DIR}/clean_sumstats.R"
echo "Primary cleaner output: ${OUT_DIR}"
Rscript "${SCRIPT_DIR}/clean_sumstats.R" \
  "${MAPPING_CSV}" "${OUT_DIR}" "${ALREADY_MUNGED_CSV}" "Primary"
