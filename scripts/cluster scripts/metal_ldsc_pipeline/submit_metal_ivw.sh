#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=metal_ivw
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=metal_ivw_%j.out
#SBATCH --error=metal_ivw_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
module load R/4.4.1-foss-2022b
module load Metal/2020-05-05-GCC-10.2.0
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
OUT_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw"
METAL_BIN="${METAL_BIN:-$(command -v metal || command -v METAL || true)}"
METAL_CHUNK_ROWS="${METAL_CHUNK_ROWS:-250000}"
METAL_FINGEN_CHUNK_ROWS="${METAL_FINGEN_CHUNK_ROWS:-100000}"
export METAL_CHUNK_ROWS METAL_FINGEN_CHUNK_ROWS

mkdir -p "${OUT_DIR}"
[[ -x "${METAL_BIN}" ]] || { echo "Cannot find an executable METAL binary after loading its module. Set METAL_BIN to its full path." >&2; exit 1; }
echo "METAL executable: ${METAL_BIN}"
echo "PanUKB batch size: ${METAL_CHUNK_ROWS} rows"
echo "FinnGen batch size: ${METAL_FINGEN_CHUNK_ROWS} rows"
echo "METAL preparation script: ${SCRIPT_DIR}/prepare_metal_ivw.R"
Rscript "${SCRIPT_DIR}/prepare_metal_ivw.R" "${MAPPING_CSV}" "${OUT_DIR}" "${METAL_BIN}"
