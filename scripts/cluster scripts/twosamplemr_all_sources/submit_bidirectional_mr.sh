#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=tsmr_all_sources
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=tsmr_all_sources_%j.out
#SBATCH --error=tsmr_all_sources_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
module load R/4.4.1-foss-2022b
export R_LIBS_USER="/vast/palmer/pi/polimanti/ellen_m/immune/r_libs/R-4.4.1"

module list 2>&1
echo "PATH=${PATH}"
command -v R || true
command -v Rscript || true
command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript unavailable after loading R/4.4.1-foss-2022b" >&2
  exit 1
}
Rscript --version

PROJECT_DIR="/vast/palmer/pi/polimanti/ellen_m/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/twosamplemr_all_sources"
OUT_DIR="${PROJECT_DIR}/results/two_sample_mr/bidirectional_all_sources"
CACHE_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/twosamplemr_all_sources_cache"

mkdir -p "${OUT_DIR}" "${CACHE_DIR}"
if [[ ! -r "${SCRIPT_DIR}/run_bidirectional_mr.R" ]]; then
  echo "Cannot read R runner: ${SCRIPT_DIR}/run_bidirectional_mr.R" >&2
  exit 1
fi

Rscript "${SCRIPT_DIR}/run_bidirectional_mr.R" "${OUT_DIR}" "${CACHE_DIR}"
