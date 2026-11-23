#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=cause_metal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=cause_metal_%j.out
#SBATCH --error=cause_metal_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail

# Use the R 4.2 environment that contains CAUSE and loo 2.9.0.
module load R/4.2.0-foss-2020b
export R_LIBS_USER="/vast/palmer/pi/polimanti/ellen_m/immune/renv/library/R-4.2/x86_64-pc-linux-gnu"

module list 2>&1
echo "PATH=${PATH}"
command -v R || true
command -v Rscript || true
command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript unavailable after loading R/4.2.0-foss-2020b" >&2
  exit 1
}
Rscript --version

PROJECT_DIR="/vast/palmer/pi/polimanti/ellen_m/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/cause_analysis"
OUT_DIR="${PROJECT_DIR}/results/two_sample_mr/cause/metal_sensitivity"

mkdir -p "${OUT_DIR}"

if [[ ! -r "${SCRIPT_DIR}/cause_loneliness_metal.R" ]]; then
  echo "Cannot read R runner: ${SCRIPT_DIR}/cause_loneliness_metal.R" >&2
  exit 1
fi

Rscript "${SCRIPT_DIR}/cause_loneliness_metal.R" "${OUT_DIR}"
