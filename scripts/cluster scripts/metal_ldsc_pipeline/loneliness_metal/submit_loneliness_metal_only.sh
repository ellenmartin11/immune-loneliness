#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=loneliness_metal_only
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=loneliness_metal_only_%j.out
#SBATCH --error=loneliness_metal_only_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
module load R/4.4.1-foss-2022b
module load Metal/2020-05-05-GCC-10.2.0
export R_LIBS_USER="/gpfs/gibbs/project/polimanti/ejm96/immune/r_libs/R-4.4.1"

module list 2>&1
echo "PATH=${PATH}"
command -v Rscript >/dev/null 2>&1 || {
  echo "Rscript unavailable after loading R/4.4.1-foss-2022b" >&2
  exit 1
}
Rscript --version

PROJECT_DIR="/vast/palmer/pi/polimanti/ellen_m/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline/loneliness_metal"
OUT_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/loneliness_metal"
METAL_BIN="${METAL_BIN:-$(command -v metal || command -v METAL || true)}"
[[ -x "${METAL_BIN}" ]] || { echo "Cannot find METAL executable; set METAL_BIN." >&2; exit 1; }

Rscript "${SCRIPT_DIR}/rerun_loneliness_metal_only.R" "${OUT_DIR}" "${METAL_BIN}"
