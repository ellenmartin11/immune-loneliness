#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=lcv_pairs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=lcv_pairs_%j.out
#SBATCH --error=lcv_pairs_%j.err

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
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/lcv_analysis"
INPUT_RDS="/vast/palmer/scratch/polimanti/ejm96/lcv/lcv_inputs_aligned_ordered.rds"
LCV_R_DIR="/vast/palmer/pi/polimanti/ellen_m/immune/LCV/R"
OUTPUT_DIR="${PROJECT_DIR}/results/two_sample_mr/lcv"

mkdir -p "${OUTPUT_DIR}"

for required_file in \
  "${SCRIPT_DIR}/run_lcv.R" \
  "${INPUT_RDS}" \
  "${LCV_R_DIR}/RunLCV.R" \
  "${LCV_R_DIR}/MomentFunctions.R"; do
  if [[ ! -r "${required_file}" ]]; then
    echo "Required file is missing or unreadable: ${required_file}" >&2
    exit 1
  fi
done

echo "Aligned, ordered LCV inputs: ${INPUT_RDS}"
echo "LCV scripts: ${LCV_R_DIR}"
echo "Results: ${OUTPUT_DIR}"

Rscript "${SCRIPT_DIR}/run_lcv.R" \
  "${INPUT_RDS}" \
  "${OUTPUT_DIR}" \
  "${LCV_R_DIR}"
