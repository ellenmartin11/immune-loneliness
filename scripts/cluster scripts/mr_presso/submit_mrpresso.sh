#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=mr_presso
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=mr_presso_%j.out
#SBATCH --error=mr_presso_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail

if [[ "$#" -ne 1 ]]; then
  echo "Usage: sbatch submit_mrpresso.sh <pair_name>" >&2
  exit 1
fi

PAIR_NAME="$1"
PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
SCRIPT="${PROJECT_DIR}/scripts/cluster scripts/mr_presso/run_mrpresso.R"

module load R/4.2.0-foss-2020b

cd "${PROJECT_DIR}"

echo "Running MR-PRESSO for pair: ${PAIR_NAME}"
echo "Rscript: $(command -v Rscript)"
Rscript --version

if [[ ! -r "${SCRIPT}" ]]; then
  echo "Cannot read R script: ${SCRIPT}" >&2
  exit 1
fi

Rscript "${SCRIPT}" "${PAIR_NAME}"