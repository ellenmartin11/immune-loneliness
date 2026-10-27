#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=ldsc_rg_aggregate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=ldsc_rg_aggregate_%j.out
#SBATCH --error=ldsc_rg_aggregate_%j.err

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
echo "PATH=${PATH}"
command -v R || true
command -v Rscript || true
command -v Rscript >/dev/null 2>&1 || { echo "Rscript unavailable after loading R/4.4.1-foss-2022b" >&2; exit 1; }
Rscript --version

PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline"
BASE="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline"
MUNGED_DIR="${BASE}/munged_ldsc"
RESULTS_DIR="${PROJECT_DIR}/results/pairwise_ldsc/metal_pipeline"

Rscript "${SCRIPT_DIR}/aggregate_metal_primary_ldsc.R" "${MUNGED_DIR}" "${RESULTS_DIR}"
