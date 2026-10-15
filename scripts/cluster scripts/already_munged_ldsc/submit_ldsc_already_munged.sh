#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=ldsc_munged_check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00-08:00:00
#SBATCH --output=ldsc_munged_check_%j.out
#SBATCH --error=ldsc_munged_check_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail

# This standalone job uses only the final munged paths in already_munged.csv.
# It does not clean or munge raw sumstats and does not change the main workflow.
PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/already_munged_ldsc"
CSV="${SCRIPT_DIR}/already_munged.csv"
RUNNER="${SCRIPT_DIR}/run_ldsc_already_munged.R"
LD_ENV="/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env"
LDSC_PY="/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-CBIIT/ldsc.py"
REF_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
W_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"

# Keep this sanity-check run in its own job-specific project-results folder.
RESULTS_DIR="${PROJECT_DIR}/results/pairwise_ldsc/already_munged_sanity_${SLURM_JOB_ID:-manual}"
mkdir -p "${RESULTS_DIR}"

# Use the Python 3.10 runtime libraries required by the existing LDSC environment.
unset PYTHONPATH PYTHONHOME
export LD_LIBRARY_PATH="${PYTHON310_LIB}:${LIBFFI_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
source "${LD_ENV}/bin/activate"
PYTHON="${LD_ENV}/bin/python"

for required_file in "${CSV}" "${RUNNER}" "${PYTHON}" "${LDSC_PY}" \
  "${REF_LD_CHR}1.l2.ldscore.gz" "${W_LD_CHR}1.l2.ldscore.gz"; do
  if [[ ! -r "${required_file}" ]]; then
    echo "Required file is missing or unreadable: ${required_file}" >&2
    exit 1
  fi
done

echo "Using already-munged inputs from: ${CSV}"
echo "Results will be written to: ${RESULTS_DIR}"
echo "Python: $("${PYTHON}" --version 2>&1)"

Rscript "${RUNNER}" "${CSV}" "${RESULTS_DIR}" "${PYTHON}" \
  "${LDSC_PY}" "${REF_LD_CHR}" "${W_LD_CHR}"
