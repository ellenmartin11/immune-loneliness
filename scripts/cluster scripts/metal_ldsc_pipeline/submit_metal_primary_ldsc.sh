#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=metal_primary_ldsc
#SBATCH --array=1-100%4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=metal_primary_ldsc_%A_%a.out
#SBATCH --error=metal_primary_ldsc_%A_%a.err

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
COMBINED_MANIFEST="${BASE}/metal_ivw/combined_manifest.tsv"
CLEAN_MANIFEST="${BASE}/primary_clean/clean_manifest.tsv"
PRIMARY_MANIFEST="${BASE}/primary_clean/primary_manifest.tsv"
ALREADY_MUNGED_CSV="${SCRIPT_DIR}/already_munged.csv"
MUNGED_DIR="${BASE}/munged_ldsc"
RESULTS_DIR="${PROJECT_DIR}/results/pairwise_ldsc/metal_pipeline"

LD_ENV="/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env"
LDSC_DIR="/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-CBIIT"
PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"
PYTHON="${LD_ENV}/bin/python"
MUNGE_PY="${LDSC_DIR}/munge_sumstats.py"
LDSC_PY="${LDSC_DIR}/ldsc.py"
HM3_SNPLIST="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/w_hm3.snplist"
REF_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
W_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"

mkdir -p "${MUNGED_DIR}" "${RESULTS_DIR}"
unset PYTHONPATH PYTHONHOME
export LD_LIBRARY_PATH="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib:/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
source "${LD_ENV}/bin/activate"
for f in "${COMBINED_MANIFEST}" "${CLEAN_MANIFEST}" "${PRIMARY_MANIFEST}" "${ALREADY_MUNGED_CSV}" \
         "${PYTHON}" "${MUNGE_PY}" "${LDSC_PY}" "${HM3_SNPLIST}" "${SCRIPT_DIR}/run_metal_primary_ldsc.R" \
         "${REF_LD_CHR}1.l2.ldscore.gz" "${W_LD_CHR}1.l2.ldscore.gz"; do
  [[ -r "${f}" ]] || { echo "Required file is missing or unreadable: ${f}" >&2; exit 1; }
done
echo "Python: $("${PYTHON}" --version 2>&1)"
echo "METAL combined manifest: ${COMBINED_MANIFEST}"
echo "Primary manifest: ${PRIMARY_MANIFEST}"
echo "Results will be saved under: ${RESULTS_DIR}"
echo "This array runs one loneliness-versus-trait comparison per task (max 4 at once)."

Rscript "${SCRIPT_DIR}/run_metal_primary_ldsc.R" "rg" \
  "${COMBINED_MANIFEST}" "${CLEAN_MANIFEST}" "${PRIMARY_MANIFEST}" \
  "${ALREADY_MUNGED_CSV}" "${MUNGED_DIR}" "${RESULTS_DIR}" \
  "${PYTHON}" "${MUNGE_PY}" "${LDSC_PY}" "${HM3_SNPLIST}" \
  "${REF_LD_CHR}" "${W_LD_CHR}"
