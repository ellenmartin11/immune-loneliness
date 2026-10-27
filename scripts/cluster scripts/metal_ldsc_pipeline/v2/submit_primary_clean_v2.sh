#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=metal_primary_clean_v2
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
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline/v2"
MAPPING_CSV="https://docs.google.com/spreadsheets/d/1tJhUVb1uH2Eku4P_pMAYwJCagoFiO_nYxyqmdHUgAFk/export?format=csv"
ALREADY_MUNGED_CSV="${SCRIPT_DIR}/already_munged.csv"
OUT_DIR="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean"
QC_PASSED_TSV="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/qc_deep_passed.tsv"
CLEAN_CHUNK_ROWS="${CLEAN_CHUNK_ROWS:-250000}"
export CLEAN_CHUNK_ROWS

# Multicore parallelization across phenotypes. The cleaner forks one worker
# per trait (mc.preschedule=FALSE), so each worker holds one input file in
# memory: keep CLEAN_CORES <= --cpus-per-task, and raise --mem if two large
# sources (FinnGen, GWAS Catalog) land in the same wave. Set CLEAN_CORES=2
# to fall back to fully sequential processing.
CLEAN_CORES="${CLEAN_CORES:-${SLURM_CPUS_PER_TASK:-1}}"  # defaults to --cpus-per-task (2)
if [[ "${CLEAN_CORES}" -gt "${SLURM_CPUS_PER_TASK:-${CLEAN_CORES}}" ]]; then
  echo "CLEAN_CORES (${CLEAN_CORES}) exceeds --cpus-per-task; clamping." >&2
  CLEAN_CORES="${SLURM_CPUS_PER_TASK}"
fi
export CLEAN_CORES
echo "Parallel workers (CLEAN_CORES): ${CLEAN_CORES}"

mkdir -p "${OUT_DIR}/tmp"
export TMPDIR="${OUT_DIR}/tmp"

if [[ ! -r "${SCRIPT_DIR}/clean_sumstats_v2.R" ]]; then
  echo "Pipeline cleaner not found: ${SCRIPT_DIR}/clean_sumstats_v2.R" >&2
  exit 1
fi
[[ -r "${ALREADY_MUNGED_CSV}" ]] || { echo "Cannot read ${ALREADY_MUNGED_CSV}" >&2; exit 1; }
echo "Primary cleaner script: ${SCRIPT_DIR}/clean_sumstats_v2.R"
echo "Primary cleaner output: ${OUT_DIR}"
echo "Coordinate-lookup batch size: ${CLEAN_CHUNK_ROWS} rows"

# Fifth argument is optional: the QC-passed registry written by the deep QC
# script. Phenotypes listed there are skipped entirely (never reprocessed).
QC_PASSED_ARG=""
if [[ -r "${QC_PASSED_TSV}" ]]; then
  QC_PASSED_ARG="${QC_PASSED_TSV}"
  echo "QC-passed registry: ${QC_PASSED_TSV}"
else
  echo "QC-passed registry not found (run qc_sumstats_deep.R first); all traits will be evaluated."
fi

Rscript "${SCRIPT_DIR}/clean_sumstats_v2.R" \
  "${MAPPING_CSV}" "${OUT_DIR}" "${ALREADY_MUNGED_CSV}" "Primary" ${QC_PASSED_ARG}
