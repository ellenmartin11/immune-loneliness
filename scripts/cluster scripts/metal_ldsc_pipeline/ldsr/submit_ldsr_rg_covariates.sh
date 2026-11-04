#!/bin/bash
# submit_ldsr_rg.sh -- all pairwise LDSC genetic correlations (poster traits)
#
# Submit via the SLURM UI or:  sbatch submit_ldsr_rg.sh
# Smoke test from the terminal:  LDSR_FILTER="Loneliness|CRP" sbatch submit_ldsr_rg.sh
#
# Reads trait list from POSTER_MANIFEST (all phenotypes by default);
# munged sumstats from MUNGED_DIR; results (rg_long.tsv, rg_matrix.tsv)
# in $ROOT/ldsr_rg/.

#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --job-name=ldsr_rg
#SBATCH --output=ldsr_rg_%j.out
#SBATCH --error=ldsr_rg_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail

ROOT=/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline
SCRIPTS="/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline/ldsr"

export MUNGED_DIR="${MUNGED_DIR:-$ROOT/munged}"
export LDSC_DIR="${LDSC_DIR:-/gpfs/gibbs/project/polimanti/ejm96/immune/ldsc_py3}"
export POSTER_MANIFEST="${POSTER_MANIFEST:-$SCRIPTS/munge_manifest_poster.tsv}"

export OUT_DIR="$ROOT/ldsr_rg_group4"
mkdir -p "$OUT_DIR"

TRAITS_FILE="$SCRIPTS/group4_sdoh_traits.txt"
printf '%s
' \
  "Loneliness isolation PanUKB" \
  "Income" \
  "Educational Attainment" \
  "BMI GIANT" \
  "Drinks Per Week GSCAN" \
  "Neuroticism GPC" \
  "Sleep Duration PanUKB" \
  "Cigarettes per Day GSCAN" \
  > "$TRAITS_FILE"

export LDSR_TRAIT_SET=group4

# --- Python 3.10 LDSC environment (same recipe as submit_munge.sh) ---
module load Python/3.10.8-GCCcore-12.2.0
PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"
export LD_LIBRARY_PATH="${PYTHON310_LIB}:${LIBFFI_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
unset PYTHONPATH PYTHONHOME
source /gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/activate
export MUNGE_PYTHON="${MUNGE_PYTHON:-/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/python}"

"$MUNGE_PYTHON" -c "import numpy, pandas, scipy, bitarray, parse" \
  || { echo "ldsc python deps missing in $MUNGE_PYTHON"; exit 1; }

module load R/4.4.1-foss-2022b

Rscript "$SCRIPTS/run_ldsr_rg_covariates.R"