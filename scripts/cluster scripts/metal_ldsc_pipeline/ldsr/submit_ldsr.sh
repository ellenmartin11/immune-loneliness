#!/bin/bash
# submit_ldsr.sh -- sbatch wrapper for run_ldsr_h2.R / run_ldsr_rg.R
#
# Usage:
#   sbatch submit_ldsr.sh h2          # univariate h2 for all poster traits
#   sbatch submit_ldsr.sh rg          # all pairwise genetic correlations
#   MODE=h2 bash submit_ldsr.sh       # interactive (no sbatch)
#
# Env vars pass through to the R scripts, e.g.:
#   LDSR_FILTER="Loneliness" sbatch submit_ldsr.sh h2
#   LDSR_DRY_RUN=1 bash submit_ldsr.sh rg
#
# Layout (same as submit_munge.sh): scripts + poster_traits.txt +
# munge_manifest_poster.tsv live with this script; munged/ and ldsc_py3/
# are at the paths below; ldsr_h2/ and ldsr_rg/ are created under ROOT.

#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --job-name=ldsr
#SBATCH --output=ldsr_%j.out
#SBATCH --error=ldsr_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail
MODE="${1:-${MODE:-h2}}"

ROOT=/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline
SCRIPTS="/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline/ldsr"

export MUNGED_DIR="${MUNGED_DIR:-$ROOT/munged}"
export LDSC_DIR="${LDSC_DIR:-/gpfs/gibbs/project/polimanti/ejm96/immune/ldsc_py3}"
export LDSR_TRAITS="${LDSR_TRAITS:-$SCRIPTS/poster_traits.txt}"
export POSTER_MANIFEST="${POSTER_MANIFEST:-$SCRIPTS/munge_manifest_poster.tsv}"
export OUT_DIR="${OUT_DIR:-$ROOT/ldsr_$MODE}"

# --- Python 3.10 LDSC environment (same recipe as submit_munge.sh) ---
module load Python/3.10.8-GCCcore-12.2.0
PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"
export LD_LIBRARY_PATH="${PYTHON310_LIB}:${LIBFFI_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"


unset PYTHONPATH PYTHONHOME
source /gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/activate
export MUNGE_PYTHON="${MUNGE_PYTHON:-/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/python}"

# fail fast if the env is broken
"$MUNGE_PYTHON" -c "import numpy, pandas, scipy, bitarray, parse" \
  || { echo "ldsc python deps missing in $MUNGE_PYTHON"; exit 1; }

module load R/4.4.1-foss-2022b

case "$MODE" in
  h2) Rscript "$SCRIPTS/run_ldsr_h2.R" ;;
  rg) Rscript "$SCRIPTS/run_ldsr_rg.R" ;;
  *)  echo "Usage: submit_ldsr.sh [h2|rg]"; exit 1 ;;
esac
