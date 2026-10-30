#!/bin/bash
# submit_ldsr_h2.sh -- univariate LDSC h2 for the poster trait set
#
# Submit via the SLURM UI or:  sbatch submit_ldsr_h2.sh
# Smoke test from the terminal:  LDSR_FILTER="Loneliness" sbatch submit_ldsr_h2.sh
#
# Reads trait list from POSTER_MANIFEST (all phenotypes by default);
# munged sumstats from MUNGED_DIR; results in $ROOT/ldsr_h2/.

#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --job-name=ldsr_h2
#SBATCH --output=ldsr_h2_%j.out
#SBATCH --error=ldsr_h2_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################
set -euo pipefail

ROOT=/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline
SCRIPTS="/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline/munging"

export MUNGED_DIR="${MUNGED_DIR:-$ROOT/munged}"
export LDSC_DIR="${LDSC_DIR:-/gpfs/gibbs/project/polimanti/ejm96/immune/ldsc_py3}"
export POSTER_MANIFEST="${POSTER_MANIFEST:-$SCRIPTS/munge_manifest_poster.tsv}"
export OUT_DIR="${OUT_DIR:-$ROOT/ldsr_h2}"

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

Rscript "$SCRIPTS/run_ldsr_h2.R"
