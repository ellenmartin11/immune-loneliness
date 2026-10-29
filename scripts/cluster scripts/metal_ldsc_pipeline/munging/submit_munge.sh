#!/bin/bash
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --job-name=munge_all
#SBATCH --output=munge_all_%j.out
#SBATCH --error=munge_all_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

# Munge all 40 traits via LDSC munge_sumstats.py, N policy from munge_manifest.tsv.
# Env-overridable; defaults match the metal_ldsc_pipeline layout.
#
# Smoke test first (recommended): run 2 traits interactively or via
#   MUNGE_FILTER="Schizophrenia PGC,Loneliness isolation PanUKB" sbatch submit_munge.sh


ROOT=/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline
SCRIPTS=/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster\ scripts/metal_ldsc_pipeline/munging

# --- Python 3.10 LDSC environment (same as previous munging/LDSC runs) ---
module load Python/3.10.8-GCCcore-12.2.0
PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"
export LD_LIBRARY_PATH="${PYTHON310_LIB}:${LIBFFI_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
# cluster modules inject a py3.8 PYTHONPATH (SciPy-bundle) that shadows the
# venv's own py3.10 packages with incompatible builds — clear it
unset PYTHONPATH PYTHONHOME
source /gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/activate

export MUNGE_MANIFEST="${MUNGE_MANIFEST:-$ROOT/munge_manifest/munge_manifest.tsv}"
export MUNGED_DIR="${MUNGED_DIR:-$ROOT/munged}"
export LDSC_DIR="${LDSC_DIR:-/gpfs/gibbs/project/polimanti/ejm96/immune/ldsc_py3}"
export HM3="${HM3:-/gpfs/gibbs/pi/polimanti/LDSC_files/LDSCORE_w_hm3.snplist}"
export MUNGE_PYTHON="${MUNGE_PYTHON:-/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/python}"
export MUNGE_CORES="${MUNGE_CORES:-6}"

# fail fast if the env is broken
"$MUNGE_PYTHON" -c "import numpy, pandas, scipy, bitarray, parse" \
  || { echo "ldsc python deps missing in $MUNGE_PYTHON"; exit 1; }

module load R/4.4.1-foss-2022b

Rscript "$SCRIPTS/run_munge.R"