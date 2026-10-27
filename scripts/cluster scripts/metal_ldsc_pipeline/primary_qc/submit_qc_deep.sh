#!/bin/bash
#SBATCH --job-name=qc_deep
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=08:00:00
#SBATCH --output=qc_deep_%j.log

# Deep QC of cleaned summary statistics (parallel across traits, QC_CORES workers).
# Defaults read primary_clean/clean_manifest.tsv and write primary_clean/qc_deep_*.
# Override paths: sbatch --export=ALL,QC_MANIFEST=/path,QC_OUT=/prefix submit_qc_deep.sh

export SLURM_EXPORT_ENV=ALL
module load R/4.4.1-foss-2022b
export R_LIBS_USER=/gpfs/gibbs/project/polimanti/ejm96/immune/r_libs/R-4.4.1

PIPELINE_DIR=/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline
SCRIPT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline"

QC_MANIFEST="${QC_MANIFEST:-${PIPELINE_DIR}/primary_clean/clean_manifest.tsv}"
QC_OUT="${QC_OUT:-${PIPELINE_DIR}/primary_clean/qc_deep}"
export QC_CORES="${QC_CORES:-4}"

Rscript "${SCRIPT_DIR}/qc_sumstats_deep.R" "${QC_MANIFEST}" "${QC_OUT}"
