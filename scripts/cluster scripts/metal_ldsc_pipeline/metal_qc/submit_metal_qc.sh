#!/bin/bash
#SBATCH --job-name=metal_qc
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH -o /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/logs/metal_qc_%j.out
#SBATCH -e /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/logs/metal_qc_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

# Deep QC of the 7 METAL IVW-combined traits (PanUKB + FinnGen).
# Outputs: <prefix>_table.tsv (one row per trait) and <prefix>_flags.tsv

module load R/4.4.1-foss-2022b  # adjust if your modules differ

PIPE=/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw
SCRIPT="/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline/metal_qc/qc_metal_deep.R"

mkdir -p "$PIPE/logs"

Rscript "$SCRIPT" \
    "$PIPE/metal_manifest.tsv" \
    "$PIPE/metal_ivw/qc_metal_deep"

echo "Done. Check metal_ivw/qc_metal_deep_table.tsv and qc_metal_deep_flags.tsv"
