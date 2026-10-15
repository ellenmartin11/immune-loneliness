#!/bin/bash

#SBATCH --partition=day
#SBATCH --job-name=clean_finngen
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-04:00:00
#SBATCH --mail-type=ALL

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

#SBATCH --output=logs/clean_finngen_%j.out
#SBATCH --error=logs/clean_finngen_%j.err

echo "Starting FinnGen Cleaning"

module load R

Rscript clean_finngen_files.R