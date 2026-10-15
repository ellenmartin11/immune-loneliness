#!/bin/bash

#SBATCH --partition=day
#SBATCH --job-name=finngen_ldsc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00-08:00:00
#SBATCH --output=logs/finngen_ldsc_%j.out
#SBATCH --error=logs/finngen_ldsc_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

# Avoid inherited Python 3.8 paths
unset PYTHONPATH
unset PYTHONHOME

# Libraries required by the Python 3.10 LDSC environment
PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"

export LD_LIBRARY_PATH="${PYTHON310_LIB}:${LIBFFI_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

source /gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/activate

python --version

set -euo pipefail

PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
CSV="${PROJECT_DIR}/data/GWAS_sumstats_sources_munged.csv"

LDSC_PY="/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-python3/ldsc.py"
REF_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
W_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
OUT_DIR="${PROJECT_DIR}/results/ldsc_finngen_results"

# Create the results directory before writing sumstats_files.txt.
mkdir -p "${OUT_DIR}"


echo "Python: $(command -v python)"
python --version

# Check the input CSV and LDSC script exist.
if [[ ! -f "${CSV}" ]]; then
    echo "Cannot find CSV: ${CSV}" >&2
    exit 1
fi

if [[ ! -f "${LDSC_PY}" ]]; then
    echo "Cannot find LDSC script: ${LDSC_PY}" >&2
    exit 1
fi

# Read the CSV, select FinnGen rows, and verify each listed sumstats file.
python - "${CSV}" > "${OUT_DIR}/sumstats_files.txt" <<'PY'
import csv
import os
import re
import sys

csv_path = sys.argv[1]

def normalize(value):
    return re.sub(r"[^a-z0-9]", "", value.strip().lower())

with open(csv_path, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)

    if not reader.fieldnames:
        raise SystemExit("CSV has no header row.")

    columns = {normalize(name): name for name in reader.fieldnames}

    phenotype_col = columns.get("phenotype")
    source_col = columns.get("source")
    path_col = columns.get("path")

    if not phenotype_col or not source_col or not path_col:
        raise SystemExit(
            "Could not find Phenotype, Source, and Path columns. "
            f"CSV headers are: {reader.fieldnames}"
        )

    files = []
    seen = set()

    for row in reader:
        source = (row.get(source_col) or "").strip()
        phenotype = (row.get(phenotype_col) or "").strip()
        path = (row.get(path_col) or "").strip()

        if "finngen" not in source.lower():
            continue

        if not path:
            print(f"Skipping {phenotype}: Path is blank", file=sys.stderr)
            continue

        # Resolve relative CSV paths relative to the CSV's directory.
        if not os.path.isabs(path):
            path = os.path.join(os.path.dirname(os.path.abspath(csv_path)), path)

        path = os.path.abspath(path)

        if path not in seen:
            if not os.path.isfile(path):
                raise SystemExit(
                    f"File listed for '{phenotype}' does not exist: {path}"
                )
            files.append((phenotype, path))
            seen.add(path)

if len(files) < 2:
    raise SystemExit(
        f"Found {len(files)} existing FinnGen file(s); LDSC rg needs at least two."
    )

for phenotype, path in files:
    print(path)
    print(f"Selected: {phenotype} -> {path}", file=sys.stderr)
PY

mapfile -t SUMSTATS_FILES < "${OUT_DIR}/sumstats_files.txt"

if [[ "${#SUMSTATS_FILES[@]}" -lt 2 ]]; then
    echo "Fewer than two sumstats files were selected; LDSC rg requires at least two." >&2
    exit 1
fi

# Check that the reference LD-score and weight files exist.
if [[ ! -f "${REF_LD_CHR}1.l2.ldscore.gz" ]]; then
    echo "Cannot find reference LD scores; checked: ${REF_LD_CHR}1.l2.ldscore.gz" >&2
    exit 1
fi

if [[ ! -f "${W_LD_CHR}1.l2.ldscore.gz" ]]; then
    echo "Cannot find LD-score weights; checked: ${W_LD_CHR}1.l2.ldscore.gz" >&2
    exit 1
fi

# LDSC expects the sumstats paths as a comma-separated argument.
RG_INPUT=$(IFS=,; printf '%s' "${SUMSTATS_FILES[*]}")
echo "Running LDSC genetic correlations for ${#SUMSTATS_FILES[@]} FinnGen files."

python "${LDSC_PY}" \
    --rg "${RG_INPUT}" \
    --ref-ld-chr "${REF_LD_CHR}" \
    --w-ld-chr "${W_LD_CHR}" \
    --out "${OUT_DIR}/finngen_rg"