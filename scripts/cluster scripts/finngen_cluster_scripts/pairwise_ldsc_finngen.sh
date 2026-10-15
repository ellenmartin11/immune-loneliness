#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=finngen_ldsc_pairs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/finngen_ldsc_pairs_%j.out
#SBATCH --error=logs/finngen_ldsc_pairs_%j.err

#################################################
# please DO NOT remove the following two commands
#################################################
module load StdEnv
export SLURM_EXPORT_ENV=ALL
#################################################

set -euo pipefail

# Python libraries required by the LDSC virtual environment
unset PYTHONPATH
unset PYTHONHOME

PYTHON310_LIB="/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
LIBFFI_LIB="/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"
export LD_LIBRARY_PATH="${PYTHON310_LIB}:${LIBFFI_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

source /gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env/bin/activate

PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
CSV="${PROJECT_DIR}/data/GWAS_sumstats_sources_munged.csv"

LDSC_PY="/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-CBIIT/ldsc.py"
REF_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"
W_LD_CHR="/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/"

OUT_DIR="${PROJECT_DIR}/results/ldsc_finngen_pairwise"
PAIR_DIR="${OUT_DIR}/pair_runs"

mkdir -p "${PAIR_DIR}"

echo "Python: $(command -v python)"
python --version

# Confirm required input files exist.
if [[ ! -f "${CSV}" ]]; then
    echo "Cannot find CSV: ${CSV}" >&2
    exit 1
fi

if [[ ! -f "${LDSC_PY}" ]]; then
    echo "Cannot find LDSC script: ${LDSC_PY}" >&2
    exit 1
fi

if [[ ! -f "${REF_LD_CHR}1.l2.ldscore.gz" ]]; then
    echo "Cannot find reference LD scores: ${REF_LD_CHR}1.l2.ldscore.gz" >&2
    exit 1
fi

if [[ ! -f "${W_LD_CHR}1.l2.ldscore.gz" ]]; then
    echo "Cannot find LD-score weights: ${W_LD_CHR}1.l2.ldscore.gz" >&2
    exit 1
fi

# Check that this LDSC checkout has the genetic-correlation function.
python - "${LDSC_PY}" <<'PY'
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[1])))
import ldscore.sumstats as sumstats

print("Imported LDSC sumstats:", sumstats.__file__)
if not hasattr(sumstats, "estimate_rg"):
    raise SystemExit("This LDSC installation does not contain estimate_rg.")
PY

# Read FinnGen paths from the CSV and verify they exist.
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
        if not os.path.isabs(path):
            path = os.path.join(os.path.dirname(os.path.abspath(csv_path)), path)

        path = os.path.abspath(path)
        if path not in seen:
            if not os.path.isfile(path):
                raise SystemExit(
                    f"File listed for '{phenotype}' does not exist: {path}"
                )
            files.append(path)
            seen.add(path)

if len(files) < 2:
    raise SystemExit(
        f"Found {len(files)} FinnGen files; at least two are required."
    )

for path in files:
    print(path)

print(f"Selected {len(files)} FinnGen sumstats files.", file=sys.stderr)
PY

mapfile -t SUMSTATS_FILES < "${OUT_DIR}/sumstats_files.txt"
N="${#SUMSTATS_FILES[@]}"

if (( N < 2 )); then
    echo "Fewer than two sumstats files were selected." >&2
    exit 1
fi

# Run each unique pair once: anchor i against all files after i.
for ((i = 0; i < N - 1; i++)); do
    RG_FILES=( "${SUMSTATS_FILES[i]}" "${SUMSTATS_FILES[@]:i+1}" )
    RG_INPUT=$(IFS=,; printf '%s' "${RG_FILES[*]}")

    ANCHOR=$(basename "${SUMSTATS_FILES[i]}")
    ANCHOR="${ANCHOR%.sumstats.gz}"
    printf -v INDEX "%02d" "$i"

    echo "Running anchor ${i}/$((N - 2)): ${ANCHOR} against ${#RG_FILES[@]} traits"

    python "${LDSC_PY}" \
        --rg "${RG_INPUT}" \
        --ref-ld-chr "${REF_LD_CHR}" \
        --w-ld-chr "${W_LD_CHR}" \
        --out "${PAIR_DIR}/rg_${INDEX}_${ANCHOR}"
done

# Extract full-precision rg, SE, z, and p from the CBIIT LDSC logs.
python - "${PAIR_DIR}" "${OUT_DIR}/finngen_rg_pairwise.tsv" <<'PY'
import csv
import glob
import os
import re
import sys

pair_dir, output_path = sys.argv[1:3]
rows = []

def trait_label(path):
    name = os.path.basename(path)
    return re.sub(r"\.sumstats\.gz$", "", name)

for log_path in sorted(glob.glob(os.path.join(pair_dir, "rg_*.log"))):
    p1 = None
    current = None

    with open(log_path, encoding="utf-8") as log:
        for raw_line in log:
            line = raw_line.strip()

            if "Reading summary statistics from " in line:
                path = line.split("Reading summary statistics from ", 1)[1]
                path = path.removesuffix(" ...").strip()

                if current is None and p1 is None:
                    p1 = path
                elif current is not None and current["p2"] is None:
                    current["p2"] = path

            elif line.startswith("Computing rg for phenotype "):
                current = {
                    "p2": None,
                    "rg": None,
                    "se": None,
                    "z": None,
                    "p": None,
                }

            elif current is not None and line.startswith("Genetic Correlation:"):
                match = re.search(
                    r"Genetic Correlation:\s*(\S+)\s+\(([^)]+)\)", line
                )
                if match:
                    current["rg"] = match.group(1)
                    current["se"] = match.group(2)

            elif current is not None and line.startswith("Z-score:"):
                current["z"] = line.split(":", 1)[1].strip().split()[0]

            elif current is not None and line.startswith("P:"):
                current["p"] = line.split(":", 1)[1].strip().split()[0]

                if p1 and current["p2"]:
                    rows.append([
                        trait_label(p1),
                        trait_label(current["p2"]),
                        current["rg"],
                        current["se"],
                        current["z"],
                        current["p"],
                    ])
                current = None

if not rows:
    raise SystemExit(
        f"No genetic-correlation estimates were extracted from logs in {pair_dir}"
    )

with open(output_path, "w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["p1", "p2", "rg", "se", "z", "p"])
    writer.writerows(rows)

print(f"Wrote {len(rows)} pairwise estimates to {output_path}")
PY