# METAL + LDSC workflow

This folder contains the self-contained scripts for combining the specified METAL traits, preparing Primary and combined traits for LDSC, and estimating each trait’s genetic correlation with `Loneliness isolation PanUKB`.

Copy the folder to:

`/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline/`

The pipeline uses raw sumstats for METAL. PanUKB METAL inputs are lifted from GRCh37 to GRCh38 before combining with FinnGen. Primary traits follow the existing LDSC path and may reuse exact-name pre-munged files from this folder’s `already_munged.csv`.

## Run order

1. Run `submit_metal_ivw.sh` and `submit_primary_clean.sh` (either order or concurrently).
2. Run `qc_clean_sumstats.R` on both manifests as described below. Review QC reports before proceeding; this script reports but does not filter.
3. Run `submit_metal_primary_munge.sh` after both QC reviews. It prepares or reuses LDSC-munged inputs and records sign-consistency diagnostics.
4. Run `submit_metal_primary_ldsc.sh`; after the array completes, run `submit_metal_primary_ldsc_aggregate.sh` to check pair completeness and apply BH-FDR.

METAL and Primary cleaning can run in either order or concurrently. LDSC preparation requires both manifests; the correlation array requires LDSC preparation to finish; aggregation requires the full array to finish.

## Slurm job → R script map

| Slurm submission script | R script it runs | Purpose / where to edit |
|---|---|---|
| `submit_metal_ivw.sh` | `prepare_metal_ivw.R` | Selects `METAL` spreadsheet rows, prepares the two input datasets per pair, runs METAL IVW, and writes combined traits. Edit this R script for input parsing, liftover, duplicate handling, pair definitions, or METAL output handling. |
| `submit_primary_clean.sh` | `clean_sumstats.R` | Selects `Primary` rows, reuses eligible exact-name pre-munged files, and cleans raw Primary inputs. Edit this R script for source formats, column mapping, cleaning, or batch/cache behavior. |
| No dedicated Slurm wrapper currently | `qc_clean_sumstats.R` | Runs separately on each manifest to produce diagnostic QC summary and flag tables. Edit this R script to add or revise QC checks; it does not modify the sumstats. |
| `submit_metal_primary_munge.sh` | `run_metal_primary_ldsc.R` with `prepare` mode | Combines the Primary and METAL manifests, reuses or creates munged files, and records sign checks. Edit `prepare` logic in the R script for LDSC input preparation. |
| `submit_metal_primary_ldsc.sh` | `run_metal_primary_ldsc.R` with `rg` mode | Array job: runs one LDSC genetic correlation between the loneliness-isolation exposure and a comparison trait per task. Edit `rg` logic or the Slurm array/concurrency settings here. |
| `submit_metal_primary_ldsc_aggregate.sh` | `aggregate_metal_primary_ldsc.R` | Checks expected pair outputs, combines results, and applies BH-FDR. Edit this R script for aggregation or multiple-testing behavior. |

The `.sh` files define the Slurm resources, modules, paths, and which R script/mode is launched. The R files contain the analysis logic. Keep the mapping explicit if the workflow is reorganized.

## QC before LDSC

Load the same R module used by the jobs, then run `qc_clean_sumstats.R` once per manifest:

```bash
PROJECT_DIR="/gpfs/gibbs/project/polimanti/ejm96/immune"
SCRIPT_DIR="${PROJECT_DIR}/scripts/cluster scripts/metal_ldsc_pipeline"
BASE="/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline"
module load R/4.4.1-foss-2022b
Rscript "${SCRIPT_DIR}/qc_clean_sumstats.R" \
  "${BASE}/primary_clean/clean_manifest.tsv" \
  "${BASE}/primary_clean/pre_ldsc_qc" \
  "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/w_hm3.snplist"

Rscript "${SCRIPT_DIR}/qc_clean_sumstats.R" \
  "${BASE}/metal_ivw/combined_manifest.tsv" \
  "${BASE}/metal_ivw/pre_ldsc_qc" \
  "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/w_hm3.snplist"
```

Set `SCRIPT_DIR` to the copied pipeline folder and `BASE` to `/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline`. Each run writes `pre_ldsc_qc_summary.tsv` and `pre_ldsc_qc_flags.tsv`. The QC script checks IDs/duplicates, alleles, effect and association fields, sample size, coordinates, available INFO/frequency fields, sign agreement when both BETA/SE and Z are present, and HapMap3 overlap. It reports findings but does not decide universal thresholds or filter variants.

## Standardized QC: proposed staged approach

Treat QC as a documented review with source-specific thresholds, not a single pass/fail number. Apply the same framework to cleaned Primary and METAL-combined traits, while retaining distinct checks for how each was produced.

1. **Provenance and compatibility:** record source, ancestry, genome build, phenotype definition, sample/case-control counts, effect scale, and any cohort overlap. Confirm build and ancestry are appropriate for the intended LDSC reference; document any exceptions.
2. **Variant and field integrity:** check required fields, missingness, valid alleles, finite effects/SE/P values, plausible N, unique SNP IDs, and coordinate/rsID consistency. Review allele orientation and any available beta-versus-Z sign checks.
3. **Association-statistic QC:** inspect P-value distributions and precision, effect/SE consistency, allele frequencies and INFO when available, extreme effects, and variant counts before/after filtering. Compare against source documentation rather than applying one arbitrary cutoff to every study.
4. **Combination-specific QC:** for METAL traits, report input and shared-variant counts, duplicate handling, allele alignment, heterogeneity statistics, and sample-size assumptions. Check that combined effects and SE/P values are valid and that the combined file has unique IDs.
5. **LDSC readiness and diagnostics:** confirm usable HapMap3 overlap, adequate SNP counts, correct munging/sample-size settings, and inspect LDSC SNP-h2, intercept, attenuation ratio, and standard errors before interpreting genetic correlations. For later GenomicSEM, also assess the genetic covariance matrix for convergence/positive definiteness and document liability-scale inputs for binary traits.

The present QC script covers much of stages 2 and 5’s input readiness, but not every item above (notably source-specific threshold decisions, METAL heterogeneity review, or LDSC model diagnostics). The next useful addition is a single QC summary that joins the per-trait QC reports, METAL overlap/heterogeneity outputs, and LDSC h2 diagnostics, while leaving filtering decisions explicit and reviewable.

## Key outputs

- `metal_ivw/combined_manifest.tsv`, combined sumstats, and `rsid_overlap_qc.tsv`.
- `primary_clean/clean_manifest.tsv`, `primary_manifest.tsv`, and cleaned sumstats.
- `primary_clean/pre_ldsc_qc/` and `metal_ivw/pre_ldsc_qc/`: QC summaries and flags.
- `munged_ldsc/ldsc_inputs_manifest.tsv` and `ldsc_sign_consistency.tsv`.
- `results/pairwise_ldsc/metal_pipeline/`: per-pair results and the final loneliness-isolation correlation table with BH-FDR.

## Operational notes

- `prepare_metal_ivw.R` uses PanUKB’s matching EUR effect/SE/P fields when available; it does not mix EUR effects with meta-analysis P values. For `ref`/`alt` inputs, ALT is treated as the effect allele. PanUKB data are lifted to GRCh38; FinnGen coordinate-only inputs are mapped using GRCh38.
- METAL inputs are restricted to shared rsIDs. Duplicated rsIDs are excluded rather than arbitrarily selecting one row. The combined N is the sum of cohort Ns, assuming no sample overlap.
- PanUKB and coordinate-only Primary files are processed in batches with resumable scratch caches. Keep those caches between reruns. `METAL_CHUNK_ROWS`, `METAL_FINGEN_CHUNK_ROWS`, and `CLEAN_CHUNK_ROWS` adjust batch sizes.
- LDSC preparation requests 64 GB RAM; each correlation task requests 16 GB. The Slurm array concurrency is capped at four; adjust only if your allocation permits.
