# METAL + LDSC workflow

This folder contains the self-contained scripts for combining the specified METAL traits, preparing Primary and combined traits for LDSC, and estimating each trait’s genetic correlation with `Loneliness isolation PanUKB`.

The pipeline uses raw sumstats for METAL. PanUKB METAL inputs are lifted from GRCh37 to GRCh38 before combining with FinnGen. Primary traits follow the existing LDSC path and may reuse exact-name pre-munged files from this folder’s `already_munged.csv`.

## Run order

1. Run `submit_metal_ivw.sh` and `submit_primary_clean.sh` (either order or concurrently).
2. In a cluster RStudio session, run `2_check_sumstats_deep.R` on the Primary manifest. Review the table and flags; this script reports QC but does not filter variants.
3. Rerun `submit_primary_clean_v2.sh`. It reads `primary_clean/qc_deep_passed.tsv`, reuses completed outputs for passed traits, and retries flagged or otherwise unfinished traits. Run deep QC again on the updated manifest to review the results of those retries.
4. Run deep QC on the METAL combined manifest in the same RStudio session and review those reports. The Primary passed-trait registry is only used by the Primary cleaner.
5. Once the QC reports are satisfactory, run `submit_metal_primary_munge.sh`. It prepares or reuses LDSC-munged inputs and records sign-consistency diagnostics.
6. Run `submit_metal_primary_ldsc.sh`; after the array completes, run `submit_metal_primary_ldsc_aggregate.sh` to check pair completeness and apply BH-FDR.

METAL and Primary cleaning can run in either order or concurrently. LDSC preparation requires both manifests; the correlation array requires LDSC preparation to finish; aggregation requires the full array to finish.

## Slurm job → R script map

| Slurm submission script | R script it runs | Purpose / where to edit |
|---|---|---|
| `submit_metal_ivw.sh` | `prepare_metal_ivw.R` | Selects `METAL` spreadsheet rows, prepares the two input datasets per pair, runs METAL IVW, and writes combined traits. |
| `submit_primary_clean_v2.sh` | `clean_sumstats_v2.R` | Selects `Primary` rows, reuses eligible exact-name pre-munged files, and cleans raw Primary inputs. |
| Run interactively in cluster RStudio | `qc_clean_sumstats.R` | Deep-checks a manifest, writing a QC table, flags, and a passed-phenotype registry. Run it after cleaning and before LDSC preparation. It does not modify sumstats. |
| `submit_metal_primary_munge.sh` | `run_metal_primary_ldsc.R` with `prepare` mode | Combines the Primary and METAL manifests, reuses or creates munged files, and records sign checks. Edit `prepare` logic in the R script for LDSC input preparation. |
| `submit_metal_primary_ldsc.sh` | `run_metal_primary_ldsc.R` with `rg` mode | Array job: runs one LDSC genetic correlation between the loneliness-isolation exposure and a comparison trait per task. Edit `rg` logic or the Slurm array/concurrency settings here. |
| `submit_metal_primary_ldsc_aggregate.sh` | `aggregate_metal_primary_ldsc.R` | Checks expected pair outputs, combines results, and applies BH-FDR. Edit this R script for aggregation or multiple-testing behavior. |

The `.sh` files define the Slurm resources, modules, paths, and which R script/mode is launched. The R files contain the analysis logic. Keep the mapping explicit if the workflow is reorganized.


## Key outputs

- `metal_ivw/combined_manifest.tsv`, combined sumstats, and `rsid_overlap_qc.tsv`.
- `primary_clean/clean_manifest.tsv`, `primary_manifest.tsv`, and cleaned sumstats.
- `primary_clean/qc_deep_table.tsv`, `qc_deep_flags.tsv`, and `qc_deep_passed.tsv`; `metal_ivw/qc_deep_table.tsv`, `qc_deep_flags.tsv`, and `qc_deep_passed.tsv`.
- `munged_ldsc/ldsc_inputs_manifest.tsv` and `ldsc_sign_consistency.tsv`.
- `results/pairwise_ldsc/metal_pipeline/`: per-pair results and the final loneliness-isolation correlation table with BH-FDR.

## Operational notes

- `prepare_metal_ivw.R` uses PanUKB’s matching EUR effect/SE/P fields when available; it does not mix EUR effects with meta-analysis P values. For `ref`/`alt` inputs, ALT is treated as the effect allele. PanUKB data are lifted to GRCh38; FinnGen coordinate-only inputs are mapped using GRCh38.
- METAL inputs are restricted to shared rsIDs. Duplicated rsIDs are excluded rather than arbitrarily selecting one row. The combined N is the sum of cohort Ns, assuming no sample overlap.
- PanUKB and coordinate-only Primary files are processed in batches with resumable scratch caches. Keep those caches between reruns. `METAL_CHUNK_ROWS`, `METAL_FINGEN_CHUNK_ROWS`, and `CLEAN_CHUNK_ROWS` adjust batch sizes.
- LDSC preparation requests 64 GB RAM; each correlation task requests 16 GB. The Slurm array concurrency is capped at four.
