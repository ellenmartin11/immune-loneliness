# METAL IVW + LDSC workflow

This is a separate three-job workflow for rows marked `METAL` and `Primary` in the Google Sheet.
Copy the files in this folder to:

`/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/cluster scripts/metal_ldsc_pipeline/`

This folder contains the scripts and local registry needed by all three jobs, including its own `clean_sumstats.R` and `already_munged.csv`. Update or replace that CSV in this folder when your pre-munged dataset registry changes; the jobs do not depend on `already_munged_ldsc`. The cleaner accepts an optional analysis argument to select `Primary` rows and downloads linked public `.gz` inputs to scratch. The METAL jobs use raw summary statistics; LDSC-munged files do not have the effect/SE/allele detail needed for METAL. For METAL only, PanUKB inputs are lifted from GRCh37 to GRCh38 before combining with FinnGen. The Primary/LDSC branch is unchanged and may reuse the exact-name pre-munged files regardless of build, as requested.

## Job order

1. Submit `submit_metal_ivw.sh`. It prepares all rows marked `METAL`, runs the seven specified IVW combinations, and saves `combined_manifest.tsv` plus combined trait files under `/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw/`.
2. Submit `submit_primary_clean.sh`. It processes only rows marked `Primary`, skips exact-match datasets in `already_munged.csv`, downloads linked files as needed, and writes the cleaning manifests under `/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/`.
3. After both jobs succeed, submit `submit_metal_primary_ldsc.sh`. It treats the seven METAL-combined traits as members of the Primary analysis set, reuses eligible exact-name pre-munged Primary files, munges remaining inputs, and compares every trait with `Loneliness isolation PanUKB` under the project results directory.

The first two jobs may run in either order or concurrently. The LDSC job must wait for both to finish successfully. METAL can be selected by setting `METAL_BIN` to its executable path before submission; the default is `METAL` on `PATH`.

## Outputs

- `metal_ivw/combined_manifest.tsv`: paths and sample sizes for the seven METAL outputs.
- `primary_clean/primary_manifest.tsv` and `clean_manifest.tsv`: Primary trait registry and cleaned inputs.
- `loneliness_isolation_rg_all_traits_fdr.tsv`: LDSC correlations for every Primary and METAL-combined trait against `Loneliness isolation PanUKB`, with BH-FDR values across this set of tests.
- `ldsc_sign_consistency.tsv`: allele-aware comparison of input BETA signs against LDSC-munged Z signs for newly munged traits. Pre-munged inputs are marked as not checked because this step does not have their raw effect/allele columns.

The METAL script maps coordinate-only FinnGen inputs to rsIDs using their GRCh38 coordinates. PanUKB files are first mapped from GRCh37 coordinates to rsIDs, then lifted to GRCh38 before combination. It recognizes PanUKB fields such as `beta_meta`, `se_meta`, and `pval_meta`. If all three EUR-specific effect fields (`beta_EUR`, `se_EUR`, and `pval_EUR`) are present, it uses that matching EUR trio; it does not mix an EUR P value with meta-analysis beta and SE. For `ref`/`alt` inputs, ALT is treated as the effect allele and REF as the other allele. Coordinate lookup and liftover require MungeSumstats, the GRCh37/GRCh38 SNPlocs and BSgenome reference packages, and access to the Ensembl chain file on the cluster. The job checks for these packages and prints their versions before processing. The script writes `metal_ivw/rsid_overlap_qc.tsv` with shared-ID counts and overlap fractions. METAL aligns effect direction using the alleles; duplicate rsIDs stop the job rather than being silently double-counted.

Coordinate-only inputs are mapped against the source build’s reference (GRCh38 for FinnGen; GRCh37 for PanUKB) before the shared-rsID comparison. The combined LDSC sample size is set to the sum of the two cohort sample sizes, so this assumes the cohorts do not overlap. The LDSC job runs one exposure-versus-trait comparison at a time and does not calculate correlations among the other traits.
