# METAL IVW + LDSC workflow

This is a separate three-job workflow for rows marked `METAL` and `Primary` in the Google Sheet containing sumstats for the Immune Biomarkers and Loneliness project.

- This folder contains the scripts and local registry needed by all three jobs (1: cleaning and combining specific sumstats with METAL, 2: cleaning and munging sumstats for LDSC, 3: LDSC with loneliness). 
- Update or replace that CSV in this folder when the pre-munged dataset registry changes; the jobs do not depend on `already_munged_ldsc`. 
- Also don't need to manually update the csv to include munged data that occurred through using this script. The script will automatically check that. 
- The cleaner accepts an optional analysis argument to select `Primary` rows and downloads linked public `.gz` inputs to scratch. Otherwise, it'll use local endpoints. The METAL jobs use raw summary statistics; LDSC-munged files do not have the effect/SE/allele detail needed for METAL. 
- For METAL only, PanUKB inputs are lifted from GRCh37 to GRCh38 before combining with FinnGen. The Primary/LDSC branch is unchanged and may reuse the exact-name pre-munged files regardless of build, as requested.

## Job order

1. Submit `submit_metal_ivw.sh`. It prepares all rows marked `METAL`, runs the seven specified IVW combinations, and saves `combined_manifest.tsv` plus combined trait files under `/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/metal_ivw/`.
2. Submit `submit_primary_clean.sh`. It processes only rows marked `Primary`, skips exact-match datasets in `already_munged.csv`, downloads linked files as needed, and writes the cleaning manifests under `/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/`.
3. After both jobs succeed, submit `submit_metal_primary_ldsc.sh`. It treats the seven METAL-combined traits as members of the Primary analysis set, reuses eligible exact-name pre-munged Primary files, munges remaining inputs, and compares every trait with `Loneliness isolation PanUKB` under the project results directory.

The first two jobs may run in either order or concurrently. The LDSC job must wait for both to finish successfully.

The Primary cleaner processes the coordinate-only GRCh37 GWAS Catalog `Lymphocyte Count` file in batches (default 250,000 rows) and caches completed rsID-lookup batches under `primary_clean/batch_cache/`.
Keep the cache on scratch between attempts; completed batches are reused, while an interrupted/incomplete output is rebuilt from cached batches. The source file is still streamed from the beginning to find batch boundaries.

## Outputs

- `metal_ivw/combined_manifest.tsv`: paths and sample sizes for the seven METAL outputs.
- `primary_clean/primary_manifest.tsv` and `clean_manifest.tsv`: Primary trait registry and cleaned inputs.
- `loneliness_isolation_rg_all_traits_fdr.tsv`: LDSC correlations for every Primary and METAL-combined trait against `Loneliness isolation PanUKB`, with BH-FDR values across this set of tests.
- `ldsc_sign_consistency.tsv`: allele-aware comparison of input BETA signs against LDSC-munged Z signs for newly munged traits. Pre-munged inputs are marked as not checked because this step does not have their raw effect/allele columns.

The METAL script maps coordinate-only FinnGen inputs to rsIDs using their GRCh38 coordinates. PanUKB files are streamed in batches (default 250,000 rows); each batch is mapped from GRCh37 coordinates to rsIDs and lifted to GRCh38 before being written. The finished batches are joined into the usual cleaned input file, so downstream METAL and LDSC manifests do not change. 
`METAL_CHUNK_ROWS` before submission to use a different batch size. 

### Caching for Large Sumstats
Completed PanUKB batches are cached under `metal_ivw/batch_cache/` on scratch and reused on reruns when the input file metadata, selected fields, sample size, and batch size match. The GRCh37-to-GRCh38 Ensembl chain is also cached there. 
Keep this scratch cache to resume; changing the input, sample size, or batch size creates a separate cache. 
The cache avoids repeating mapping and liftover for completed batches, although the source file is still streamed to find the batch boundaries.

### Sumstats Fields and Genome References for METAL Lifting
- The script recognizes PanUKB fields such as `beta_meta`, `se_meta`, and `pval_meta`. 
- If all three EUR-specific effect fields (`beta_EUR`, `se_EUR`, and `pval_EUR`) are present, it uses that matching EUR trio; it does not mix an EUR P value with meta-analysis beta and SE. 
- For `ref`/`alt` inputs, ALT is treated as the effect allele and REF as the other allele. 
- Coordinate lookup and liftover require MungeSumstats, the GRCh37/GRCh38 SNPlocs and BSgenome reference packages, and access to the Ensembl chain file on the cluster. 
- The job checks for these packages and prints their versions before processing. The script writes `metal_ivw/rsid_overlap_qc.tsv` with shared-ID counts and overlap fractions. METAL aligns effect direction using the alleles.
- For PanUKB only, when an rsID occurs more than once, all rows for that duplicated ID are excluded and the log reports the number and examples; this avoids arbitrarily selecting one multiallelic/ambiguous effect.
- Coordinate-only inputs are mapped against the source build’s reference (GRCh38 for FinnGen; GRCh37 for PanUKB) before the shared-rsID comparison.

### LDSC and Final Sample Size
- The combined LDSC sample size is set to the sum of the two cohort sample sizes, so this assumes the cohorts do not overlap - this was checked when the sumstats sheet was prepared and is why we're mostly combining PanUKB and FinnGen (which are independent cohorts). 
- The LDSC job runs one exposure-versus-trait comparison at a time and does not calculate correlations among the other traits.
- Might extend/add a new job that lets us compute the entire correlation matrix and build a heatmap. 
