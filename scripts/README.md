# Immune Conditions/Biomarkers and Social Isolation

Date Edited: 06/10/2026
Written By: Ellen Martin (e.martin@yale.edu)

### File Structure and Directories:
- AoU: contains scripts and documentation relevant to AoU processing
- cluster scripts: contain scripts specifically for submission to Yale's SLURM cluster
    - finngen_cluster_scripts: scripts for cleaning and running LDSC on R12 FinnGen traits (not used for main project)
    - already_munged_ldsc: scripts for running LDSC on already munged traits (used as sanity check and expected genetic correlations)
    - generic_ldsc_cluster: scripts for cleaning sumstats from various sources and running LDSC against loneliness exposure (buggy - haven't fixed)
    - metal_ldsc_pipeline: scripts for preparing traits for METAL, cleaning sumstats, and combining for LDSC analysis
- functions: various custom R functions which may be used in cluster scripts/general sanity scripts
- R: contains environment set-up for work within the R interactive app on Yale's HPC 
- extra_scripts: contains other scripts from other projects that may be useful (for implementing MR, LDSC heatmaps etc)

### Analytic Pipeline:
0) First, run 0_init.R to initialize the R environment after activating the immune.Rproj
1a) then run SLURM jobs to clean population-level summary statistics from the GWAS_sumstats_sources google sheet
1b) Simultaneously (if you wan) run METAL to combine biomarker measures (e.g., monocytes PanUKB + monocytes FinnGen)
2a) check quality of SLURM cleaning using 2_check_sumstats_deep.R, which should return which sumstats need double-checking, and which are good to go!
2b) Check quality of METAL analysis and readiness for LDSC
3) build munge manifest once sumstats look good (gets correct Ns) - for both METAL and PRIMARY traits
### Not yet done
3) Munge all sumstats (including new combined METAL phenotypes) for LDSC (in parallel and in batches) - TO WRITE SCRIPT
4) run LDSC (pairwise against Loneliness isolation PanUKB) and generate forest plot - TO WRITE SCRIPT
5) run all pairwise combinations and generate correlation heatmap (in parallel and in batches without repeating pairs) - TO WRITE SCRIPT

### Other Scripts
- 0_basophils_rg_sanity.R is used to check the genetic correlation between PanUKB and FinnGen basophil count (measured differently) to assess suitability of METAL
- 01_preprocess_sumstats.R is another sanity check used to get the column names in the different datasets, which then informs the way the batch/generic cleaning script works

### Script Info

2_check_sumstats_deep.R

- quantifies: 
    - (i) variant ID integrity, including the proportion of rsID-formatted identifiers, duplicate rsIDs into position conflicts, multi-allelic sites, and identical rows; 
    - (ii) allele integrity (proportion of valid ACGT allele pairs; A1 = A2 rows); 
    - (iii) allele frequency validity (proportion outside (0,1)); 
    - (iv) internal consistency of signed statistics, testing Z = BETA/SE (≥ 95% agreement required) and the correlation between |Z| and the inverse-normal-transformed P-value (≥ 0.95 required); 
    - (v) P-value range validity; 
    - (vi) coordinate plausibility against GRCh37 chromosome lengths to detect genome-build mismatches;
    - (vii) sample-size distribution, including effective sample size Neff = 4/(1/N_case + 1/N_control) for case–control phenotypes.

Issues are classified as action-level (blocking munging until resolved) or informational (expected for whole-genome sources and handled by munging); phenotypes with no action-level flags were registered as QC-passed

