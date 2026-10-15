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

### Analytic Pipeline:
0) First, run 0_init.R to initialize the R environment after activating the immune.Rproj
1) then run SLURM jobs to clean population-level summary statistics from the GWAS_sumstats_sources google sheet
2) after cleaning, run METAL to combine biomarker measures (e.g., monocytes PanUKB + monocytes FinnGen)
3) after METAL, check columns and readiness for LDSC
4) run LDSC (pairwise against Loneliness isolation PanUKB) and generate forest plot
5) run all pairwise combinations and generate correlation heatmap

### Other Scripts
- 0_basophils_rg_sanity.R is used to check the genetic correlation between PanUKB and FinnGen basophil count (measured differently) to assess suitability of METAL
- 01_preprocess_sumstats.R is another sanity check used to get the column names in the different datasets, which then informs the way the batch/generic cleaning script works
- setup_ldsc.R contains cluster-specific information about the location of Python and LDSC packages also used to setup SLURM jobs using LDSC

