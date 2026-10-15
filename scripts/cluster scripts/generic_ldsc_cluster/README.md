README.md

### Cluster Scripts
- contains several different scripts to be used for preprocessing raw summary statistics listed in the GWAS_sumstats_sources google sheet: 
https://docs.google.com/spreadsheets/d/1tJhUVb1uH2Eku4P_pMAYwJCagoFiO_nYxyqmdHUgAFk/edit?usp=sharing


### Running on SLURM
- submit_clean_sumstats.sh will run clean_sumstats.R, which retrieves sumstats files from the corresponding location in the sheet, and for FinnGen R13 is designed to pull from the public URL 
- after cleaning sumstats, check all sumstats using the 2_check_sumstats.R script, which should return the names and locations of munged, LDSC-ready sumstats and their column names. 
- after sumstats look cleaned and ready for LDSC, submit_ldsc_exposures.sh will run several pairwise LDSC analyses using the PAN UKB isolation loneliness phenotype using prepare_and_run_ldsc.R and also run_ldsc_exposures.R
- after that, we can use output to plot a forest plot (to look at genetic correlations) or make a heatmap. 
- LDSC will also give us SNP-h2 and SNP-h2-z-score estimates. 


