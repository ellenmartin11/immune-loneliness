#2_check_sumstats.R

# this script is used to check how well sumstats were prepared and preprocessed for LDSC
# this is a sanity check prior to running LDSC via slurm. If this all looks good, proceed to cluster scripts and run submit_ldsc_exposures.sh in generic_ldsc_cluster

cleaned_path <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/"

# to do: list all munged datasets in the munged_path that match the phenotype in Primary_Population_Sumstats tab
# loop through and print column names, should all return TRUE if matches expected LdSC input
# if not, print FALSE sumstats
# veridy that all phenotyupes are there. 

# ADHD
ADHD <- fread(paste0(cleaned_path, "001_ADHD_PGC_clean.tsv.gz"))
colnames(ADHD)

# Allergic Rhinitis
AllergicRhinitis <- fread(paste0(cleaned_path, "002_Allergic_rhinitis_FinnGen_clean.tsv.gz"))
colnames(AllergicRhinitis)

# Anorexia
Anorexia <- fread(paste0(cleaned_path, "003_Anorexia_Nervosa_PGC_clean.tsv.gz"))
colnames(Anorexia)

# Asthma
Asthma <- fread(paste0(cleaned_path, "004_Asthma_Multi_clean.tsv.gz"))
colnames(Asthma)
