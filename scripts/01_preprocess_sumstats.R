# 1_preprocess_sumstats.R is used to check the formatting of downloaded summary statisics, which is used to design the SLURM batch cleaning scripts
# Part of sanity process, but does not actually need to be run for analyses to work. 
# Author: Ellen Martin



# location of raw downloaded sumstats in scratch
raw_sumstats_dir <- "/vast/palmer/scratch/polimanti/ejm96/raw_sumstats" 

# finngen directory
# finngen_dir <- "/gpfs/gibbs/pi/polimanti/FinnGen_R12/sumstats/"

# location for processed sumstats files
out_sumstats_dir <- "/gpfs/gibbs/project/polimanti/ejm96/immune/data/1_preprocess_sumstats"

# csv containing phenotype names and associated locations in the Yale cluster
# sumstats_file_mapping <- "/gpfs/gibbs/project/polimanti/ejm96/immune/data/GWAS_sumstats_sources.csv"

# google drive link to csv
sumstats_file_mapping <- "https://docs.google.com/spreadsheets/d/1tJhUVb1uH2Eku4P_pMAYwJCagoFiO_nYxyqmdHUgAFk/export?format=csv"

library(data.table)

ref <- fread(sumstats_file_mapping)

paths <- setNames(
  ref$`file location`,
  ref$phenotype
)

############################
# Checking Albumin FinnGen #
############################
library(dplyr)
albumin_fg <- ref %>%
  filter(phenotype == "Albumin FinnGen")

albumin_panukb <- ref %>%
  filter(phenotype == "Albumin PanUKB")

albumin_data <- fread(albumin_panukb$`file location`)

head(albumin_data)
colnames(albumin_data)

#################################
# Checking GSCAN Substance Data #
#################################

cigs_gscan <- ref %>%
  filter(phenotype == "Cigarettes per Day GSCAN")

cigs_data <- fread(cigs_gscan$`file location`)

head(cigs_data)
colnames(cigs_data)

#################################
# Checking Lymphocytes Data #####
#################################
lymph <- ref %>%
  filter(phenotype == "Lymphocyte Count")

lymph_data <- fread(lymph$`file location`)

head(lymph_data)
colnames(lymph_data)

###########################
# CLEANING PGC PHENOTYPES #
###########################

library(dplyr)
pgc_phenotypes <- ref %>%
  filter(cohort== "PGC")

print(pgc_phenotypes$phenotype)

pgc_primary <- c(
  'ADHD PGC',
  'Anorexia Nervosa PGC',
  'Anxiety Disorder PGC',
  'Bipolar PGC',
  "Depression PGC",
  "OCD",
  "PTSD"
  )

# check paths
for (phenotype in pgc_primary) {
  print(paths[[phenotype]])
}

# check column names
for (phenotype in pgc_primary) {
  path <- paths[[phenotype]]
  header <- fread(path, nrows = 0)
  
  cat("\n", phenotype, "\n", sep = "")
  print(names(header))
}

# checking Anorexia, Depression and PTSD

ANX_PGC <- fread("/vast/palmer/scratch/polimanti/ejm96/raw_sumstats/pgcAN2.2019-07.vcf.tsv.gz")

head(ANX_PGC)


############################
# CHECKING PHENOTYPES #
############################
library(dplyr)

check_phenotypes <- ref
check_phenotypes$phenotype

income <- check_phenotypes %>%
  filter(phenotype=="Income")

income_path <- income$`file location`

# checking one example with allergic rhinitis
income_raw <- data.table::fread(income_path)
colnames(income_raw)




#################################
# CLEANING ALL OTHER PHENOTYPES #
#################################

other_phenotypes <- ref %>%
  filter(!cohort %in% c("UKB2", "PGC"))



other_primary <- c(
  'Allergic rhintis',
  'Asthma',
  'C-reactive protein',
  'Crohn disease',
  'Glycoprotein acetylation', 
  'Lymphocyte Count', 
  #'Obesity', #still need to download
  #'Platelet Count', #still need to download
  'White blood cell count'
)

# loop to print file locations 
for (phenotype in other_primary) {
  print(paths[[phenotype]])
}

# loop to check column names
for (phenotype in other_primary) {
  path <- paths[[phenotype]]
  header <- fread(path, nrows = 0)
  
  cat("\n", phenotype, "\n", sep = "")
  print(names(header))
}


#####################
# UKBPPP Phenotypes #
#####################
library(dplyr)
ukbppp <- ref %>%
  filter(cohort== "UKBPPP")

ukbppp_primary <- ukbppp$phenotype

# loop to print file locations 
for (phenotype in ukbppp_primary) {
  print(paths[[phenotype]])
}

# extracting tar file gives each chromosomes .gz sumstats data
out_dir <- "ukb_ppp_extracted"
dir.create(out_dir, showWarnings = FALSE)

untar("/vast/palmer/scratch/polimanti/ejm96/raw_sumstats/TNF_P01375_OID20473_v1_Inflammation.tar", exdir = out_dir)
list.files(out_dir, recursive = TRUE, full.names = TRUE)

# inspect just one
library(data.table)

tnf_dir <- "ukb_ppp_extracted/TNF_P01375_OID20473_v1_Inflammation"

tnf_files <- list.files(
  tnf_dir,
  pattern = "\\.gz$",
  full.names = TRUE
)

# Inspect column names
fread(tnf_files[1], nrows = 0)

# Inspect first five variants
fread(tnf_files[1], nrows = 5)

