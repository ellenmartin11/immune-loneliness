library(data.table)

source("/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/functions/cleaning_functions.R") # loads functions for cleaning

# location of raw sumstats
raw_sumstats_dir <- "/vast/palmer/scratch/polimanti/ejm96/raw_sumstats" 

# csv containing phenotype names and associated locations in the Yale cluster
sumstats_file_mapping <- "/gpfs/gibbs/project/polimanti/ejm96/immune/data/GWAS_sumstats_sources.csv"

# clean output directory
out_dir <- "/vast/palmer/scratch/polimanti/ejm96/clean_sumstats"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##################################################
# load mapping csv and select primary phenotypes #
##################################################

ref <- fread(sumstats_file_mapping)

paths <- setNames(
  ref$`Download Location`,
  ref$Phenotype
)

# specify primary phenotypes from finngen (refer to csv sheet for accurate phenotype names)
finngen_primary <- c(
  'Allergic rhinitis',
  'Allergy Asthma',
  'Depression FinnGen',
  'BMI FinnGen',
  'Covid-19', 
  'Crohn disease FinnGen',
  'Hypothyroidism FinnGen',
  'Rheumatoid Arthritis',
  'Sleep Disorders',
  'Type 2 Diabetes FinnGen',
  'Current Smoking',
  'Ulcerative colitis FinnGen'
)

##########################
# load and clean columns #
##########################

for (phenotype in finngen_primary) {
  
  message("Processing: ", phenotype)
  
  input_path <- paths[[phenotype]]
  
  if (!file.exists(input_path)) {
    warning("Missing file: ", input_path)
    next
  }
  
  # Load only one large file at a time
  dt <- fread(input_path, showProgress = TRUE)
  
  # Rename and standardize columns using custom function
  dt_clean <- clean_finngen(dt)
  
  # Create a safe output filename
  output_name <- gsub("[^A-Za-z0-9]+", "_", phenotype)
  output_path <- file.path(
    out_dir,
    paste0(output_name, "_clean.tsv.gz")
  )
  
  fwrite(
    dt_clean,
    output_path,
    sep = "\t"
  )
  
  # Explicitly release memory before the next phenotype
  rm(dt, dt_clean)
  gc()
  
  message("Saved: ", output_path)
}
