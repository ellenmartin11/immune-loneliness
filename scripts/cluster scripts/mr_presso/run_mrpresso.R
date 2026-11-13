library(MRPRESSO)

# For submitting mrpresso simulations that take ages. Submit in terminal sbatch submit_mrpresso.sh tg_to_t2d_stringent
# And it will run for that comparison and save the output as an rds in the results/two_sample_mr folder for next steps
# Can repeat with multiple other comparisons with presso

project_dir <- "/gpfs/gibbs/project/polimanti/ejm96/immune"
input_file <- file.path(
  project_dir,
  "results/two_sample_mr/tsmr_pairs_stringent.rds"
)
results_dir <- file.path(project_dir, "results/two_sample_mr/presso")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop("Usage: Rscript run_mrpresso.R <pair_name>")
}

pair_name <- args[[1]]

# Keep the pair key safe for use in an output filename.
if (!grepl("^[A-Za-z0-9_-]+$", pair_name)) {
  stop("Pair name may only contain letters, numbers, underscores, and hyphens.")
}

run_presso_pair <- function(hdat, nb_distribution = 10000, seed = 123) {
  d <- as.data.frame(hdat)
  
  required <- c(
    "SNP", "beta.exposure", "beta.outcome",
    "se.exposure", "se.outcome"
  )
  missing <- setdiff(required, names(d))
  if (length(missing)) {
    stop("Harmonized data are missing: ", paste(missing, collapse = ", "))
  }
  
  # Use the same variants retained for the main MR, if this flag is present.
  if ("mr_keep" %in% names(d)) {
    d <- d[d$mr_keep %in% TRUE, ]
  }
  
  d <- d[
    is.finite(d$beta.exposure) &
      is.finite(d$beta.outcome) &
      is.finite(d$se.exposure) & d$se.exposure > 0 &
      is.finite(d$se.outcome) & d$se.outcome > 0,
  ]
  
  d$SNP <- as.character(d$SNP)
  d <- d[!duplicated(d$SNP), ]
  
  if (anyNA(d$SNP) || any(!nzchar(d$SNP))) {
    stop("Missing SNP IDs in the harmonized data.")
  }
  
  rownames(d) <- d$SNP
  
  cat("Usable SNPs:", nrow(d), "\n")
  
  if (nrow(d) < 4L) {
    stop(
      "MR-PRESSO needs at least 4 usable SNP instruments; found ",
      nrow(d)
    )
  }
  
  # NbDistribution must exceed the number of instruments.
  nb_distribution <- max(nb_distribution, nrow(d) + 1L)
  
  cat("MR-PRESSO simulations:", nb_distribution, "\n")
  
  MRPRESSO::mr_presso(
    BetaOutcome = "beta.outcome",
    BetaExposure = "beta.exposure",
    SdOutcome = "se.outcome",
    SdExposure = "se.exposure",
    data = d,
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    NbDistribution = nb_distribution,
    SignifThreshold = 0.05,
    seed = seed
  )
}

cat("Reading:", input_file, "\n")
tsmr_pairs_stringent <- readRDS(input_file)

pair_data <- tsmr_pairs_stringent[[pair_name]]
if (is.null(pair_data)) {
  stop("Pair not found in tsmr_pairs_stringent: ", pair_name)
}
if (is.null(pair_data$harmonised)) {
  stop("No harmonized data found for pair: ", pair_name)
}

cat("Running MR-PRESSO for:", pair_name, "\n")

presso_result <- run_presso_pair(
  pair_data$harmonised,
  nb_distribution = 10000,
  seed = 123
)

print(presso_result)

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(results_dir, paste0("mrpresso_", pair_name, ".rds"))

saveRDS(presso_result, output_file)
cat("Saved MR-PRESSO result to:", output_file, "\n")