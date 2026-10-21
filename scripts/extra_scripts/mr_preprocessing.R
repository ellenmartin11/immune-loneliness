# harmonize for MR
harmonize_mr <- function(
    exposure_dat,
    outcome_dat,
    exposure_name,
    outcome_name,
    p_threshold = 5e-8,
    action = 2
) {
  
  # Filter exposure to requested p-value threshold
  exposure_dat <- exposure_dat[
    exposure_dat$pval < p_threshold &
      !is.na(exposure_dat$pval),
  ]
  
  # Convert data.table to data.frame for TwoSampleMR
  exposure_df <- as.data.frame(exposure_dat)
  outcome_df  <- as.data.frame(outcome_dat)
  
  # Format exposure
  exposure_formatted <- format_data(
    exposure_df,
    type = "exposure",
    snps = exposure_df$SNP,
    phenotype_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "EA",
    other_allele_col = "OA",
    pval_col = "pval"
  )
  
  # Format outcome
  outcome_formatted <- format_data(
    outcome_df,
    type = "outcome",
    snps = outcome_df$SNP,
    phenotype_col = "outcome",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "EA",
    other_allele_col = "OA",
    pval_col = "pval"
  )
  
  # Set phenotype names/IDs
  exposure_formatted$id.exposure <- exposure_name
  exposure_formatted$exposure <- exposure_name
  
  outcome_formatted$id.outcome <- outcome_name
  outcome_formatted$outcome <- outcome_name
  
  # Harmonize
  harmonized <- harmonise_data(
    exposure_dat = exposure_formatted,
    outcome_dat = outcome_formatted,
    action = action
  )
  
  # Add analysis information
  attr(harmonized, "p_threshold") <- p_threshold
  attr(harmonized, "exposure_name") <- exposure_name
  attr(harmonized, "outcome_name") <- outcome_name
  
  return(harmonized)
}

# clumping for MR
clump_mr <- function( 
    mr_dat, 
    clump_kb = 10000, 
    clump_r2 = 0.001, 
    clump_p = 1, 
    pop = "EUR" ) { 
  # Clump instruments using LD reference panel 
  clumped <- clump_data( mr_dat, 
                         clump_kb = clump_kb, 
                         clump_r2 = clump_r2, 
                         clump_p = clump_p, 
                         pop = pop ) 
  return(clumped) }


# MR pipeline
run_mr_pipeline <- function(
    clumped_dat,
    analysis_name,
    methods = c(
      "mr_egger_regression",
      "mr_ivw",
      "mr_weighted_mode",
      "mr_weighted_median",
      "mr_egger_regression_bootstrap"
    )
) {
  
  # -----------------------------
  # 1. Main MR analysis
  # -----------------------------
  
  mr_results <- mr(
    clumped_dat,
    method_list = methods
  )
  
  
  # -----------------------------
  # 2. Heterogeneity
  # -----------------------------
  
  heterogeneity <- mr_heterogeneity(
    clumped_dat
  )
  
  
  # -----------------------------
  # 3. Horizontal pleiotropy
  # -----------------------------
  
  pleiotropy <- mr_pleiotropy_test(
    clumped_dat
  )
  
  
  # -----------------------------
  # 4. Leave-one-out analysis
  # -----------------------------
  
  leaveoneout_ivw <- mr_leaveoneout(
    clumped_dat
  )
  
  leaveoneout_egger <- mr_leaveoneout(
    clumped_dat,
    method = mr_egger_regression
  )
  
  leaveoneout_egger_bootstrap <- mr_leaveoneout(
    clumped_dat,
    method = mr_egger_regression_bootstrap
  )
  
  
  # -----------------------------
  # 5. Store results
  # -----------------------------
  
  results <- list(
    analysis_name = analysis_name,
    mr = mr_results,
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy,
    leaveoneout_ivw = leaveoneout_ivw,
    leaveoneout_egger = leaveoneout_egger,
    leaveoneout_egger_bootstrap = leaveoneout_egger_bootstrap
  )
  
  return(results)
}

