# cleaning / preprocessing functions

library(data.table)

### Clean FinnGen
clean_finngen <- function(dt, add_effect_alleles = TRUE) {
  
  if (is.character(dt) && length(dt) == 1 && file.exists(dt)) {
    dt <- data.table::fread(dt)
  }
  
  if(!data.table::is.data.table(dt)) {
    data.table::setDT(dt)
  }
  
  required <- c(
    "#chrom",
    "pos",
    "ref",
    "alt",
    "rsids",
    "pval",
    "beta",
    "sebeta"
  )
  
  # check for missing columns first
  missing_cols <- setdiff(required, names(dt))

  if (length(missing_cols) > 0) {
    stop(
      "Missing expected FinnGen columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # otherwise progress with renaming
  setnames(
    dt,
    old = c(
      "#chrom",
      "pos",
      "ref",
      "alt",
      "rsids",
      "pval",
      "beta",
      "sebeta"
    ),
    new = c(
      "CHR",
      "POS",
      "REF",
      "ALT",
      "SNP",
      "P",
      "BETA",
      "SE"
    ),
    skip_absent = TRUE
  )
  
  # renaming columns for logistic traits
  optional_renames <- list(
    nearest_genes = "NEAREST_GENE",
    mlogp = "NEG_LOG10_P",
    af_alt = "EAF",
    af_alt_cases = "EAF_CASES",
    af_alt_controls = "EAF_CONTROLS"
  )
  
  for (old_name in names(optional_renames)) {
    new_name <- optional_renames[[old_name]]
    
    if (old_name %in% names(dt)) {
      setnames(dt, old_name, new_name)
    }
  }
  
  # add optional columns when absent
  
  optional_columns <- c(
    "NEAREST_GENE",
    "NEG_LOG10_P",
    "EAF",
    "EAF_CASES",
    "EAF_CONTROLS"
  )
  
  for (column in optional_columns) {
    if (!(column %in% names(dt))) {
      dt[, (column) := NA_real_]
    }
  }
  
  
  if (add_effect_alleles) {
    dt[, EA := ALT]
    dt[, OA := REF]
  }
  
  # standardising datatypes
  dt[, CHR := as.character(CHR)]
  dt[, POS := as.integer(POS)]
  dt[, P := as.numeric(P)]
  dt[, BETA := as.numeric(BETA)]
  dt[, SE := as.numeric(SE)]
  dt[, EAF := as.numeric(EAF)]
  dt[, EAF_CASES := as.numeric(EAF_CASES)]
  dt[, EAF_CONTROLS := as.numeric(EAF_CONTROLS)]
  
  if ("NEG_LOG10_P" %in% names(dt)) {
    dt[, NEG_LOG10_P := as.numeric(NEG_LOG10_P)]
  }
  
  # order
  preferred_order <- c(
    "CHR", "POS", "SNP", "REF", "ALT", "EA", "OA",
    "BETA", "SE", "P", "NEG_LOG10_P", 
    "EAF", "EAF_CASES", "EAF_CONTROLS",
    "NEAREST_GENE"
  ) 
  
  remaining <- setdiff(names(dt), preferred_order)
  
  setcolorder(
    dt,
    c(intersect(preferred_order, names(dt)), remaining)
  )
  
  return(dt)
}



### Clean PanUKB
prepare_panukb <- function(input_path, output_path) {
  
  message("Preprocessing PanUKB European GWAS statistics")
  
  required_cols <- c(
    "chr", "pos", "ref", "alt",
    "beta_EUR", "se_EUR", "pval_EUR", "af_EUR"
  )
  
  dt <- fread(
    input_path,
    select = required_cols,
    showProgress = TRUE
  )
  
  setnames(
    dt,
    old = required_cols,
    new = c(
      "CHR", "BP", "A2", "A1",
      "BETA", "SE", "P", "FRQ"
    )
  )
  
  message("Original variants: ", nrow(dt))
  
  # Remove variants with missing or invalid essential statistics
  dt <- dt[
    !is.na(CHR) &
      !is.na(BP) &
      !is.na(A1) &
      !is.na(A2) &
      is.finite(BETA) &
      is.finite(SE) &
      SE > 0 &
      is.finite(P) &
      P >= 0 &
      P <= 1
  ]
  
  # Remove invalid allele frequencies
  dt <- dt[
    is.finite(FRQ) &
      FRQ >= 0 &
      FRQ <= 1
  ]
  
  # Restrict to autosomes for LDSC
  dt <- dt[
    CHR %in% as.character(1:22)
  ]
  
  message("Variants after preprocessing: ", nrow(dt))
  
  fwrite(
    dt,
    output_path,
    sep = "\t"
  )
  
  message("Clean PanUKB file saved: ", output_path)
  
  output_path
}


### Formatting for LDSC
format_for_ldsc <- function(path, label, source_build, output_path) {
  
  message(
    "Formatting ", label,
    " (input build: ", source_build, ")"
  )
  
  # Only perform liftover when necessary
  target_build <- if (source_build == "GRCh37") {
    "GRCh38"
  } else {
    NULL
  }
  
  result <- MungeSumstats::format_sumstats(
    path = path,
    ref_genome = source_build,
    convert_ref_genome = target_build,
    save_path = output_path,
    return_data = FALSE,
    force_new = TRUE,
    nThread = 2,
    compute_z = TRUE,
    snp_ids_are_rs_ids = has_rsid_column(path),
    on_ref_genome = TRUE,
    allele_flip_check = TRUE,
    allele_flip_drop = TRUE,
    log_mungesumstats_msgs = TRUE,
    log_folder = work_dir
  )
  
  if (!file.exists(output_path)) {
    stop(
      "MungeSumstats did not create output: ",
      output_path
    )
  }
  
  message("Formatted file saved: ", output_path)
  
  output_path
}

### LDSC Munging
run_ldsc_munge <- function(input_path, output_prefix, N_override) {
  header <- names(fread(input_path, nrows = 0L, showProgress = FALSE, check.names = FALSE))
  header_norm <- tolower(gsub("[^a-z0-9]", "", header))
  has_n_column <- any(header_norm %in% c("n", "samplesize", "ncas", "ncon", "ncase", "ncontrol"))
  if (!is.finite(N_override) && !has_n_column) {
    stop("No sample-size column was found in ", input_path,
         ". Enter the study's total sample size in the corresponding *_N setting at the top.")
  }
  cmd_args <- c(
    munge_py,
    "--sumstats", input_path,
    "--out", output_prefix,
    "--merge-alleles", hm3_snplist,
    "--ignore", "BETA" # ignores BETA because z-scores are present - having both will cause issues
  )
  
  if (is.finite(N_override)) cmd_args <- c(cmd_args, "--N", as.character(N_override))
  status <- system2(ldsc_python, args = cmd_args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status")
  if (is.null(exit_status)) exit_status <- 0L
  if (exit_status != 0L) {
    stop("LDSC munge failed for ", input_path, ":\n", paste(status, collapse = "\n"))
  }
  sumstats_path <- paste0(output_prefix, ".sumstats.gz")
  if (!file.exists(sumstats_path)) stop("Expected LDSC-munged file was not created: ", sumstats_path)
  sumstats_path
}

############################################
# CHECKING QUALITY OF PRIMARY_CLEAN SCRIPT #
############################################

qc_sumstats <- function(path) {
  
  message("Checking: ", path)
  
  dt <- fread(path)
  
  # check basic dim and headers
  
  n_rows <- nrow(dt)
  cols <- names(dt)
  
  missing_required <- setdiff(required_cols, cols)
  
  # check SNP id formatting
  
  n_missing_snp <- if ("SNP" %in% cols) sum(is.na(dt$SNP) | dt$SNP == "") else NA_integer_
  
  n_duplicate_snp <- if ("SNP" %in% cols) {
    sum(duplicated(dt$SNP))
  } else {
    NA_integer_
  }
  
  prop_rsids <- if ("SNP" %in% cols) {
    mean(grepl("^rs[0-9]+$", dt$SNP), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  # check valid alleles
  
  valid_alleles <- c("A", "C", "G", "T")
  
  invalid_A1 <- if ("A1" %in% cols) {
    sum(!dt$A1 %in% valid_alleles, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  invalid_A2 <- if ("A2" %in% cols) {
    sum(!dt$A2 %in% valid_alleles, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  same_allele <- if (all(c("A1", "A2") %in% cols)) {
    sum(dt$A1 == dt$A2, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  # check palindromic alleles
  palindromic <- if (all(c("A1", "A2") %in% cols)) {
    
    pair <- paste0(dt$A1, dt$A2)
    
    sum(pair %in% c(
      "AT", "TA",
      "CG", "GC"
    ), na.rm = TRUE)
    
  } else {
    NA_integer_
  }
  
  # chromosome sanity check
  
  invalid_chr <- if ("CHR" %in% cols) {
    
    chr <- suppressWarnings(as.integer(dt$CHR))
    
    sum(
      is.na(chr) |
        chr < 1 |
        chr > 22
    )
    
  } else {
    NA_integer_
  }
  
  # check positions of base pairs
  
  invalid_bp <- if ("BP" %in% cols) {
    
    bp <- suppressWarnings(as.numeric(dt$BP))
    
    sum(
      is.na(bp) |
        bp <= 0
    )
    
  } else {
    NA_integer_
  }
  
  # check effect sizes (BETA)
  
  missing_beta <- if ("BETA" %in% cols) {
    sum(is.na(dt$BETA) | !is.finite(dt$BETA))
  } else {
    NA_integer_
  }
  
  extreme_beta <- if ("BETA" %in% cols) {
    sum(abs(dt$BETA) > 10, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  # check standard errors
  
  invalid_se <- if ("SE" %in% cols) {
    sum(
      is.na(dt$SE) |
        !is.finite(dt$SE) |
        dt$SE <= 0
    )
  } else {
    NA_integer_
  }
  
  # check p-vals
  
  invalid_p <- if ("P" %in% cols) {
    sum(
      is.na(dt$P) |
        !is.finite(dt$P) |
        dt$P < 0 |
        dt$P > 1
    )
  } else {
    NA_integer_
  }
  
  zero_p <- if ("P" %in% cols) {
    sum(dt$P == 0, na.rm = TRUE)
  } else {
    NA_integer_
  }
  
  # check N
  
  invalid_n <- if ("N" %in% cols) {
    sum(
      is.na(dt$N) |
        !is.finite(dt$N) |
        dt$N <= 0
    )
  } else {
    NA_integer_
  }
  
  median_n <- if ("N" %in% cols) {
    median(dt$N, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  
  # Z-score consistency
  # Compare stored Z against BETA / SE if both available.
  # This is useful for catching sign flips / bad cleaning.
  # ----------------------------------------------------------
  
  z_correlation <- NA_real_
  
  if (all(c("BETA", "SE", "Z") %in% cols)) {
    
    z_calc <- dt$BETA / dt$SE
    
    valid <- (
      is.finite(z_calc) &
        is.finite(dt$Z)
    )
    
    if (sum(valid) > 100) {
      z_correlation <- cor(
        z_calc[valid],
        dt$Z[valid],
        use = "complete.obs"
      )
    }
  }
  
  # P-value consistency
  # Check whether P approximately agrees with BETA / SE.
  
  p_correlation <- NA_real_
  
  if (all(c("BETA", "SE", "P") %in% cols)) {
    
    z_calc <- dt$BETA / dt$SE
    
    p_calc <- 2 * pnorm(
      abs(z_calc),
      lower.tail = FALSE
    )
    
    valid <- (
      is.finite(p_calc) &
        is.finite(dt$P) &
        dt$P > 0 &
        p_calc > 0
    )
    
    if (sum(valid) > 100) {
      
      p_correlation <- cor(
        -log10(dt$P[valid]),
        -log10(p_calc[valid]),
        use = "complete.obs"
      )
    }
  }
  
  # Test-statistic distribution
  
  median_chisq <- NA_real_
  lambda_gc <- NA_real_
  
  if (all(c("BETA", "SE") %in% cols)) {
    
    z <- dt$BETA / dt$SE
    
    chisq <- z^2
    
    median_chisq <- median(
      chisq[is.finite(chisq)],
      na.rm = TRUE
    )
    
    # Median chi-square under H0 for 1 df ~= 0.4549364
    lambda_gc <- median_chisq / 0.4549364
  }
  
  # Summary
  
  data.table(
    file = path,
    
    n_rows = n_rows,
    n_columns = length(cols),
    
    missing_required_cols =
      ifelse(
        length(missing_required) == 0,
        "",
        paste(missing_required, collapse = ",")
      ),
    
    n_missing_snp = n_missing_snp,
    n_duplicate_snp = n_duplicate_snp,
    prop_rsids = prop_rsids,
    
    invalid_A1 = invalid_A1,
    invalid_A2 = invalid_A2,
    same_allele = same_allele,
    palindromic = palindromic,
    
    invalid_chr = invalid_chr,
    invalid_bp = invalid_bp,
    
    missing_beta = missing_beta,
    extreme_beta = extreme_beta,
    
    invalid_se = invalid_se,
    
    invalid_p = invalid_p,
    zero_p = zero_p,
    
    invalid_n = invalid_n,
    median_n = median_n,
    
    z_correlation = z_correlation,
    p_correlation = p_correlation,
    
    median_chisq = median_chisq,
    lambda_gc = lambda_gc
  )
}