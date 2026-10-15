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
