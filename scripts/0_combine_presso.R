# 0_combine_presso.R

# Usage: This script takes the output .rds objects from the cluster job (submit_mrpresso.sh) for the 12 comparison pairs. 
# Results are combined into csv showing outlier-corrected MR estimate, distortion test results and outlier test results. 


library(data.table)

presso_dir <- "/vast/palmer/pi/polimanti/ellen_m/immune/results/two_sample_mr/presso"

rds_files <- list.files(
  presso_dir,
  pattern = "^mrpresso_.*\\.rds$",
  full.names = TRUE
)

if (!length(rds_files)) {
  stop("No mrpresso_*.rds files found in: ", presso_dir)
}

get_presso_output <- function(x) {
  # Handle either the direct MR-PRESSO result or a wrapper containing it.
  if (!is.null(x[["MR-PRESSO results"]])) return(x)
  
  for (field in c("result", "presso_result", "mr_presso", "presso")) {
    if (is.list(x) && !is.null(x[[field]])) {
      found <- get_presso_output(x[[field]])
      if (!is.null(found)) return(found)
    }
  }
  
  NULL
}

get_value <- function(x, field) {
  if (is.null(x) || is.null(x[[field]]) || !length(x[[field]])) {
    return(NA_character_)
  }
  as.character(x[[field]][1])
}

presso_summary_rows <- list()
outlier_detail_rows <- list()

for (path in rds_files) {
  object <- readRDS(path)
  result <- get_presso_output(object)
  
  if (is.null(result)) {
    warning("Could not find MR-PRESSO results in: ", basename(path))
    next
  }
  
  comparison <- basename(path)
  comparison <- sub("^mr_presso_", "", comparison)
  comparison <- sub("\\.rds$", "", comparison)
  comparison <- sub("_stringent$", "", comparison)
  comparison <- gsub("_to_", " \u2192 ", comparison)
  
  main_mr <- as.data.table(result[["Main MR results"]])
  presso <- result[["MR-PRESSO results"]]
  
  raw <- main_mr[`MR Analysis` == "Raw"]
  corrected <- main_mr[`MR Analysis` == "Outlier-corrected"]
  
  global <- presso[["Global Test"]]
  distortion <- presso[["Distortion Test"]]
  outliers <- presso[["Outlier Test"]]
  
  outlier_n <- 0L
  outlier_sig_n <- NA_integer_
  outlier_min_p <- NA_character_
  
  if (!is.null(outliers) && nrow(as.data.table(outliers))) {
    outliers <- as.data.table(as.data.frame(outliers))
    p_text <- as.character(outliers$Pvalue)
    p_bound <- suppressWarnings(as.numeric(sub("^<", "", p_text)))
    
    outlier_n <- nrow(outliers)
    outlier_sig_n <- sum(p_bound <= 0.05, na.rm = TRUE)
    
    if (any(is.finite(p_bound))) {
      outlier_min_p <- p_text[which.min(p_bound)]
    }
    
    # Do not include row names or SNP IDs in this file.
    outlier_detail_rows[[length(outlier_detail_rows) + 1L]] <- data.table(
      comparison = comparison,
      outlier_test_row = seq_len(nrow(outliers)),
      RSSobs = as.numeric(outliers$RSSobs),
      Pvalue = p_text
    )
  }
  
  presso_summary_rows[[length(presso_summary_rows) + 1L]] <- data.table(
    comparison = comparison,
    
    raw_beta = if (nrow(raw)) as.numeric(raw[["Causal Estimate"]][1]) else NA_real_,
    raw_se = if (nrow(raw)) as.numeric(raw[["Sd"]][1]) else NA_real_,
    raw_p = if (nrow(raw)) as.character(raw[["P-value"]][1]) else NA_character_,
    
    outlier_adjusted_beta = if (nrow(corrected)) {
      as.numeric(corrected[["Causal Estimate"]][1])
    } else NA_real_,
    outlier_adjusted_se = if (nrow(corrected)) {
      as.numeric(corrected[["Sd"]][1])
    } else NA_real_,
    outlier_adjusted_p = if (nrow(corrected)) {
      as.character(corrected[["P-value"]][1])
    } else NA_character_,
    
    global_test_RSSobs = get_value(global, "RSSobs"),
    global_test_p = get_value(global, "Pvalue"),
    
    distortion_coefficient = get_value(distortion, "Distortion Coefficient"),
    distortion_p = get_value(distortion, "Pvalue"),
    
    outlier_tests_run = outlier_n > 0L,
    outlier_tests_n = outlier_n,
    outlier_p_le_0_05_n = outlier_sig_n,
    smallest_outlier_p = outlier_min_p
  )
}

summary_csv <- file.path(presso_dir, "mr_presso_summary.csv")
fwrite(rbindlist(presso_summary_rows, fill = TRUE), summary_csv)

outlier_csv <- file.path(presso_dir, "mr_presso_outlier_test_statistics_no_snp_ids.csv")
outlier_details <- rbindlist(outlier_detail_rows, fill = TRUE)
fwrite(outlier_details, outlier_csv)

cat("Saved summary: ", summary_csv, "\n", sep = "")
cat("Saved outlier-test statistics without SNP IDs: ", outlier_csv, "\n", sep = "")