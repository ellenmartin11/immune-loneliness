############################################################################
### genomicsem_1_interactive.R
###
### Interactive fitting of the univariate adjustment models:
###   LONELY ~ P_i + INCOME + BMI   for each psychiatric trait P_i
###
### Use AFTER genomicsem_1_ldsc_setup.R has produced LDSCoutput.rds.
###
### In an interactive R session (RStudio on the cluster):
###
###   setwd("/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/genomicsem")
###   library(GenomicSEM)
###   covstruc <- readRDS("LDSCoutput.rds")
###   source("genomicsem_2_interactive.R")
###
### Then fit models one at a time, e.g.:
###
###   f_DEP <- fit_one("DEP")                      # primary: LONELY (PanUKB)
###   f_DEP$path                      # inspect the DEP -> LONELY path
###   f_DEP$modelfit                  # chisq, p, df, AIC, CFI, SRMR
###   f_DEP$results                   # full parameter table
###
### Sensitivity analysis with the AoU UCLA loneliness GWAS:
###
###   s_DEP <- fit_one("DEP", outcome = "LONELY_AOU")
###
### When you have run the ones you want, collect them into one table:
###
###   fits <- list(ANX = f_ANX, DEP = f_DEP, PTSD = f_PTSD, ADHD = f_ADHD,
###                BIP = f_BIP, ANOREXIA = f_ANOREXIA, OCD = f_OCD, SCZ = f_SCZ)
###   tab <- collect_results(fits)
###   write.csv(tab, "lonely_adjusted_paths.csv", row.names = FALSE)
###
###   fits_aou <- lapply(fits, function(f) fit_one(f$trait, outcome = "LONELY_AOU"))
###   tab_aou <- collect_results(fits_aou)
###   write.csv(tab_aou, "lonely_aou_adjusted_paths.csv", row.names = FALSE)
############################################################################

PSYCH_TRAITS <- c("ANX", "DEP", "PTSD", "ADHD", "BIP", "ANOREXIA", "OCD", "SCZ")
COVARIATES   <- c("INCOME", "BMI")
OUTCOME      <- "LONELY"          # primary (PanUKB loneliness isolation)
OUTCOME_AOU  <- "LONELY_AOU"      # sensitivity (AoU UCLA loneliness scale)

## Fit one model: OUTCOME ~ P + INCOME + BMI
## covstruc: the ldsc() output (list with V, S, I, N, m). If omitted, the
## global `covstruc` is used (set by readRDS in your session).
fit_one <- function(p, covstruc = NULL,
                    outcome = OUTCOME, covariates = COVARIATES,
                    estimation = "DWLS", std.lv = FALSE) {
  if (is.null(covstruc)) {
    covstruc <- get0("covstruc", envir = globalenv())
    if (is.null(covstruc)) {
      stop("covstruc not found - either readRDS it into your global environment or pass it: fit_one(p, covstruc)")
    }
  }
  stopifnot(outcome %in% colnames(covstruc$S))
  model <- paste(outcome, "~", paste(c(p, covariates), collapse = " + "))
  cat("\n== model:", model, "==\n")
  um <- usermodel(covstruc, estimation = estimation, model = model, std.lv = std.lv)
  
  res <- um$results
  ## Defensive: converged branch has Unstand_Est/Unstand_SE/STD_Genotype/STD_All;
  ## error branch has Unstandardized_Estimate/SE only.
  est_col <- if ("Unstand_Est" %in% names(res)) "Unstand_Est" else "Unstandardized_Estimate"
  se_col  <- if ("Unstand_SE"  %in% names(res)) "Unstand_SE"  else "SE"
  std_col <- if ("STD_Genotype" %in% names(res)) "STD_Genotype" else NA_character_
  
  keep <- c("lhs", "op", "rhs", est_col, se_col,
            if (!is.na(std_col)) c(std_col, "STD_Genotype_SE"),
            if (!is.na(std_col)) "STD_All",
            if ("p_value" %in% names(res)) "p_value")
  
  path <- res[res$lhs == outcome & res$op == "~" & res$rhs == p, keep, drop = FALSE]
  rownames(path) <- NULL
  
  ## Flag implausible SEs (ill-conditioned V / non-convergence artifacts)
  se_num <- suppressWarnings(as.numeric(path[[se_col]]))
  if (any(!is.finite(se_num)) || any(se_num > 1e3, na.rm = TRUE)) {
    warning(sprintf(
      "SE for the %s ~ %s path is non-finite or > 1e3 - the sampling covariance matrix is likely ill-conditioned for this model. Inspect f$results and consider whether the trait set is too low-powered.", outcome, p))
  }
  
  list(trait = p, outcome = outcome, model = model, path = path,
       modelfit = um$modelfit, results = res)
}

## Collect fitted models into one table with BH (FDR) adjustment across paths
collect_results <- function(fits) {
  rows <- lapply(fits, function(f) {
    x <- f$path
    x$trait   <- f$trait
    x$outcome <- f$outcome
    x$model   <- f$model
    x
  })
  tab <- do.call(rbind, rows)
  if ("p_value" %in% names(tab)) {
    tab$p_FDR_BH <- p.adjust(tab$p_value, method = "BH")
  }
  ## order rows as PSYCH_TRAITS where possible
  tab <- tab[order(match(tab$trait, PSYCH_TRAITS)), , drop = FALSE]
  rownames(tab) <- NULL
  tab
}
