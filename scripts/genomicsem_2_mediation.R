############################################################################
### genomicsem_2_mediation.R
###
### GenomicSEM mediation models:
###   LONELY -> M -> P   controlling for INCOME + BMI
###
###   M ~ a*LONELY + INCOME + BMI
###   P ~ b*M + c*LONELY + INCOME + BMI
###   ind := a*b          (indirect effect, delta-method SE)
###
### Source this after covstruc is loaded (same session as the adjustment
### models). Requires GenomicSEM's usermodel(), which computes SEs for
### ":=" defined parameters via the delta method.
############################################################################

## Fit one mediation model and return a one-row summary plus full results.
## mediator: short name in covstruc (e.g. "TG", "NEUT", "CRP")
## outcome:  short name in covstruc (e.g. "DEP", "PTSD", "ASTHMA", "T2D")
## exposure: "LONELY" (primary, PanUKB) or "LONELY_AOU" (sensitivity, AoU)
fit_med <- function(mediator, outcome, covstruc = NULL,
                    covariates = c("INCOME", "BMI"), estimation = "DWLS",
                    exposure = "LONELY", std.lv = FALSE) {
  if (is.null(covstruc)) covstruc <- get0("covstruc", envir = globalenv())
  if (is.null(covstruc)) {
    stop("covstruc not found - readRDS it into your global environment or pass it")
  }
  need <- c(exposure, mediator, outcome, covariates)
  miss <- need[!need %in% colnames(covstruc$S)]
  if (length(miss)) {
    stop("not in covstruc S matrix: ", paste(miss, collapse = ", "))
  }
  
  cov_str <- paste0(" + ", paste(covariates, collapse = " + "))
  model <- sprintf("%s ~ a*%s%s\n%s ~ b*%s + c*%s%s\nind := a*b\n",
                   mediator, exposure, cov_str,
                   outcome, mediator, exposure, cov_str)
  cat("\n== mediation model ==\n", model, "\n")
  um <- usermodel(covstruc, estimation = estimation, model = model,
                  std.lv = std.lv)
  res <- um$results
  
  est_col <- if ("Unstand_Est" %in% names(res)) "Unstand_Est" else "Unstandardized_Estimate"
  se_col  <- if ("Unstand_SE"  %in% names(res)) "Unstand_SE"  else "SE"
  num <- function(x) suppressWarnings(as.numeric(x))
  
  grab <- function(lhs, op, rhs) {
    r <- res[res$lhs == lhs & res$op == op & res$rhs == rhs, , drop = FALSE]
    if (nrow(r) != 1) {
      return(setNames(rep(NA_real_, 5), c("est", "se", "std", "std_se", "p")))
    }
    c(est = num(r[[est_col]][1]), se = num(r[[se_col]][1]),
      std = if ("STD_Genotype" %in% names(r)) num(r$STD_Genotype[1]) else NA_real_,
      std_se = if ("STD_Genotype_SE" %in% names(r)) num(r$STD_Genotype_SE[1]) else NA_real_,
      p = if ("p_value" %in% names(r)) suppressWarnings(as.numeric(r$p_value[1])) else NA_real_)
  }
  
  a   <- grab(mediator, "~", exposure)   # exposure -> M
  b   <- grab(outcome,  "~", mediator)   # M -> P (adjusted)
  cp  <- grab(outcome,  "~", exposure)   # direct effect
  ## indirect effect: defined parameter row (op ":=", rhs may be "" or NA)
  ind_r <- res[res$op == ":=", , drop = FALSE]
  ind <- if (nrow(ind_r) >= 1) {
    c(est = num(ind_r[[est_col]][1]), se = num(ind_r[[se_col]][1]),
      std = if ("STD_Genotype" %in% names(ind_r)) num(ind_r$STD_Genotype[1]) else NA_real_,
      std_se = if ("STD_Genotype_SE" %in% names(ind_r)) num(ind_r$STD_Genotype_SE[1]) else NA_real_,
      p = if ("p_value" %in% names(ind_r)) suppressWarnings(as.numeric(ind_r$p_value[1])) else NA_real_)
  } else {
    setNames(rep(NA_real_, 5), c("est", "se", "std", "std_se", "p"))
  }
  
  summary_row <- data.frame(
    exposure = exposure, mediator = mediator, outcome = outcome,
    a_est = a[["est"]],   a_se = a[["se"]],   a_p = a[["p"]],
    b_est = b[["est"]],   b_se = b[["se"]],   b_p = b[["p"]],
    direct_est = cp[["est"]], direct_se = cp[["se"]], direct_p = cp[["p"]],
    direct_std = cp[["std"]], direct_std_se = cp[["std_se"]],
    ind_est = ind[["est"]],   ind_se = ind[["se"]],   ind_p = ind[["p"]],
    total_est = as.numeric(cp[["est"]]) + as.numeric(ind[["est"]]),
    stringsAsFactors = FALSE
  )
  
  ## flag implausible SEs (ill-conditioned V)
  if (any(!is.finite(c(summary_row$a_se, summary_row$b_se,
                       summary_row$direct_se, summary_row$ind_se))) ||
      any(c(summary_row$a_se, summary_row$b_se, summary_row$direct_se,
            summary_row$ind_se) > 1e3, na.rm = TRUE)) {
    warning(sprintf("Implausible SE in mediation %s -> %s - model likely ill-conditioned; inspect $results",
                    mediator, outcome))
  }
  
  list(mediator = mediator, outcome = outcome, summary = summary_row,
       modelfit = um$modelfit, results = res)
}

## Assemble a table of mediation fits with BH-FDR across the indirect p-values
collect_mediation <- function(fits) {
  tab <- do.call(rbind, lapply(fits, function(f) f$summary))
  tab$ind_p_FDR_BH <- p.adjust(tab$ind_p, method = "BH")
  tab
}
