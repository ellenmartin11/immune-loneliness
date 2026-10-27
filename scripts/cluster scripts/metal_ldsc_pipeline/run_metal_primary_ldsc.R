#!/usr/bin/env Rscript

# Run either the input-preparation stage or the resumable LDSC correlation stage.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L || !args[[1L]] %chin% c("prepare", "rg")) {
  stop(paste("Usage: run_metal_primary_ldsc.R <prepare|rg> <combined_manifest> <clean_manifest> <primary_manifest>",
             "<already_munged.csv> <munged_dir> <results_dir> <python> <munge_sumstats.py>",
             "<ldsc.py> <hm3_snplist> <ref_ld_chr_prefix> <w_ld_chr_prefix>"))
}
stage <- args[[1L]]
args <- args[-1L]
combined_manifest_path <- args[[1L]]
clean_manifest_path <- args[[2L]]
primary_manifest_path <- args[[3L]]
already_munged_path <- args[[4L]]
munged_dir <- args[[5L]]
results_dir <- args[[6L]]
python <- args[[7L]]
munge_py <- args[[8L]]
ldsc_py <- args[[9L]]
hm3_snplist <- args[[10L]]
ref_ld_chr <- args[[11L]]
w_ld_chr <- args[[12L]]
dir.create(munged_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(combined_manifest_path, clean_manifest_path, primary_manifest_path,
               already_munged_path, python, munge_py, ldsc_py, hm3_snplist,
               paste0(ref_ld_chr, "1.l2.ldscore.gz"), paste0(w_ld_chr, "1.l2.ldscore.gz"))) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
}
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
find_col <- function(nms, candidate) {
  hit <- which(norm(nms) == norm(candidate))
  if (length(hit)) nms[[hit[[1L]]]] else NULL
}
safe_label <- function(x) {
  vapply(as.character(x), function(value) {
    if (is.na(value)) value <- ""
    value <- gsub("[^A-Za-z0-9]+", "_", value)
    value <- gsub("^_+|_+$", "", value)
    if (!nzchar(value)) "trait" else value
  }, character(1), USE.NAMES = FALSE)
}

# The correlation stage consumes the stable manifest created by stage one.
# Existing valid pair logs are reused, so a restarted job only runs unfinished pairs.
if (stage == "rg") {
  input_manifest_path <- file.path(munged_dir, "ldsc_inputs_manifest.tsv")
  if (!file.exists(input_manifest_path)) {
    stop("Prepared LDSC-input manifest not found: ", input_manifest_path,
         ". Complete submit_metal_primary_munge.sh first.")
  }
  traits <- fread(input_manifest_path, sep = "\t", showProgress = FALSE)
  if (!all(c("phenotype", "source", "origin", "path") %in% names(traits))) {
    stop("Prepared LDSC-input manifest is missing required fields.")
  }
  if (!any(traits$phenotype == "Loneliness isolation PanUKB")) {
    stop("Loneliness isolation PanUKB is not available in the prepared LDSC inputs.")
  }
  for (input_path in traits$path) {
    if (!file.exists(input_path) || file.access(input_path, 4L) != 0L) {
      stop("Prepared LDSC input is missing or unreadable: ", input_path)
    }
  }

  sign_check_path <- file.path(munged_dir, "ldsc_sign_consistency.tsv")
  if (file.exists(sign_check_path)) {
    file.copy(sign_check_path, file.path(results_dir, "ldsc_sign_consistency.tsv"), overwrite = TRUE)
  }

  extract_rg_result <- function(log_path) {
    if (!file.exists(log_path)) stop("LDSC did not create its log: ", log_path)
    lines <- sub("^.*?INFO\\s+-\\s+", "", readLines(log_path, warn = FALSE), perl = TRUE)
    marker <- grep("Summary of Genetic Correlation Results", lines, fixed = TRUE)
    if (!length(marker)) stop("No genetic-correlation table found in: ", log_path)
    block <- lines[seq.int(marker[[1L]] + 1L, length(lines))]
    words <- strsplit(trimws(block), "\\s+")
    header <- which(vapply(words, function(x) all(c("p1", "p2", "rg", "se", "z", "p") %in% tolower(x)), logical(1)))
    if (!length(header)) stop("Could not identify the LDSC result-table header in: ", log_path)
    h <- header[[1L]]
    rows <- block[seq.int(h + 1L, length(block))]
    row_words <- strsplit(trimws(rows), "\\s+")
    keep <- vapply(row_words, function(x) length(x) >= 6L && grepl("^[-+0-9.eE]+$|^NA$|^nan$", x[[3L]], ignore.case = TRUE), logical(1))
    if (!any(keep)) stop("No correlation row found in: ", log_path)
    tab <- read.table(text = paste(c(block[[h]], rows[keep]), collapse = "\n"),
                      header = TRUE, sep = "", check.names = FALSE, stringsAsFactors = FALSE,
                      na.strings = c("NA", "nan"))
    if (nrow(tab) != 1L) stop("Expected one LDSC pair in ", log_path, "; got ", nrow(tab), ".")
    tab
  }

  anchor_i <- which(traits$phenotype == "Loneliness isolation PanUKB")
  others <- setdiff(seq_len(nrow(traits)), anchor_i)
  array_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (nzchar(array_id)) {
    slot <- suppressWarnings(as.integer(array_id))
    if (is.na(slot) || slot < 1L) stop("Invalid SLURM_ARRAY_TASK_ID: ", array_id)
    if (slot > length(others)) {
      message("Array task ", slot, " has no phenotype pair; nothing to do.")
      quit(save = "no", status = 0L, runLast = FALSE)
    }
    slots <- slot
  } else {
    slots <- seq_along(others)
  }

  for (slot in slots) {
    target_i <- others[[slot]]
    out_prefix <- file.path(results_dir, sprintf("loneliness_%03d_%s", target_i,
                                                  safe_label(traits$phenotype[[target_i]])))
    signature <- c(
      normalizePath(traits$path[[anchor_i]], mustWork = TRUE),
      as.character(file.info(traits$path[[anchor_i]])$size),
      as.character(file.info(traits$path[[anchor_i]])$mtime),
      normalizePath(traits$path[[target_i]], mustWork = TRUE),
      as.character(file.info(traits$path[[target_i]])$size),
      as.character(file.info(traits$path[[target_i]])$mtime),
      normalizePath(ldsc_py, mustWork = TRUE)
    )
    signature_path <- paste0(out_prefix, ".signature")
    log_path <- paste0(out_prefix, ".log")
    signature_matches <- file.exists(signature_path) &&
      identical(readLines(signature_path, warn = FALSE), signature)
    tab <- if (signature_matches && file.exists(log_path)) {
      tryCatch(extract_rg_result(log_path), error = function(e) NULL)
    } else NULL

    if (is.null(tab)) {
      cmd <- c(ldsc_py, "--rg", paste(normalizePath(traits$path[c(anchor_i, target_i)], mustWork = TRUE), collapse = ","),
               "--ref-ld-chr", ref_ld_chr, "--w-ld-chr", w_ld_chr, "--out", out_prefix)
      console <- system2(python, cmd, stdout = TRUE, stderr = TRUE)
      status <- attr(console, "status")
      if (is.null(status)) status <- 0L
      writeLines(console, paste0(out_prefix, ".console.log"))
      if (status != 0L) {
        stop("LDSC failed for ", traits$phenotype[[target_i]], "; see ", out_prefix, ".console.log")
      }
      tab <- extract_rg_result(log_path)
      writeLines(signature, signature_path)
      message("Completed loneliness correlation with ", traits$phenotype[[target_i]])
    } else {
      message("Reusing completed loneliness correlation with ", traits$phenotype[[target_i]])
    }

    result <- data.table(
      exposure = "Loneliness isolation PanUKB", exposure_source = traits$source[[anchor_i]],
      trait = traits$phenotype[[target_i]], trait_source = traits$source[[target_i]],
      trait_origin = traits$origin[[target_i]], rg = as.numeric(tab$rg), se = as.numeric(tab$se),
      z = as.numeric(tab$z), p = as.numeric(tab$p))
    fwrite(result, paste0(out_prefix, ".result.tsv"), sep = "\t", na = "")
  }
  message("LDSC task complete. Pairwise result files are saved in ", results_dir,
          ". Run the aggregation job after all array tasks finish.")
  quit(save = "no", status = 0L, runLast = FALSE)
}

combined <- fread(combined_manifest_path, sep = "\t", strip.white = FALSE, showProgress = FALSE)
clean <- fread(clean_manifest_path, sep = "\t", strip.white = FALSE, showProgress = FALSE)
primary <- fread(primary_manifest_path, sep = "\t", strip.white = FALSE, showProgress = FALSE)
if (!all(c("phenotype", "path") %in% names(combined))) stop("Combined manifest must contain phenotype and path.")
if (!all(c("phenotype", "source", "clean_path") %in% names(clean))) stop("Clean manifest must contain phenotype, source, and clean_path.")
if (!all(c("phenotype", "source", "input_path") %in% names(primary))) stop("Primary manifest must contain phenotype, source, and input_path.")
if (anyDuplicated(primary$phenotype)) stop("Primary manifest has duplicate exact phenotype names.")

# Use only exact-name pre-munged files that belong to this Primary analysis.
pre <- fread(already_munged_path, strip.white = FALSE, na.strings = c("", "NA", "N/A"), showProgress = FALSE)
if (!all(c("phenotype", "source", "path") %in% names(pre))) stop("already_munged.csv needs phenotype, source, and path columns.")
pre[, `:=`(phenotype = as.character(phenotype), source = trimws(as.character(source)), path = trimws(as.character(path)))]
pre <- pre[phenotype %chin% primary$phenotype]
r13_primary <- primary[grepl("^https?://", input_path, ignore.case = TRUE) &
                         grepl("finngen-public-data-r13", input_path, ignore.case = TRUE) &
                         grepl("finngen", paste(source, cohort), ignore.case = TRUE), phenotype]
pre <- pre[!phenotype %chin% r13_primary]
if (anyDuplicated(pre$phenotype)) stop("already_munged.csv has multiple exact-name paths for a Primary phenotype.")
clean <- clean[phenotype %chin% primary$phenotype & !phenotype %chin% pre$phenotype]

# Existing LDSC-munged inputs are used directly; full summary stats are munged only when needed.
inputs <- list()
sign_checks <- list()
complement_allele <- function(x) {
  chartr("ACGT", "TGCA", toupper(as.character(x)))
}

# Confirm LDSC preserved the input effect direction after allele alignment.
check_sign_consistency <- function(input_path, munged_path, phenotype) {
  raw <- fread(input_path, select = c("SNP", "A1", "A2", "BETA"), showProgress = FALSE)
  munged <- fread(munged_path, select = c("SNP", "A1", "A2", "Z"), showProgress = FALSE)
  joined <- merge(raw, munged, by = "SNP", suffixes = c("_input", "_ldsc"))
  joined <- joined[is.finite(BETA) & is.finite(Z) & BETA != 0 & Z != 0]
  if (!nrow(joined)) stop("No variants available for LDSC sign check: ", phenotype)

  ia1 <- toupper(trimws(joined$A1_input))
  ia2 <- toupper(trimws(joined$A2_input))
  oa1 <- toupper(trimws(joined$A1_ldsc))
  oa2 <- toupper(trimws(joined$A2_ldsc))
  alleles_present <- !is.na(ia1) & !is.na(ia2) & !is.na(oa1) & !is.na(oa2)
  direct <- alleles_present & ia1 == oa1 & ia2 == oa2
  swapped <- alleles_present & ia1 == oa2 & ia2 == oa1
  ambiguous <- alleles_present & paste0(pmin(ia1, ia2), pmax(ia1, ia2)) %chin% c("AT", "CG")
  comp_direct <- alleles_present & !ambiguous & complement_allele(ia1) == oa1 & complement_allele(ia2) == oa2
  comp_swapped <- alleles_present & !ambiguous & complement_allele(ia1) == oa2 & complement_allele(ia2) == oa1
  orientation <- rep(NA_real_, nrow(joined))
  orientation[direct] <- 1
  orientation[!direct & swapped] <- -1
  orientation[!direct & !swapped & comp_direct] <- 1
  orientation[!direct & !swapped & !comp_direct & comp_swapped] <- -1

  comparable <- !is.na(orientation)
  if (sum(comparable) < 1000L) {
    stop("Fewer than 1,000 allele-aligned variants for LDSC sign check: ", phenotype,
         " (", sum(comparable), " comparable variants)")
  }
  expected_sign <- sign(joined$BETA[comparable]) * orientation[comparable]
  concordant <- sign(joined$Z[comparable]) == expected_sign
  rate <- mean(concordant)
  result <- data.table(
    phenotype = phenotype,
    variants_compared = sum(comparable),
    variants_concordant = sum(concordant),
    concordance = rate,
    discordance = 1 - rate,
    status = if (rate >= 0.99) "PASS" else "FAIL"
  )
  if (rate < 0.99) {
    stop("LDSC sign consistency failed for ", phenotype, ": ",
         sprintf("%.3f%% concordant", 100 * rate),
         ". See the input and munged alleles; this may indicate a sign reversal.")
  }
  result
}

munge_one <- function(path, phenotype, source, origin, index) {
  if (!file.exists(path)) stop("Input file missing for ", phenotype, ": ", path)
  prefix <- file.path(munged_dir, sprintf("%s_%03d_%s__%s", origin, index, safe_label(phenotype), safe_label(source)))
  expected <- paste0(prefix, ".sumstats.gz")
  munged_is_complete <- file.exists(expected) && !is.na(file.info(expected)$size) &&
    file.info(expected)$size > 0 &&
    identical(as.integer(system2("gzip", c("-t", expected), stdout = FALSE, stderr = FALSE)), 0L)
  if (!munged_is_complete) {
    if (file.exists(expected)) unlink(expected)
    headers <- names(fread(path, nrows = 0L, showProgress = FALSE))
    sample_size_args <- if ("N" %in% headers) {
      c("--N-col", "N")
    } else if (all(c("N_CASE", "N_CONTROL") %in% headers)) {
      c("--N-cas-col", "N_CASE", "--N-con-col", "N_CONTROL")
    } else {
      stop("No usable sample-size columns for ", phenotype,
           ". Expected N or both N_CASE and N_CONTROL in ", path)
    }
    output <- system2(python, c(munge_py, "--sumstats", path, "--out", prefix,
                                "--merge-alleles", hm3_snplist,
                                "--snp", "SNP", "--a1", "A1", "--a2", "A2",
                                "--p", "P", "--signed-sumstats", "BETA,0",
                                sample_size_args),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    if (status != 0L) stop("munge_sumstats.py failed for ", phenotype, ":\n", paste(output, collapse = "\n"))
  }
  if (!file.exists(expected) || file.info(expected)$size == 0) stop("Munging did not create a complete output: ", expected)
  sign_checks[[length(sign_checks) + 1L]] <<- check_sign_consistency(path, expected, phenotype)
  data.table(phenotype = phenotype, source = source, origin = origin,
             path = normalizePath(expected, mustWork = TRUE))
}

for (i in seq_len(nrow(combined))) {
  inputs[[length(inputs) + 1L]] <- munge_one(as.character(combined$path[[i]]),
    as.character(combined$phenotype[[i]]), "METAL_IVW", "METAL_combined", i)
}
for (i in seq_len(nrow(clean))) {
  inputs[[length(inputs) + 1L]] <- munge_one(as.character(clean$clean_path[[i]]),
    as.character(clean$phenotype[[i]]), as.character(clean$source[[i]]), "Primary_clean", i)
}
for (i in seq_len(nrow(pre))) {
  path <- as.character(pre$path[[i]])
  if (!file.exists(path) || file.access(path, 4L) != 0L) stop("Pre-munged file missing for ", pre$phenotype[[i]], ": ", path)
  headers <- names(fread(path, nrows = 0L, showProgress = FALSE))
  if (!all(c("SNP", "N", "Z") %in% headers)) stop("Not LDSC-munged for ", pre$phenotype[[i]], ": expected SNP, N, Z columns.")
  sign_checks[[length(sign_checks) + 1L]] <- data.table(
    phenotype = as.character(pre$phenotype[[i]]), variants_compared = NA_integer_,
    variants_concordant = NA_integer_, concordance = NA_real_, discordance = NA_real_,
    status = "NOT CHECKED: pre-munged file has no raw BETA/allele comparison in this step")
  inputs[[length(inputs) + 1L]] <- data.table(phenotype = as.character(pre$phenotype[[i]]),
    source = as.character(pre$source[[i]]), origin = "already_munged", path = normalizePath(path, mustWork = TRUE))
}
traits <- rbindlist(inputs, fill = TRUE)
if (!nrow(traits)) stop("No LDSC inputs were assembled.")
if (anyDuplicated(traits$phenotype)) stop("Duplicate trait names across combined METAL and Primary inputs: ",
  paste(traits$phenotype[duplicated(traits$phenotype)], collapse = ", "))
if (!any(traits$phenotype == "Loneliness isolation PanUKB")) {
  stop("Loneliness isolation PanUKB is not available as an LDSC-ready Primary input.")
}

# Persist the prepared list and QC report on scratch. Stage two reads this
# manifest directly and can be rerun without repeating successful munging.
traits[, `:=`(ldsc_name = sprintf("trait_%03d_%s__%s", .I,
                                  safe_label(phenotype), safe_label(source)),
              ldsc_path = path)]
input_manifest_path <- file.path(munged_dir, "ldsc_inputs_manifest.tsv")
fwrite(traits, input_manifest_path, sep = "\t", na = "")
fwrite(rbindlist(sign_checks, fill = TRUE),
       file.path(munged_dir, "ldsc_sign_consistency.tsv"), sep = "\t", na = "")
message("LDSC preparation complete: ", nrow(traits), " inputs listed in ", input_manifest_path)
message("Submit submit_metal_primary_ldsc.sh next to run resumable correlations.")
quit(save = "no", status = 0L, runLast = FALSE)
