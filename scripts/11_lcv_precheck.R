
# 11_LCV_precheck.R
# Check the munged GWAS inputs, align alleles pairwise, attach European LD
# scores, sort variants by genomic position, and save the LCV job input RDS.

suppressPackageStartupMessages(library(data.table))

# Munged GWAS files used in the LCV comparisons.
gwas_files <- c(
  Loneliness = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/munged/Loneliness_isolation_PanUKB.sumstats.gz",
  Triglycerides = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/munged/Triglycerides.sumstats.gz",
  Depression = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/munged/Depression_PGC.sumstats.gz",
  T2D = "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/munged/Type_2_Diabetes.sumstats.gz"
)

ld_score_dir <- "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr"
output_rds <- "/vast/palmer/scratch/polimanti/ejm96/lcv/lcv_inputs_aligned_ordered.rds"


# Read the LCV fields, report data quality, and return usable SNP rows.
read_lcv_gwas <- function(trait, path) {
  cat("\n========== ", trait, " ==========\n", sep = "")
  cat("File: ", path, "\n", sep = "")
  
  if (!file.exists(path)) stop("GWAS file not found: ", path)
  
  header <- names(fread(path, nrows = 0, showProgress = FALSE))
  required <- c("SNP", "A1", "A2", "Z", "N")
  missing <- setdiff(required, header)
  if (length(missing)) {
    stop(trait, " is missing required columns: ", paste(missing, collapse = ", "))
  }
  
  cat("Columns: ", paste(header, collapse = ", "), "\n", sep = "")
  x <- fread(path, select = required, showProgress = TRUE)
  x[, `:=`(
    SNP = as.character(SNP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    Z = as.numeric(Z),
    N = as.numeric(N)
  )]
  
  valid <- !is.na(x$SNP) & nzchar(x$SNP) &
    x$A1 %chin% c("A", "C", "G", "T") &
    x$A2 %chin% c("A", "C", "G", "T") &
   is.finite(x$Z) & is.finite(x$N) & x$N > 0

  cat("Rows: ", format(nrow(x), big.mark = ","),
     "\nUsable rows: ", format(sum(valid), big.mark = ","),
     "\nExcluded rows: ", format(sum(!valid), big.mark = ","), "\n", sep = "")

  x <- x[valid]

  # Remove all records for repeated SNP IDs rather than selecting an arbitrary one.
  duplicated_ids <- x[, .N, by = SNP][N > 1L, SNP]
 if (length(duplicated_ids)) {
   cat("Removing all records for ", length(duplicated_ids),
        " duplicated SNP IDs.\n", sep = "")
    x <- x[!SNP %chin% duplicated_ids]
  }

  if (!nrow(x)) stop("No usable unique SNPs remain for ", trait)

  cat("Usable unique SNPs: ", format(nrow(x), big.mark = ","),
     "\nN column: min=", signif(min(x$N), 7),
     "; median=", signif(median(x$N), 7),
     "; max=", signif(max(x$N), 7),
     "; distinct values=", uniqueN(x$N), "\n", sep = "")
  cat("Example usable rows:\n")
  print(x[seq_len(min(5L, nrow(x)))])

  x
}




valid_lcv_rows <- lapply(names(gwas_files), function(trait) {
  read_lcv_gwas(trait, gwas_files[[trait]])
})
names(valid_lcv_rows) <- names(gwas_files)

# Match a pair by rsID and count same, swapped, or incompatible allele pairs.
check_pair_alignment <- function(trait1, trait2) {
  x <- copy(valid_lcv_rows[[trait1]])
  y <- copy(valid_lcv_rows[[trait2]])
  
  pair <- merge(x, y, by = "SNP", suffixes = c("_1", "_2"))
  pair[, allele_match := fcase(
    A1_1 == A1_2 & A2_1 == A2_2, "same",
    A1_1 == A2_2 & A2_1 == A1_2, "swapped",
    default = "different or strand issue"
  )]
  
  cat("\n", trait1, " vs ", trait2,
      "\nOverlapping SNPs: ", format(nrow(pair), big.mark = ","), "\n", sep = "")
  print(pair[, .N, by = allele_match][order(-N)])
  pair
}

# Align trait 2's Z scores to trait 1's effect-allele orientation.
align_pair_for_lcv <- function(pair, label) {
  pair <- copy(pair)
  pair <- pair[allele_match %chin% c("same", "swapped")]
  n_swapped <- pair[allele_match == "swapped", .N]
  pair[allele_match == "swapped", Z_2 := -Z_2]
  
  cat(label, ": retained ", format(nrow(pair), big.mark = ","),
      " allele-matched SNPs; flipped Z for ",
      format(n_swapped, big.mark = ","), " swapped-allele SNPs.\n", sep = "")
  
  pair[, .(SNP, Z_1, Z_2, N_1, N_2)]
}

# Load all 22 chromosomes of the EUR per-SNP LD-score reference.
ld_files <- file.path(ld_score_dir, paste0(1:22, ".l2.ldscore.gz"))
missing_ld <- ld_files[!file.exists(ld_files)]
if (length(missing_ld)) {
  stop("Missing LD-score files: ", paste(basename(missing_ld), collapse = ", "))
}

ld_scores <- rbindlist(lapply(ld_files, function(path) {
  fread(path, select = c("CHR", "SNP", "BP", "MAF", "L2"),
        showProgress = FALSE)
}))
ld_scores[, `:=`(CHR = as.integer(CHR), BP = as.numeric(BP))]

if (anyDuplicated(ld_scores$SNP)) {
  stop("LD-score reference contains duplicated SNP IDs.")
}

# Attach reference coordinates and L2, apply MAF/MHC filters, and sort.
attach_l2_and_order <- function(x, label) {
  x <- merge(
    as.data.table(x),
    ld_scores[, .(SNP, CHR, BP, MAF, L2)],
    by = "SNP",
    all = FALSE,
    sort = FALSE
  )
  
  n_matched <- nrow(x)
  x <- x[is.finite(MAF) & MAF > 0.05 &
           !(CHR == 6L & BP >= 25e6 & BP <= 34e6)]
  setorder(x, CHR, BP)
  
  if (!nrow(x)) stop("No variants remain after LD-score matching for ", label)
  if (anyDuplicated(x$SNP)) stop("Duplicate SNPs remain for ", label)
  if (any(!is.finite(x$L2)) || any(x$L2 < 0)) {
    stop("Invalid L2 values remain for ", label)
  }
  
  cat(label, ": matched ", format(n_matched, big.mark = ","),
      "; retained after MAF/MHC filters ", format(nrow(x), big.mark = ","),
      " variants.\n", sep = "")
  x
}

# Specify one input direction per unique trait pair.
pair_specs <- list(
  T2D_to_TG = c("T2D", "Triglycerides"),
  TG_to_Loneliness = c("Triglycerides", "Loneliness"),
  Depression_to_T2D = c("Depression", "T2D"),
  T2D_to_Loneliness = c("T2D", "Loneliness"),
  Loneliness_to_Depression = c("Loneliness", "Depression"),
  Depression_to_TG = c("Depression", "Triglycerides")
)

lcv_inputs <- lapply(names(pair_specs), function(pair_name) {
  traits <- pair_specs[[pair_name]]
  pair <- check_pair_alignment(traits[1], traits[2])
  aligned <- align_pair_for_lcv(pair, pair_name)
  attach_l2_and_order(aligned, pair_name)
})
names(lcv_inputs) <- names(pair_specs)

dir.create(dirname(output_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(lcv_inputs, output_rds)

cat("\nSaved LCV inputs to: ", output_rds, "\n", sep = "")
cat("Input columns by comparison:\n")
print(lapply(lcv_inputs, names))