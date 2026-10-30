#!/usr/bin/env Rscript
# run_ldsr_rg.R -- pairwise LDSC genetic correlations for the poster trait set
#
# Runs ldsc.py --rg for every pair of traits in the traits file, then parses
# the logs into a long table and a symmetric rg matrix.
# Pairs containing >=1 case-control trait get --samp-prev 0.5 --pop-prev 0.5
# (liability conversion is identity under the effective-N convention, and
# ldsc requires both flags together).
#
# Usage (interactive, on the cluster):
#   Rscript run_ldsr_rg.R
#   LDSR_FILTER="Loneliness|CRP" Rscript run_ldsr_rg.R   # quick smoke test
#   LDSR_DRY_RUN=1 Rscript run_ldsr_rg.R
#
# Env vars: same as run_ldsr_h2.R (LDSR_TRAITS optional; defaults to the
#                  manifest phenotypes), plus
#   LDSR_MAX_PAIRS  cap on number of pairs, after lexicographic sort
#
# Outputs in OUT_DIR:
#   <trait1>__<trait2>.log          one per pair
#   rg_long.tsv                     trait1, trait2, rg, se, z, p
#   rg_matrix.tsv                   n x n symmetric matrix (diagonal = 1)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

here        <- function(x) file.path(getwd(), x)
traits_file <- Sys.getenv("LDSR_TRAITS",    unset = "")  # empty = use POSTER_MANIFEST phenotypes
manifest_f  <- Sys.getenv("POSTER_MANIFEST",unset = "munge_manifest_poster.tsv")
munged_dir  <- sub("/$", "", Sys.getenv("MUNGED_DIR",  unset = here("munged")))
ldsc_dir    <- sub("/$", "", Sys.getenv("LDSC_DIR",    unset = here("ldsc_py3")))
ref_ld      <- sub("/$", "", Sys.getenv("REF_LD",      unset = "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr"))
w_ld        <- sub("/$", "", Sys.getenv("W_LD",        unset = ref_ld))
out_dir     <- sub("/$", "", Sys.getenv("OUT_DIR",     unset = here("ldsr_rg")))
py          <- Sys.getenv("MUNGE_PYTHON", unset = "python3")
n_cores     <- as.integer(Sys.getenv("LDSR_CORES", unset = "6"))
dry_run     <- Sys.getenv("LDSR_DRY_RUN",  unset = "0") == "1"
filter_str  <- Sys.getenv("LDSR_FILTER",   unset = "")
max_pairs   <- Sys.getenv("LDSR_MAX_PAIRS", unset = "")

if (nchar(traits_file) > 0 && !file.exists(traits_file)) {
  stop("LDSR_TRAITS file not found: ", traits_file)
}
if (!file.exists(manifest_f)) {
  stop("Required file not found: ", manifest_f, "\n(Set POSTER_MANIFEST)")
}
man <- fread(manifest_f)
stopifnot(all(c("phenotype", "n_policy") %in% names(man)))
traits <- if (nchar(traits_file) > 0) fread(traits_file, header = FALSE)$V1 else man$phenotype

if (nchar(filter_str) > 0) {
  keep <- strsplit(filter_str, "[ ,]+")[[1]]
  traits <- traits[grepl(paste(keep, collapse = "|"), traits)]
}
traits <- traits[file.exists(file.path(munged_dir, paste0(gsub("\\s+", "_", traits), ".sumstats.gz")))]
if (length(traits) < 2) stop("Need >= 2 traits with munged sumstats in ", munged_dir)

## ---- environment guard: R startup can set PYTHONPATH/PYTHONHOME and
##      break the spawned python (py3.8 bundle shadows the venv); sanitize ----
Sys.unsetenv("PYTHONPATH")
Sys.unsetenv("PYTHONHOME")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
s_names  <- sanitize(traits)
is_cc    <- man[phenotype %in% traits, .(phenotype, is_cc = n_policy[1] == "ncase_ncontrol")]

pairs <- utils::combn(length(traits), 2)
ord <- order(pairs[1, ], pairs[2, ])
pairs <- pairs[, ord, drop = FALSE]
if (nchar(max_pairs) > 0) {
  pairs <- pairs[, seq_len(min(as.integer(max_pairs), ncol(pairs))), drop = FALSE]
}

cmds <- apply(pairs, 2, function(idx) {
  i <- idx[1]; j <- idx[2]
  cc <- is_cc[phenotype %in% traits[c(i, j)], any(is_cc)]
  prev <- if (isTRUE(cc)) " --samp-prev 0.5,0.5 --pop-prev 0.5,0.5" else ""
  sprintf("%s %s --rg %s,%s --ref-ld-chr %s/ --w-ld-chr %s/ --out %s/%s__%s%s",
          py, file.path(ldsc_dir, "ldsc.py"),
          file.path(munged_dir, paste0(s_names[i], ".sumstats.gz")),
          file.path(munged_dir, paste0(s_names[j], ".sumstats.gz")),
          ref_ld, w_ld, out_dir, s_names[i], s_names[j], prev)
})

cat("== run_ldsr_rg:", length(cmds), "pairs,", n_cores, "cores,",
    if (dry_run) "DRY RUN" else "live", "==\n\n")
if (length(cmds) <= 10) cat(paste0(cmds, "\n"), sep = "") else
  cat(paste0(c(head(cmds, 5), "...", tail(cmds, 5)), "\n"), sep = "")

if (dry_run) quit(save = "no", status = 0)

# parse "Genetic Correlation: RG (SE)" from a log; NA if absent
parse_rg <- function(txt) {
  if (!grepl("Genetic Correlation: -?[0-9]", txt)) return(c(NA_real_, NA_real_))
  m <- sub(".*Genetic Correlation: (-?[0-9.eE+-]+) \\((-?[0-9.eE+-]+)\\).*", "\\1 \\2", txt)
  as.numeric(strsplit(m, " ")[[1]])
}

res <- mclapply(seq_along(cmds), function(i) {
  i1 <- pairs[1, i]; i2 <- pairs[2, i]
  logf <- file.path(out_dir, paste0(s_names[i1], "__", s_names[i2], ".log"))
  status <- system(paste0(cmds[i], " > ", shQuote(paste0(logf, ".stdout")), " 2>&1"))
  txt <- paste(readLines(logf, warn = FALSE), collapse = "\n")
  rgse <- parse_rg(txt)
  zp <- function(pat) {
    m <- regmatches(txt, regexpr(paste0(pat, "-?[0-9.eE+-]+"), txt))
    if (length(m) == 0) return(NA_real_)
    as.numeric(sub(pat, "", m))
  }
  data.table(trait1 = traits[i1], trait2 = traits[i2], exit_code = status,
             rg = rgse[1], se = rgse[2], z = zp("Z-score: "), p = zp("P: "))
}, mc.cores = n_cores)

long <- rbindlist(res, fill = TRUE)
setkey(long, trait1, trait2)
cat("\n== rg run summary ==\n")
cat("Pairs run:", nrow(long), "| failed:", sum(long$exit_code != 0),
    "| unparsed rg:", sum(is.na(long$rg)), "\n")
if (any(long$exit_code != 0) || any(is.na(long$rg))) {
  cat("Failed/unparsed pairs:\n")
  print(long[exit_code != 0 | is.na(rg), .(trait1, trait2, exit_code)])
}

fwrite(long, file.path(out_dir, "rg_long.tsv"), sep = "\t")

n <- length(traits)
mat <- matrix(NA_real_, n, n, dimnames = list(traits, traits))
diag(mat) <- 1
ok <- long[!is.na(rg) & exit_code == 0]
if (nrow(ok) > 0) {
  mat[cbind(ok$trait1, ok$trait2)] <- ok$rg
  mat[cbind(ok$trait2, ok$trait1)] <- ok$rg
}
fwrite(as.data.table(mat, keep.rownames = "trait"),
       file.path(out_dir, "rg_matrix.tsv"), sep = "\t")
cat("\nWritten:", file.path(out_dir, "rg_long.tsv"), "and",
    file.path(out_dir, "rg_matrix.tsv"), "\n")