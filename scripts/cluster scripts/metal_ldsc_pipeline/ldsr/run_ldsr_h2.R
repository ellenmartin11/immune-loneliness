#!/usr/bin/env Rscript
# run_ldsr_h2.R -- univariate LDSC h2 for the poster trait set
#
# Runs ldsc.py --h2 for every trait listed in the traits file, using the
# munged sumstats produced by run_munge.R. Case-control traits (n_policy
# == "ncase_ncontrol" in the manifest) get --samp-prev 0.5 --pop-prev 0.5,
# which makes liability h2 == observed h2 -- the convention that matches
# munging with effective N (GenomicSEM guidance, Grotzinger et al. 2023).
#
# Usage (interactive, on the cluster):
#   Rscript run_ldsr_h2.R
#   MUNGE_FILTER="Loneliness" Rscript run_ldsr_h2.R     # subset of traits
#   LDSR_DRY_RUN=1 Rscript run_ldsr_h2.R                # print commands only
#
# Env vars (defaults in parentheses):
#   LDSR_TRAITS    optional file with one phenotype per line;
#                  omit to use all phenotypes in POSTER_MANIFEST
#   POSTER_MANIFEST  manifest with phenotype + n_policy columns
#                  (munge_manifest_poster.tsv next to this script)
#   MUNGED_DIR     dir with <Trait>.sumstats.gz  (munged/ under cwd)
#   LDSC_DIR       dir with ldsc.py              (ldsc_py3 under cwd)
#   REF_LD         ref-ld-chr prefix  (/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr/)
#   W_LD           w-ld-chr prefix    (same as REF_LD)
#   OUT_DIR        where .log/.h2 files go       (ldsr_h2/ under cwd)
#   MUNGE_PYTHON   python interpreter            (python3 on PATH)
#   LDSR_CORES     parallel jobs                 (6)
#   LDSR_DRY_RUN   1 = print commands, run nothing
#   LDSR_FILTER    space-separated trait substrings to keep

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

here        <- function(x) file.path(getwd(), x)
traits_file <- Sys.getenv("LDSR_TRAITS",    unset = "")  # empty = use POSTER_MANIFEST phenotypes
manifest_f  <- Sys.getenv("POSTER_MANIFEST",unset = "munge_manifest_poster.tsv")
munged_dir  <- sub("/$", "", Sys.getenv("MUNGED_DIR",  unset = here("munged")))
ldsc_dir    <- sub("/$", "", Sys.getenv("LDSC_DIR",    unset = here("ldsc_py3")))
ref_ld      <- sub("/$", "", Sys.getenv("REF_LD",      unset = "/gpfs/gibbs/pi/polimanti/LDSC_files/eur_w_ld_chr"))
w_ld        <- sub("/$", "", Sys.getenv("W_LD",        unset = ref_ld))
out_dir     <- sub("/$", "", Sys.getenv("OUT_DIR",     unset = here("ldsr_h2")))
py          <- Sys.getenv("MUNGE_PYTHON", unset = "python3")
n_cores     <- as.integer(Sys.getenv("LDSR_CORES", unset = "6"))
dry_run     <- Sys.getenv("LDSR_DRY_RUN",  unset = "0") == "1"
filter_str  <- Sys.getenv("LDSR_FILTER",   unset = "")

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
if (length(traits) == 0) stop("No munged .sumstats.gz found for the requested traits in ", munged_dir)

## ---- environment guard: R startup can set PYTHONPATH/PYTHONHOME and
##      break the spawned python (py3.8 bundle shadows the venv); sanitize ----
Sys.unsetenv("PYTHONPATH")
Sys.unsetenv("PYTHONHOME")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

cmds <- vapply(traits, function(tr) {
  s     <- sanitize(tr)
  pol   <- man[phenotype == tr, n_policy]
  is_cc <- length(pol) > 0 && !is.na(pol[1]) && pol[1] == "ncase_ncontrol"
  prev  <- if (is_cc) " --samp-prev 0.5 --pop-prev 0.5" else ""
  sprintf("%s %s --h2 %s --ref-ld-chr %s/ --w-ld-chr %s/ --out %s/%s%s",
          py, file.path(ldsc_dir, "ldsc.py"),
          file.path(munged_dir, paste0(s, ".sumstats.gz")),
          ref_ld, w_ld, out_dir, s, prev)
}, character(1), USE.NAMES = FALSE)

cat("== run_ldsr_h2:", length(cmds), "traits,", n_cores, "cores,",
    if (dry_run) "DRY RUN" else "live", "==\n\n")
cat(paste0(cmds, "\n"), sep = "")

if (dry_run) quit(save = "no", status = 0)

res <- mclapply(seq_along(cmds), function(i) {
  tr <- traits[i]
  logf <- file.path(out_dir, paste0(sanitize(tr), ".log"))
  status <- system(paste0(cmds[i], " > ", shQuote(paste0(logf, ".stdout")), " 2>&1"))
  txt <- paste(readLines(logf, warn = FALSE), collapse = "\n")
  grab_num <- function(pat, grp = 1) {
    m <- regmatches(txt, regexpr(pat, txt))
    if (length(m) == 0) return(NA_real_)
    as.numeric(sub(pat, paste0("\\", grp), m))
  }
  # with --samp-prev/--pop-prev ldsc reports "Total Liability scale h2"
  # (== observed when both are 0.5); without, "Total Observed scale h2"
  data.table(trait = tr, exit_code = status,
             h2        = grab_num("Total (Observed|Liability) scale h2: (-?[0-9.eE+-]+)", grp = 2),
             h2_se     = grab_num("Total (Observed|Liability) scale h2: -?[0-9.eE+-]+ \\((-?[0-9.eE+-]+)\\)", grp = 2),
             intercept = grab_num("Intercept: (-?[0-9.eE+-]+)"),
             lambda_gc = grab_num("Lambda GC: (-?[0-9.eE+-]+)"),
             mean_chi2 = grab_num("Mean Chi\\^2: (-?[0-9.eE+-]+)"))
}, mc.cores = n_cores)

summary <- rbindlist(res, fill = TRUE)
setkey(summary, trait)
cat("\n== h2 run summary ==\n")
print(summary, row.names = FALSE)
cat("\nFailed:", sum(summary$exit_code != 0), "of", nrow(summary),
    "| unparsed h2:", sum(is.na(summary$h2)), "\n")
if (any(summary$exit_code != 0) || any(is.na(summary$h2))) {
  cat("Check logs in", out_dir, "for the failing traits.\n")
}
fwrite(summary, file.path(out_dir, "h2_summary.tsv"), sep = "\t")
cat("Written:", file.path(out_dir, "h2_summary.tsv"), "\n")
