#!/usr/bin/env Rscript
# ============================================================================
# run_munge.R — munge all traits in munge_manifest.tsv via LDSC's
#               munge_sumstats.py, with per-trait N policy from the manifest.
#
# N-policy logic (matches n_overrides.tsv / GenomicSEM recommendation):
#   ncase_ncontrol -> --N <Neff>  (GenomicSEM: effective N = 4*Nc*Nt/(Nc+Nt);
#                     set sample prevalence = 0.5 in downstream LDSC h2 runs)
#   per_snp_n      -> no N flag; the file's per-SNP N column is used
#
# Column-collision handling (munge_sumstats.py precedence traps):
#   - N_CASE/N_CONTROL columns are ALWAYS ignored (munge would otherwise use
#     Ncas+Ncon, i.e. total N — wrong for our Neff policy, and all-NA columns
#     would zero out the file)
#   - the N column is ignored for ncase_ncontrol traits (munge uses a file N
#     column INSTEAD of --N otherwise)
#   - BETA is ignored when Z is also present (munge errors on two signed-
#     sumstat columns; Z is the preferred signed statistic)
#   - FRQ is ignored when entirely NA (e.g. Loneliness; NA FRQ would drop
#     every SNP at the MAF filter)
#
# Usage:
#   Interactive: edit CONFIG paths below, then source("run_munge.R")
#   Headless:    Rscript run_munge.R [manifest_path] [munged_dir]
#
# Env overrides:
#   MUNGE_MANIFEST, MUNGED_DIR, LDSC_DIR (dir with munge_sumstats.py),
#   HM3 (w_hm3.snplist), MUNGE_PYTHON (default: python),
#   MUNGE_CORES (default 4), MUNGE_DRY_RUN (1 = print commands only),
#   MUNGE_FILTER (comma-separated phenotypes to run a subset)
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

ROOT <- "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline"
manifest_path <- Sys.getenv("MUNGE_MANIFEST", file.path(ROOT, "munge_manifest/munge_manifest.tsv"))
munged_dir    <- Sys.getenv("MUNGED_DIR",    file.path(ROOT, "munged"))
ldsc_dir      <- Sys.getenv("LDSC_DIR",      "/gpfs/gibbs/project/polimanti/ejm96/immune/ldsc_py3")
hm3           <- Sys.getenv("HM3",           "/gpfs/gibbs/pi/polimanti/LDSC_files/LDSCORE_w_hm3.snplist")
py            <- Sys.getenv("MUNGE_PYTHON",  "python3")
n_cores       <- as.integer(Sys.getenv("MUNGE_CORES", "4"))
dry_run       <- Sys.getenv("MUNGE_DRY_RUN", "0") == "1"
filter_phenos <- Sys.getenv("MUNGE_FILTER", "")  # e.g. "Schizophrenia PGC,Loneliness isolation PanUKB"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) manifest_path <- args[1]
if (length(args) >= 2) munged_dir    <- args[2]

munge_py <- file.path(ldsc_dir, "munge_sumstats.py")

## ---- validation ----
stopifnot("Manifest not found" = file.exists(manifest_path))
stopifnot("munge_sumstats.py not found (set LDSC_DIR)" = file.exists(munge_py))
stopifnot("HapMap3 snplist not found (set HM3)" = file.exists(hm3))

m <- fread(manifest_path)
need <- c("phenotype", "clean_path", "n_policy", "Neff")
stopifnot("Manifest missing columns" = all(need %in% names(m)))
if (nzchar(filter_phenos)) m <- m[phenotype %in% strsplit(filter_phenos, ",")[[1]]]

# completeness: every ncase_ncontrol row needs Neff; every per_snp_n row
# needs a per-SNP N median (else munge cannot determine N)
bad1 <- m[n_policy == "ncase_ncontrol" & is.na(Neff), phenotype]
bad2 <- m[n_policy == "per_snp_n" & is.na(n_median), phenotype]
if (length(bad1)) stop("ncase_ncontrol traits with NA Neff: ", paste(bad1, collapse = ", "))
if (length(bad2)) stop("per_snp_n traits with NA n_median: ", paste(bad2, collapse = ", "))

missing_files <- m[!file.exists(clean_path), clean_path]
if (length(missing_files)) stop("Clean files not found:\n", paste(missing_files, collapse = "\n"))

dir.create(munged_dir, recursive = TRUE, showWarnings = FALSE)

## ---- per-trait command builder ----
build_cmd <- function(phenotype, clean_path, n_policy, Neff) {
  # header only (fast on .gz)
  hdr <- names(fread(clean_path, nrows = 0, showProgress = FALSE))
  ignore <- character(0)
  
  # case/control N columns: always drop (constant sheet-filled or real per-SNP —
  # neither may drive N under our policy). Covers both header spellings:
  # N_CASE/N_CONTROL (standard cleaner output) and N_CAS/N_CON (Crohn,
  # Hypothyroidism — munge maps these to the same fields).
  cc_cols <- intersect(c("N_CASE", "N_CONTROL", "N_CAS", "N_CON"), hdr)
  ignore <- c(ignore, cc_cols)
  
  # N column: drop when policy says constant Neff (munge prefers file N over --N)
  if ("N" %in% hdr && n_policy == "ncase_ncontrol") ignore <- c(ignore, "N")
  
  # signed sumstat: pass explicitly (version-proof — some py3 munge ports
  # fail auto-detection after --ignore BETA). Z preferred; BETA for files
  # without Z (METAL traits).
  signed <- if ("Z" %in% hdr) "Z,0" else "BETA,0"
  if ("Z" %in% hdr && "BETA" %in% hdr) ignore <- c(ignore, "BETA")
  
  # FRQ / INFO: drop if present but entirely NA (NA FRQ -> every SNP fails the
  # MAF filter; NA INFO -> every SNP fails the INFO > 0.9 filter)
  for (col in c("FRQ", "INFO")) {
    if (col %in% hdr) {
      head_vals <- fread(clean_path, nrows = 5000, select = col, showProgress = FALSE)[[col]]
      if (all(is.na(head_vals))) ignore <- c(ignore, col)
    }
  }
  
  out_prefix <- file.path(munged_dir, gsub("[^A-Za-z0-9._-]", "_", phenotype))
  cmd <- c(munge_py, "--sumstats", clean_path,
           "--out", out_prefix,
           "--merge-alleles", hm3,
           "--signed-sumstats", signed)
  if (n_policy == "ncase_ncontrol") cmd <- c(cmd, "--N", format(Neff, scientific = FALSE))
  if (length(ignore))               cmd <- c(cmd, "--ignore", paste(ignore, collapse = ","))
  list(cmd = cmd, out_prefix = out_prefix, ignore = ignore)
}

cat(sprintf("Munging %d traits -> %s (cores: %d%s)\n",
            nrow(m), munged_dir, n_cores, if (dry_run) ", DRY RUN" else ""))

jobs <- lapply(seq_len(nrow(m)), function(i) {
  b <- do.call(build_cmd, as.list(m[i, .(phenotype, clean_path, n_policy, Neff)]))
  list(phenotype = m$phenotype[i], cmd = b$cmd, out_prefix = b$out_prefix, ignore = b$ignore)
})

for (j in jobs) {
  cat(sprintf("  %-40s %s\n", j$phenotype, paste(j$cmd, collapse = " ")))
  if (length(j$ignore)) cat(sprintf("    ignore: %s\n", paste(j$ignore, collapse = ",")))
}
if (dry_run) cat("Dry run — nothing executed.\n")

## ---- execution ----
run_one <- function(j) {
  err_file <- paste0(j$out_prefix, ".stderr")
  out_gz <- paste0(j$out_prefix, ".sumstats.gz")
  status <- suppressWarnings(system2(py, j$cmd,
                                     stdout = NULL, stderr = err_file))
  list(phenotype = j$phenotype, exit = status, err_file = err_file, out_gz = out_gz)
}

res <- if (dry_run) {
  NULL
} else if (n_cores > 1) {
  mclapply(jobs, run_one, mc.cores = n_cores)
} else {
  lapply(jobs, run_one)
}

## ---- summary ----
if (!dry_run) {
  summ <- rbindlist(lapply(res, function(r) {
    err <- if (r$exit != 0 && file.exists(r$err_file)) {
      paste(tail(readLines(r$err_file, warn = FALSE), 3), collapse = " | ")
    } else ""
    data.table(phenotype = r$phenotype,
               exit_code = r$exit,
               sumstats_gz = file.exists(r$out_gz),
               size_kb = round(file.size(r$out_gz) / 1024, 1),
               note = if (nchar(err)) paste(err, collapse = " | ") else "")
  }))
  cat("\n==== munge summary ====\n")
  print(summ, width = 200)
  ok <- summ[exit_code == 0 & sumstats_gz == TRUE, .N]
  cat(sprintf("\n%d/%d traits munged successfully.\n", ok, nrow(summ)))
  cat("Per-trait logs: <munged_dir>/<phenotype>.log ; stderr: <phenotype>.stderr\n")
  cat("REMINDER: for ncase_ncontrol traits, use --samp-prev 0.5 (not the true\n")
  cat("prevalence) in downstream LDSC h2/rg runs — Neff already corrects for\n")
  cat("ascertainment (GenomicSEM recommendation, Grotzinger et al. 2023).\n")
}