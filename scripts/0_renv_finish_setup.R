############################################################################
### renv_finish_setup.R
###
### Run from INSIDE the immune project directory, under R/4.4.1:
###
###   cd /gpfs/gibbs/project/polimanti/ejm96/immune
###   module load R/4.4.1-foss-2022b
###   Rscript renv_finish_setup.R
###
### 1. Copies the prebuilt GenomicSEM stack from ~/Rlibs/4.4.1 into the
###    renv library (the path renv is actually using, .libPaths()[1]).
### 2. Copies svglite + its dependencies from the module site library.
### 3. Verifies the full closure is loadable.
###
### Then run:  Rscript -e 'renv::snapshot()'
############################################################################

lib <- .libPaths()[1]
cat("renv library:", lib, "\n")
if (!grepl("renv", lib)) {
  stop("renv does not appear to be active - run from the project dir under R/4.4.1")
}

## 1. Prebuilt GenomicSEM stack (compiled under this same R 4.4.1 module)
src_dirs <- list.dirs(path.expand("~/Rlibs/4.4.1"), recursive = FALSE)
if (length(src_dirs) == 0) stop("nothing found in ~/Rlibs/4.4.1")
for (d in src_dirs) {
  tgt <- file.path(lib, basename(d))
  if (dir.exists(tgt)) {
    cat("already   :", basename(d), "\n")
  } else {
    ok <- file.copy(d, lib, recursive = TRUE, copy.date = TRUE)
    cat(if (ok) "copied    :" else "FAILED    :", basename(d), "\n")
  }
}

## 2. svglite + dependencies from the module site library
for (p in c("svglite", "systemfonts", "textshaping", "rcmdcheck")) {
  src <- tryCatch(find.package(p), error = function(e) NULL)
  if (is.null(src)) {
    cat("NOT FOUND :", p, "(only needed if renv keeps asking)\n")
  } else if (!dir.exists(file.path(lib, p))) {
    ok <- file.copy(src, lib, recursive = TRUE, copy.date = TRUE)
    cat(if (ok) "copied    :" else "FAILED    :", p, "\n")
  } else {
    cat("already   :", p, "\n")
  }
}

## 3. Verify the closure loads from the renv library
cat("\nVerifying...\n")
ok <- requireNamespace("GenomicSEM", quietly = TRUE) &&
  requireNamespace("svglite", quietly = TRUE)
cat(if (ok) "GenomicSEM + svglite load OK\n" else "STILL FAILING - paste this output\n")
