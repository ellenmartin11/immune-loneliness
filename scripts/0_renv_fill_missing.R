############################################################################
### renv_fill_missing.R
###
### Copies packages that renv::snapshot() flagged as missing (they live in
### the R module's site/recommended library) into the project's renv R-4.4
### library, so the lockfile can pin the full dependency closure.
###
### Run from a directory OUTSIDE the immune project (e.g. your home dir),
### under R/4.4.1, so renv does not activate and find.package() can see the
### module libraries:
###
###   cd "$HOME"
###   module load R/4.4.1-foss-2022b
###   Rscript renv_fill_missing.R
###
### Then, inside the project:
###   Rscript -e 'renv::snapshot()'
############################################################################

pk <- c(
  # flagged by renv::snapshot() pre-flight
  "class", "crayon", "doParallel", "foreach", "gridBase", "gtools",
  "iterators", "MASS", "Matrix", "numDeriv", "pbivnorm", "pkgconfig",
  "plyr", "progress", "quadprog", "R.methodsS3", "tidyselect",
  # Matrix dependency + Rmd parsing (yaml warning)
  "lattice", "yaml"
)

lib <- "/gpfs/gibbs/project/polimanti/ejm96/immune/renv/library/linux-rhel-8.10/R-4.4/x86_64-pc-linux-gnu"
dir.create(lib, showWarnings = FALSE, recursive = TRUE)

cat("R:", R.version.string, "\n")
cat(".libPaths searched:\n"); print(.libPaths())

for (p in pk) {
  src <- tryCatch(find.package(p), error = function(e) NULL)
  if (is.null(src)) {
    cat("NOT FOUND :", p, "- install it manually\n")
    next
  }
  if (normalizePath(dirname(src)) == normalizePath(lib)) {
    cat("already   :", p, "\n")
    next
  }
  tgt <- file.path(lib, p)
  if (dir.exists(tgt)) unlink(tgt, recursive = TRUE)   # replace stale copies
  ok <- file.copy(src, lib, recursive = TRUE, copy.date = TRUE)
  cat(if (ok) "copied    :" else "FAILED    :", p, "from", dirname(src), "\n")
}

cat("\nDone. Now run inside the project:\n",
    "  Rscript -e 'renv::snapshot()'\n")
