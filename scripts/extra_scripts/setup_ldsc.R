# ============================================================
# LDSC environment setup
# ============================================================

# Paths
ldsc_venv <- "/gpfs/gibbs/project/polimanti/ejm96/project/ldsc_env"
ldsc_python <- file.path(ldsc_venv, "bin", "python")

ldsc_dir <- "/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-python3"
cb_ldsc_dir <- "/gpfs/gibbs/project/polimanti/ejm96/project/ldsc-CBIIT"

python310_lib <- "/vast/palmer/apps/avx2/software/Python/3.10.8-GCCcore-12.2.0/lib"
libffi_lib <- "/vast/palmer/apps/avx2/software/libffi/3.4.4-GCCcore-12.2.0/lib64"

# Remove inherited Python 3.8 package paths
Sys.unsetenv("PYTHONPATH")
Sys.unsetenv("PYTHONHOME")

# Preserve existing LD_LIBRARY_PATH
current_ld <- Sys.getenv("LD_LIBRARY_PATH")

# Add libraries needed by Python 3.10
Sys.setenv(
  LD_LIBRARY_PATH = paste(
    python310_lib,
    libffi_lib,
    current_ld,
    sep = ":"
  )
)

# ============================================================
# Sanity checks
# ============================================================

cat("\n--- LDSC setup check ---\n")

cat("\nPython executable:\n")
print(ldsc_python)

if (!file.exists(ldsc_python)) {
  stop("LDSC Python executable not found.")
}

cat("\nPython version:\n")
system2(
  ldsc_python,
  args = "--version"
)

cat("\nPython / NumPy check:\n")
status_numpy <- system2(
  ldsc_python,
  args = c(
    "-c",
    shQuote(
      paste(
        "import sys, numpy;",
        "print('Python:', sys.executable);",
        "print('NumPy:', numpy.__version__);",
        "print('NumPy path:', numpy.__file__)"
      )
    )
  )
)

if (!identical(status_numpy, 0L)) {
  stop("Python/NumPy environment check failed.")
}

cat("\nCBIIT LDSC check:\n")
status_ldsc <- system2(
  ldsc_python,
  args = c(
    file.path(cb_ldsc_dir, "ldsc.py"),
    "-h"
  ),
  stdout = FALSE,
  stderr = FALSE
)

if (!identical(status_ldsc, 0L)) {
  stop("CBIIT LDSC failed to start.")
}

cat("\nEnvironment variables:\n")
cat("PYTHONPATH =", Sys.getenv("PYTHONPATH"), "\n")
cat("PYTHONHOME =", Sys.getenv("PYTHONHOME"), "\n")

cat("\nLDSC environment successfully initialized.\n")
cat("Use:\n")
cat("  ldsc_python =", ldsc_python, "\n")
cat("  cb_ldsc_dir =", cb_ldsc_dir, "\n")
cat("----------------------------------------\n")