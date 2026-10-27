# make_missing_traits.R
suppressPackageStartupMessages(library(data.table))

url <- "https://docs.google.com/spreadsheets/d/1tJhUVb1uH2Eku4P_pMAYwJCagoFiO_nYxyqmdHUgAFk/export?format=csv"

ref <- fread(url, na.strings = c("", "NA", "N/A"))
cat("rows read:", nrow(ref), "
")
cat("columns:", paste(names(ref), collapse = " | "), "

")

cat("analysis values:
")
print(table(trimws(ref$analysis), useNA = "ifany"))

missing <- c("Allergic rhinitis FinnGen", "Hypertension FinnGen", "Rheumatoid Arthritis FinnGen",
             "Ulcerative colitis FinnGen")

ph <- trimws(ref$phenotype)
cat("
phenotypes matching the missing list:
")
print(ph[ph %in% missing])

cat("
near-matches (fuzzy) for anything not matched exactly:
")
for (m in missing[!missing %in% ph]) {
  hit <- agrep(m, ph, value = TRUE, max.distance = 2)
  cat(" -", m, "->", if (length(hit)) paste(hit, collapse = "; ") else "NO MATCH", "
")
}

sub <- ref[tolower(trimws(analysis)) == "primary" & ph %in% missing]
cat("
selected:", nrow(sub), "traits
")
fwrite(sub, "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/primary_clean/missing_traits.csv")
