Script started on 2026-11-04 08:51:39-05:00 [TERM="xterm-256color" TTY="/dev/pts/0" COLUMNS="193" LINES="18"]
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ ^C
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Error in library(GenomicSEM) : there is no package called ‘GenomicSEM’
[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[5P- The project is out-of-sync -- use `renv::status()` for details.
[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[CError in library(GenomicSEM) : there is no package called ‘GenomicSEM’
[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[C[KRscript "scripts/7_make_n_h2_table.R" "/scripts/cluster scripts/metal_ldsc_pipeline/ldsr/munge_manifest_poster.tsv" "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pi
ipeline/ldsr_h2/h2_summary.tsv" "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/group2_biomarker_traits.txt" "Immune Biomarkers"
The following package(s) were loaded before renv activated this project:
- BiocManager [~/.cache/R/renv/cache/v5/R-4.2/x86_64-pc-linux-gnu/BiocManager/1.30.19/726b3bc83ff8e0c35c7efdb74a462b0d/BiocManager]
These packages are not managed by renv and may cause `]8;;x-r-run:renv::status()renv::status()]8;;`
and `]8;;x-r-run:renv::snapshot()renv::snapshot()]8;;` to report inconsistencies. Restart R in a clean
session to recover. The paths above may help identify where each
package was loaded from.
Set `RENV_CONFIG_NAMESPACES_CHECK=FALSE` to silence this check.

- The project is out-of-sync -- use `]8;;x-r-run:renv::status()renv::status()]8;;` for details.
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ 
[K[ejm96@r209u11n01 immune]$ 
[K[ejm96@r209u11n01 immune]$ 
[K[ejm96@r209u11n01 immune]$ ^C
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Sys.setenv(RENV_CONFIG_SYSREQS_CHECK = "FALSE")
bash: syntax error near unexpected token `RENV_CONFIG_SYSREQS_CHECK'
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ # Now rerun the same renv install/restore command
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Rscript "scripts/plot_rg_heatmap.R"" "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_matrix.tsv" "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_
_pipeline/ldsr_rg/group2_biomarker_traits.txt" fig_rg_biomarker "/vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/ldsr_rg/rg_long.tsv" "Immune Biomarkers"
> 
> ^C
]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Rscript "/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/6_plot_rg_heatmap.R" /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_matrix.tsv "/vast/
/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/group2_biomarker_traits.txt" fig_rg_biomarker /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_long.tsv "Immune Bio
omarkers"
The following package(s) were loaded before renv activated this project:
- BiocManager [~/.cache/R/renv/cache/v5/R-4.2/x86_64-pc-linux-gnu/BiocManager/1.30.19/726b3bc83ff8e0c35c7efdb74a462b0d/BiocManager]
These packages are not managed by renv and may cause `]8;;x-r-run:renv::status()renv::status()]8;;`
and `]8;;x-r-run:renv::snapshot()renv::snapshot()]8;;` to report inconsistencies. Restart R in a clean
session to recover. The paths above may help identify where each
package was loaded from.
Set `RENV_CONFIG_NAMESPACES_CHECK=FALSE` to silence this check.

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25hNOTE: not in the rg matrix (skipped): Tyrosine
[?25h[?25h[?25hError in hclust(as.dist(1 - abs(sub))) : 
  NA/NaN/Inf in foreign function call (arg 10)
Execution halted
[?25h]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Rscript "/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/6_plot_rg_heatmap.R" /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_matrix.tsv "/vast/
/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/group2_biomarker_traits.txt" fig_rg_biomarker /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_long.tsv "Immune Bio
omarkers"
The following package(s) were loaded before renv activated this project:
- BiocManager [~/.cache/R/renv/cache/v5/R-4.2/x86_64-pc-linux-gnu/BiocManager/1.30.19/726b3bc83ff8e0c35c7efdb74a462b0d/BiocManager]
These packages are not managed by renv and may cause `]8;;x-r-run:renv::status()renv::status()]8;;`
and `]8;;x-r-run:renv::snapshot()renv::snapshot()]8;;` to report inconsistencies. Restart R in a clean
session to recover. The paths above may help identify where each
package was loaded from.
Set `RENV_CONFIG_NAMESPACES_CHECK=FALSE` to silence this check.

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25hError in hclust(as.dist(1 - abs(sub))) : 
  NA/NaN/Inf in foreign function call (arg 10)
Execution halted
[?25h]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Rscript "/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/6_plot_rg_heatmap.R" /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_matrix.tsv "/vast/
/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/group2_biomarker_traits.txt" fig_rg_biomarker /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_long.tsv "Immune Bio
omarkers"
The following package(s) were loaded before renv activated this project:
- BiocManager [~/.cache/R/renv/cache/v5/R-4.2/x86_64-pc-linux-gnu/BiocManager/1.30.19/726b3bc83ff8e0c35c7efdb74a462b0d/BiocManager]
These packages are not managed by renv and may cause `]8;;x-r-run:renv::status()renv::status()]8;;`
and `]8;;x-r-run:renv::snapshot()renv::snapshot()]8;;` to report inconsistencies. Restart R in a clean
session to recover. The paths above may help identify where each
package was loaded from.
Set `RENV_CONFIG_NAMESPACES_CHECK=FALSE` to silence this check.

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25hNOTE: NA rg values - clustering skipped, file order used. NA pairs: Glycoprotein acetylation-Basophil Count PanUKB, Glycoprotein acetylation-Monocytes, Glycoprotein acetylation-Neutrophils, Basophil Count PanUKB-Glycoprotein acetylation, Monocytes-Glycoprotein acetylation, Neutrophils-Glycoprotein acetylation
[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25hWritten: fig_rg_biomarker.svg and fig_rg_biomarker.png 
[?25h[?25h]0;/gpfs/gibbs/project/polimanti/ejm96/immune[ejm96@r209u11n01 immune]$ Rscript "/gpfs/gibbs/project/polimanti/ejm96/immune/scripts/6_plot_rg_heatmap.R" /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_matrix.tsv "/vast/
/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/group2_biomarker_traits.txt" fig_rg_biomarker /vast/palmer/scratch/polimanti/ejm96/metal_ldsc_pipeline/ldsr_rg/rg_long.tsv "Immune Bio
omarkers"
The following package(s) were loaded before renv activated this project:
- BiocManager [~/.cache/R/renv/cache/v5/R-4.2/x86_64-pc-linux-gnu/BiocManager/1.30.19/726b3bc83ff8e0c35c7efdb74a462b0d/BiocManager]
These packages are not managed by renv and may cause `]8;;x-r-run:renv::status()renv::status()]8;;`
and `]8;;x-r-run:renv::snapshot()renv::snapshot()]8;;` to report inconsistencies. Restart R in a clean
session to recover. The paths above may help identify where each
package was loaded from.
Set `RENV_CONFIG_NAMESPACES_CHECK=FALSE` to silence this check.

[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h[?25h