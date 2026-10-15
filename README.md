# Loneliness and Immune Biomarkers: Genomic Analyses

## Overview

This repository contains code and analysis workflows for a project investigating the relationship between loneliness, social isolation, and immune-related biomarkers using genomic and population-scale health data.

The project aims to characterize shared genetic architecture between loneliness-related traits and immune phenotypes, identify potential causal relationships, and evaluate biological pathways that may connect social and immune processes.

## Research Questions

The main questions include:

1. Do loneliness and social isolation share genetic liability with immune-related biomarkers?
2. Which immune traits show the strongest genetic overlap with loneliness-related phenotypes?
3. Are associations consistent across different cohorts and phenotype definitions?
4. Is there evidence for directional or causal relationships between loneliness and immune biomarkers?
5. To what extent may behavioral, clinical, or social factors contribute to observed associations?

## Phenotypes

### Social and Behavioral Traits

Phenotypes of interest may include:

- Loneliness
- Social isolation
- Social connection and related psychosocial traits
- Behavioral and health-related covariates where relevant

### Immune and Inflammatory Traits

Immune-related phenotypes may include:

- White blood cell count
- Basophil count and/or percentage
- Other leukocyte subtypes
- C-reactive protein (CRP)
- Additional inflammatory and immune biomarkers identified through prior literature

## Data Sources

Analyses may use genome-wide association summary statistics and/or individual-level data from resources including:

- UK Biobank
- FinnGen
- All of Us Research Program
- Published GWAS datasets and consortia

Raw and controlled-access datasets are **not included in this repository**.

Only analysis scripts, workflow files, documentation, and non-sensitive metadata required for reproducibility are version controlled.

## Analytic Approach

Depending on the phenotype and available data, analyses may include:

- GWAS summary-statistic harmonization and quality control
- SNP-based heritability estimation
- Genetic correlation using LD Score Regression (LDSC)
- Cross-cohort phenotype comparison and replication
- Mendelian randomization
- Conditional and joint genetic analyses
- Genomic structural equation modeling where appropriate
- Sensitivity analyses across phenotype definitions and cohorts

Additional observational analyses may examine relationships between immune biomarkers, loneliness, social determinants of health, and relevant clinical or behavioral covariates.

## Repository Structure

```text
.
├── scripts/          # R, Python, and shell analysis scripts
├── src/              # Reusable functions and analysis utilities
├── config/           # Configuration files and phenotype definitions
├── docs/             # Project documentation and analysis notes
├── README.md
└── .gitignore
```

Large datasets and analysis outputs are stored separately and excluded from version control.

Typical excluded directories include:

```text
data/
sumstats/
results/
output/
logs/
scratch/
```

## Reproducibility

Where possible, scripts are organized to separate:

1. Data acquisition
2. Data formatting and harmonization
3. Quality control
4. Statistical analysis
5. Sensitivity analyses
6. Result generation

Analysis-specific software versions, reference panels, and computational requirements will be documented as the project develops.

## Computing Environment

Large-scale genomic analyses are performed using the Yale University high-performance computing environment.

Common tools used in this project may include:

- R
- Python
- PLINK
- LDSC
- GenomicSEM
- METAL
- Mendelian randomization packages
- SLURM for job scheduling

## Data Privacy and Access

This repository does not contain individual-level participant data, controlled-access genomic data, or protected health information.

Access to underlying datasets is governed by the relevant data-use agreements and institutional requirements for each resource.

## Project Status

This project is under active development. Analytic workflows, phenotype definitions, and documentation may change as analyses progress.

## Author

Ellen Martin  
Polimanti Lab  
Yale University
