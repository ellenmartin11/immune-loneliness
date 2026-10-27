---
title: "AoU Observational Results"
output: html_document
---

# Examining Observational/EHR Associations between Loneliness and Psychiatric and Immune Conditions

## UCLA Loneliness Cohort Description and Sample Sizes
- N UCLA Loneliness Scale = 290,522

### N Immune Biomarker Measures with UCLA Loneliness:

| Immune Biomarker | N (within 1 year of UCLA Loneliness) |
|------------------|-------------------------------------:|
| CRP              | 16,466                               |
| Albumin          | 81,207                               |
| Basophils        | 43,382                               |
| Eosinophils      | 46,606                               |
| Leukocytes       | 75,849                               |
| Monocytes        | 45,451                               |
| Neutrophils      | 42,559                               |
| Platelets        | 80,823                               |
| Triglycerides    | 89,235                               |

### N Psychiatric Condition Diagnoses (consistent on at least 2 occasions) with UCLA Loneliness:

| Psychiatric Diagnosis (at least 2 separate classifications) | N (with UCLA Loneliness) |
|--------------------------------------------------------------|-------------------------:|
| Anxiety                                                      | 56,679                   |
| Depression                                                   | 52,388                   |
| PTSD                                                         | 9,922                    |
| ADHD                                                         | 8,756                    |
| Bipolar Disorder                                             | 6,691                    |
| Personality Disorder                                         | 2,646                    |
| Eating Disorder                                              | 2,488                    |
| OCD                                                          | 1,317                    |
| Schizophrenia                                                | 683                      |
| Autism                                                       | 620                      |

### N Immune Condition Diagnoses with UCLA Loneliness:

| Immune Condition      | N      |
|-----------------------|-------:|
| Asthma                | 32,055 |
| Allergic Rhinitis     | 42,911 |
| Rheumatoid Arthritis  | 5,420  |
| Hypothyroidism        | 30,177 |
| Ulcerative Colitis    | 2,198  |
| Type II Diabetes      | 29,183 |
| Crohn's Disease       | 2,111  |
| Atopic Dermatitis     | 5,815  |

### Demographic Covariates and Loneliness

- higher loneliness among females than males

![Loneliness by Sex](loneliness_sex (1).png)

***
- loneliness decreases as income increases
![Loneliness by Income](loneliness_income_plot.png)

## Loneliness and Psychiatric Diagnoses

- higher loneliness among individuals with any of these 10 diagnoses

![Loneliness Psych Diagnosis vs None](loneliness_any_psych_condition_barplot.png)

***

- highest loneliness among individuals with autism and personality disorder

![Loneliness Across Psych Diagnoses](loneliness_across_psych_barplot.png)

***

- loneliness (adjusting for covariates: race, sex, education, income, age, bmi) increases OR of psychiatric conditions

![Loneliness Psychiatric OR](loneliness_psych_OR.png)

***

### Loneliness and Immune Conditions

- statistically detectable but clinically not meaningful difference between individuals with immune condition diagnosis vs without

![Loneliness by Immune Condition](loneliness_any_immune_condition_barplot.png)

***

- loneliness increase OR of rheumatoid arthritis and other conditions adjusting for covariates + bmi

![Loneliness Immune Condition OR](loneliness_immune_condition_ODDS_forest (1).png)

***

### Loneliness and Immune Biomarkers

- statistically detectable but small effects of loneliness on immune biomarker counts

![Loneliness Immune Biomarker Correlation Heatmap](loneliness_immune_biomarker_correlation_heatmap(1).png)

***

- loneliness associated with increased counts of specific biomarkers, adjusting for covariates (also replicated using log-transformed values)
- adjusting for BMI, only neutrophils becomes no longer significant, but the others are attenuated down

![Loneliness Immune Biomarker Forest](loneliness_immune_biomarkers_standardized_forest.png)

***

## Loneliness, Income x Sex Interacions in Biomarker Levels


- For platelets, significant main effects of loneliness, sex and income, but no interactions
- Main effect of income on all biomarkers
- Main effect of sex on all biomarkers (but weak for neutrophils)
- Main effect of loneliness on all biomarkers except Albumin (males only)
- Sex x Income interaction for Basophils, CRP, Eosinophils, Monocytes & Triglycerides
- Sex x Loneliness interaction for Albumin, Leukocytes, Monocytes, Neutrophils & Triglycerides
- No sig. Loneliness x Income interactions
- Plotting Estimated Marginal means from interaction models

![Associations between loneliness and biomarker across income levels](biomarkers_loneliness_income/loneliness_immune_biomarker_income_sex.png)

***

### Albumin
- Loneliness associated with albumin levels only among males, not females
- Significant loneliness x male sex interaction (b = 0.063, p < 0.0001)

![Association between loneliness and albumin levels among males vs females](biomarkers_loneliness_income/albumin_loneliness_emmeans_sex.png)

***

### Basophils
- Loneliness associated with increased basophil levels.
- Male sex associated with higher basophil levels than females, but sex-specific differences reduce as income increases.
- Significant income x male sex interaction (b = -0.013, p = 0.008)

![Sex-differences basophils reduce with increased income](biomarkers_loneliness_income/basophils_loneliness_emmeans_income.png)

***

### CRP
- Loneliness associated with reduced CRP levels.
- But this pattern is reversed among low income males (relative to high income males, and females)
- Significant income x male sex interaction (b = -0.038, p < 0.001)

![Reversed effect of loneliness on CRP among low income males](biomarkers_loneliness_income/crp_loneliness_emmeans_income.png)

***

### Eosinophils
- Loneliness associated with reduced eosinophils
- Significant main effect of income on eosinophil level (b = -0.022, p < 0.001)

![Loneliness reduces eosinophil](biomarkers_loneliness_income/eosinophils_loneliness_emmeans_income.png)

***

### Leukocytes and Neutrophils
- Loneliness associated with increased leukocytes and neutrophils
- Significant loneliness x male sex interaction for leukocytes (b = -0.040, p < 0.001) and neutrophils (b = -0.029, p = 0.005)

![Loneliness x male sex interaction in Leukocytes](biomarkers_loneliness_income/leukocyte_loneliness_emmeans_sex.png)
![Loneliness x male sex interaction in Neutrophils](biomarkers_loneliness_income/neutrophils_loneliness_emmeans_sex.png)

***

### Monocytes
- Loneliness associated with increased monocyte count (slightly).
- There are opposing effects between males and females, where loneliness especially increases monocytes among females, but decreases it for males
- Significant loneliness x male sex interaction (b = -0.043, p < 0.001)

![Opposinge effects of loneliness on biomarker levels between males and females](biomarkers_loneliness_income/monocytes_loneliness_emmeans_income.png)

***

### Triglycerides
- Loneliness associated with increased triglycerides, particularly among males, and among males with higher income. 
- Significant loneliness x male sex interaction (b = 0.022, p = 0.003) and sex x income interaction (b = 0.010, p = 0.003)

![Loneliness associated with increased triglycerides, decreased with income](biomarkers_loneliness_income/triglycerides_loneliness_emmeans_income.png)








