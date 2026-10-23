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