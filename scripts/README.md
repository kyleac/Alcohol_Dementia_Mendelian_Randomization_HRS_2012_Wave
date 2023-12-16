# Relationship between alcohol consumption and dementia with Mendelian randomization approaches among older adults in the United States

## Citation Information
Ware EB, Bustamante ACM, Fu M, Bakulski KM. Type 2 diabetes and dementia in the Health and Retirement Study: A Mendelian randomization approach. Alzheimer’s & Dementia. 2020;16(S10):e041220. doi:https://doi.org/10.1002/alz.041220

This Github repository contains the data management and analytic scripts to produce the following manuscript:[Type 2 diabetes and dementia in the Health and Retirement Study: A Mendelian randomization approach. Alzheimer’s & Dementia](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8044888/)

## Abstract
Background: In observational studies, the association between alcohol consumption and dementia is mixed.
Methods: We performed two-sample Mendelian randomization (MR) using summary statistics from genome-wide association studies of weekly alcohol consumption and late-onset Alzheimer’s disease and one-sample MR in the Health and Retirement Study (HRS), wave 2012. Inverse variance weighted two-stage regression provided odds ratios of association between alcohol exposure and dementia or cognitively impaired, non-dementia relative to cognitively normal.
Results: Alcohol consumption was not associated with late-onset Alzheimer’s disease using two-sample MR (OR=1.15, 95% confidence interval (CI):[0.78, 1.72]). In HRS, doubling weekly alcohol consumption was not associated with dementia (African ancestries, n=1,322, OR=1.00, 95% CI [0.45, 2.25]; European ancestries, n=7,160, OR=1.37, 95% CI [0.53, 3.51]) or cognitively impaired, non-dementia (African ancestries, n=1,322, OR=1.17, 95% CI [0.69, 1.98]; European ancestries, n=7,160, OR=0.75, 95% CI [0.47, 1.22]).
Conclusion: Alcohol consumption was not associated with cognitively impaired, non-dementia or dementia status.

## Declarations
We thank the participants and staff of the Health and Retirement Study. This research was supported by the National Institutes of Heath: National Institute on Aging (R01 AG055406, R01 AG067592).

## Data availability
Publicly available datasets were analyzed in this study. This data can be found here: Health and Retirement Study health and covariate data are available here: https://hrs.isr.umich.edu/data-products. Health and Retirement Study genetic data are available here: dbGaP Study Accession: phs000428.v2.p2 (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000428.v2.p2).

## Script Files

alc_dem_sample_descriptions.Rmd: descriptive statistics and tables

alcohol_GSCAN_2019_LOAD_Kunkle_2019_two_sample_mr.Rmd: contains the two-sample Mendelian randomization analysis

gscan_dpw_1smr_xs_data_prep.Rmd: data prep for one-sample Mendelian randomization analysis, including inclusion/exclusion

gscan_dpw_1smr_xs_assumptions.Rmd: testing the Mendelian randomization assumptions

gscan_dpw_1smr_xs_model.Rmd: contains the actual one-sample Mendelian randomization analysis and conventional logistic regression analysis

fit_models.R: contains custom functions for model fitting

fit_relevance_models.R: contains custom functions for some Mendelian randomization assumptions testing

gscan_dpw_1smr_xs_model_functions.R: contains custom functions for one-sample Mendelian randomization analysis

formatting_functions.r: contains custom functions for data formatting and labelling


