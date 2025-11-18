# T2D_trtment_prog
This repository contains R code to accompany research titled 'Impact of depression on treatment progression in type 2 diabetes: A UK retrospective cohort study using the Clinical Practice Research Datalink Aurum database'

Abbreviations: CPRD = Clinical Practice Research Datalink; T2D = type 2 diabetes, HES-APC = Hospital Episode Statistics Admitted Patient Care, GLM = glucose-lowering medication, CKD = Chronic Kidney Disease.

---

## Data availability
This study used data from the **Clinical Practice Research Datalink (CPRD) Aurum**, obtained under license from the UK Medicines and Healthcare products Regulatory Agency (MHRA).  
These data cannot be shared publicly. Researchers may apply for access through CPRD: [https://www.cprd.com](https://www.cprd.com).

---

## Initial data set up
The raw CPRD Aurum data were processed by the [Exeter Diabetes Research Team](https://github.com/Exeter-Diabetes).
MySQL tables were created, and subsequent data management was performed using their [CPRD analysis R package](https://github.com/Exeter-Diabetes/CPRD-analysis-package).  

For further details, please see their GitHub repository or contact **Dr Katie Young** (k.young3@exeter.ac.uk).

---

## Depression codes folder
This folder contains two files relating to the depression codes used in this analysis.

1. 'depression_broad_medcodes.txt': this file contains CPRD medcode IDs, the term attached to the code and a category for the code (diagnosis, remission, partial_remission, recurrent_diag, review_invite, admin, support). This file needs to be added to the code sets table (see this [tutorial](https://exeter-diabetes.github.io/CPRD-Tutorial/) for details on setting up CPRD data for analysis).

2. 'depression_broad_medcode_terms.txt': this file contains the medcode IDs and terms again, but is used by study-specific scripts below.

--

## C1. Data Management

Data management followed the structure and code provided by the [Exeter Diabetes CPRD-Cohort-scripts](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/tree/Oct2020-download) repository, with the index date defined as the **first prescription of oral glucose-lowering monotherapy**.

### Workflow overview

To reproduce the data management workflow:

1. **Create a table with patient end dates** using the script provided here:
   - `01_create_end_dat.mysql` — creates a table called `all_patid_death_end_dat`. The scripts provided by the Exeter Research team used ONS death data linkage, and therefore had death date information from this source. This study did not have access to this linkage, and we therefore used the CPRD provided death date (cprd_ddate) as an alternative. This script creates a table that is used in the Exeter scripts instead of validDateLookup when death date is needed (e.g. column gp_death_end_date in all_patid_death_end_dat is used instead of column gp_ons_end_date in validDateLookup). This is the only script written in sql.

2. **Run the following Exeter CPRD cohort creation scripts**, as described in their repository:
   - [`all_diabetes_cohort.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/all_diabetes_cohort.R): identifies individuals with type 2 diabetes (excluding other diabetes types such as type 1).  
   - [`all_patid_ethnicity.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/all_patid_ethnicity.R): assigns ethnicity using the Exeter algorithm.  
   - [`all_patid_ckd_stages.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/all_patid_ckd_stages.R): identifies the first recorded CKD stage (1–5) in CPRD Aurum with HES-APC linkage.  

3. **Define the index date and outcome variables** using the script provided here:  
   - `02_define_index_outcomes.R` — sets the index date to the first oral GLM prescription, and creates the outcome variables (a. time to treatment intensification and it's corresponding binary event variable and b. time to insulin initiation and it's corresponding binary event variable). 

4. **Run the following Exeter scripts for variable derivation and cleaning** using the index date from item 3 and the patient end of study date from item 2:
   - [`template_alcohol.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/template_alcohol.R): Extract alcohol codes and defined alcohol status at index. 
   - [`template_smoking.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/template_smoking.R): Extract smoking codes and defined smoking status at index.
   - [`template_baseline_biomarkers.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/template_baseline_biomarkers.R): Extract and clean biomarkers, find values at index. Note. We XXX {update HbA1c?}
   - [`template_ckd_stages.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/template_ckd_stages.R): Define CKD stages at index.
   - [`template_comorbidities.R`](https://github.com/Exeter-Diabetes/CPRD-Cohort-scripts/blob/Oct2020-download/template_comorbidities.R): Extract and clean T2D comorbidities codes and identify diagnoses/ events at index.

5. **Extract preliminary depression data at index** using the scripts provided in this repository (see: `03_depression_at_index.R`).  

6. **Create the final analysis dataset** using:  
   - `04_create_final_analysis_dataset.R`


### Notes
- These scripts are intended to be run within the CPRD data environment after MySQL table setup using the [Exeter CPRD analysis package](https://github.com/Exeter-Diabetes/CPRD-analysis-package).  
- The CPRD Aurum data used in this project are not publicly available.  
- Only scripts unique to this study (i.e., those defining the index date and creating the final analysis dataset) are included in this repository.  

--

## C2. Primary Analysis

This study includes two time-to-event outcomes:  
1. **Time to insulin initiation**, and  
2. **Time to treatment intensification**.

The primary analyses for both outcomes were conducted using scripts that *perform multiple imputation and then fit the survival models* (Royston–Parmar proportional odds and proportional hazards models). These analyses were carried out separately for each outcome:

1. **Time to insulin initiation**  
   - `01_ttinsulin_MI_survival.R`  
     *Performs multiple imputation of covariates and fits the survival models for insulin initiation.*

2. **Time to treatment intensification**  
   - `02_ttintens_MI_survival.R`  
     *Performs multiple imputation of covariates and fits the survival models for treatment intensification.*

### C2.1 Standardised Survival, Risk Differences, and RMST

Following model fitting, standardised estimates were generated to facilitate interpretation of the depression–treatment progression associations. This step used the fitted models from each imputed dataset to obtain standardised survival curves, risk differences, and restricted mean survival times (RMST):

- **Time to insulin initiation**
  - `03_ttinsulin_standardised_surv.R`  
    *Generates covariate-standardised survival curves for insulin initiation across depression groups, computes standardised risk differences at prespecified time points, and estimates RMST and RMST differences. Results are pooled across imputations using Rubin’s rules.*

- **Time to treatment intensification**
  - `04_ttintens_standardised_surv.R`  
    *Produces the equivalent standardised survival, risk difference, and RMST outputs for the treatment intensification outcome.*

### C2.2 Visualisation of Results

A set of scripts were used to produce the figures included in the manuscript. These scripts take the processed outputs from the multiple imputation models and standardised survival analyses and generate plots used for interpretation and reporting.

The following figures were produced:

- **Kaplan–Meier curves**
  - `05_kmplots.R`  
    *Plots unadjusted Kaplan–Meier curves for each depression group.*

- **Adjusted effect estimates (ORs / HRs with confidence intervals)**  
  - `06_forestplot_OR.R`   
    *Produces forest-style plots of model-based estimates (odds ratios or hazard ratios) with 95% CIs.*

- **Standardised survival curves and absolute risk differences**  
  - `07_standardised_survival_plots.R`  
    *Plots covariate-standardised survival curves and pooled risk differences derived from the Royston–Parmar models.*

--

## C3. Sensitivity Analyses

To be added: Complete case. Then extras.

--

## Notes
- Scripts are provided for transparency; they cannot be executed without access to CPRD Aurum data. 
- Analyses were performed in R (v4.3.2).
