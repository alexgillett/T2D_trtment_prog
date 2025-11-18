############################################################################################
############################################################################################
### Create analysis dataset
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
# Load R libraries, link to your CPRD data via the aurum package, load in existing tables needed in script
########################################
### Section 2: Merge basic covariate information with initial analysis dataset
########################################
### Section 3: Merge with initial depression data
########################################
### Section 4: Merge with smoking, alcohol, ethnicity, biomarkers and comorbidities
########################################
### Section 5: Restrict to participants with >= 2 HbA1c >= 48 mmol/mol
########################################
### Section 6: Merge with IMD
########################################
### Section 7: Create preliminary healthcare utilisation variable(s) and merge
########################################
### Section 8: Provide non-numeric levels for (some) categorical variables
########################################
### Section 9: Collect data and reduce to those initiating medication before 2022
########################################
### Section 10: Create exposure variable: recency depression 
########################################
### Section 11: Define and add in CKD stage 3 at index
########################################
### Section 12: Extra covariate creation
########################################
### Section 13: Reduce to columns needed for multiple imputation analysis
########################################
### Section 14: Store dataset in secure approved location
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "~/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

cohort_prefix <- ""
# e.g. "mm" for treatment response (MASTERMIND) cohort
analysis_prefix <- ""
# e.g. "at1OAD" for analysis indexing at single oral antidiabetic medication

# Assign name to existing tables to be used

# T2D cohort
analysis = cprd$analysis(cohort_prefix)
t2_cohort <- t2_cohort %>% analysis$cached("t2_cohort")
# death data
analysis = cprd$analysis("all_patid")
death_end_dat <- death_end_dat %>% analysis$cached("death_end_dat")
# index dates
analysis = cprd$analysis(analysis_prefix)
index_dates <- index_dates1oad %>% analysis$cached("index_dates1oad")
# Clean HbA1c measures
analysis = cprd$analysis("t2_cohort")
hba1c <- clean_hba1c_medcodes %>% analysis$cached("clean_hba1c_medcodes")

# Initial time-to-progression analysis dataset
analysis = cprd$analysis(analysis_prefix)
first_change <- ttprogression_v1 %>% analysis$cached("ttprogression_v1")

# Depression data...
analysis = cprd$analysis(analysis_prefix)
pre1oad_date <- pre1oad_lastdate_v3 %>% analysis$cached("pre1oad_lastdate_v3")

### ethncity
analysis = cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity")
### smoking and alcohol
analysis = cprd$analysis(analysis_prefix)
smoking <- smoking %>% analysis$cached("smoking")
analysis = cprd$analysis(analysis_prefix)
alcohol <- alcohol %>% analysis$cached("alcohol")
### biomarkers and comorbidities
analysis = cprd$analysis(analysis_prefix)
baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")
comorbidities <- comorbidities %>% analysis$cached("comorbidities")

############################################################################################
### Section 2: Merge basic covariate information with initial analysis dataset
############################################################################################
### Extract relevant info such as gender, age at T2D diagnosis, T2D diagnosis date
t2_cohort2 <- t2_cohort %>%
    select(patid, gender, dob, pracid, prac_region, dm_diag_date, dm_diag_age, regstartdate, regenddate, lcd)

analysis = cprd$analysis(analysis_prefix)
med_change_cohort <- first_change %>%
    left_join(t2_cohort2) %>%
    mutate(
        age1oad = (datediff(first_start, dob))/365.25,
        disease_dur1oad_yrs = (datediff(first_start, dm_diag_date))/365.25
    ) %>% 
    analysis$cached("cohort1oad_interim", indexes=c("patid", "gender", "dob", "pracid", "dm_diag_date", "dm_diag_age", "gp_death_end_date","first_start",
        "date_ins_cens", "date_change"))
    
############################################################################################
### Section 3: Merge with depression data
############################################################################################

med_change_cohort2 <- med_change_cohort %>%
  left_join(pre1oad_date) %>%
  mutate(depr_yn = ifelse(is.na(depr_yn), 0, depr_yn)) %>%
  mutate(sum_dep_tmp = predep10y + predep5y + predep2y + predep1y + predep6m) %>%
    analysis$cached("cohort1oad_interim2", indexes=c("patid", "gender", "dob", "pracid", "dm_diag_date", "dm_diag_age", "gp_death_end_date", "first_start",
        "date_ins_cens", "date_change"))

############################################################################################
### Section 4: Merge with smoking, alcohol, ethnicity, biomarkers and comorbidities
############################################################################################

cohort2011 <- med_change_cohort2 %>%
    mutate(year1oad = year(first_start)) %>%
    filter(year1oad >= 2011)

analysis = cprd$analysis(analysis_prefix)
cohort2011 <- cohort2011 %>%
    left_join(ethnicity) %>%
    left_join(smoking) %>%
    left_join(alcohol)

cohort2011bio <- cohort2011 %>%
    left_join(baseline_biomarkers, by="patid") %>%
    analysis$cached("cohort2011_biomarkers", indexes = "patid")

cohort2011comorbs <- cohort2011bio %>%
    left_join(comorbidities, by="patid") %>%
    analysis$cached("cohort2011_bio_comorbs_int1", indexes = "patid")

############################################################################################
### Section 5: Restrict to participants with >= 2 HbA1c >= 48 mmol/mol
############################################################################################
### Find IDs where 2 hba1c values need to be >= 48 to be included

analysis = cprd$analysis(analysis_prefix)
hba1c_keep <- hba1c %>% 
    inner_join(index_dates, by="patid") %>%
    mutate(t2d_indicator = ifelse(testvalue >= 48, 1, 0)) %>%
    group_by(patid) %>%
    mutate(numt2dvalues = sum(t2d_indicator, na.rm=T)) %>%
    ungroup() %>%
    select(patid, numt2dvalues) %>% distinct() %>%
    filter(numt2dvalues >= 2) %>%
    select(patid) %>%
    analysis$cached("hba1c_t2d_ids", indexes = c("patid"))

### Restrict to only the IDs matching this HbA1c criteria
cohort2011comorbs <- cohort2011comorbs %>%
  inner_join(hba1c_t2d_ids) %>%
    analysis$cached("cohort2011_bio_comorbs_int2", indexes = "patid")

############################################################################################
### Section 6: Merge with IMD
############################################################################################
patient_imd <- cprd$tables$patientImd
practice_imd <- cprd$tables$practiceImd

patient_imd1 <- patient_imd %>% select(patid, imd_decile) %>% rename(imd_decile_patient = imd_decile)

cohort2011comorbs <- cohort2011comorbs %>%
  left_join(patient_imd1) %>%
  left_join(practice_imd) %>%
  analysis$cached("cohort2011_bio_comorbs_int3", indexes = "patid")

############################################################################################
### Section 7: Create preliminary healthcare utilisation variables and merge
############################################################################################
### point to depression count dataset...
analysis = cprd$analysis(analysis_prefix)
pret2ddep_code1y <- pret2ddep_code1y %>% analysis$cached("pret2ddep_code1y")

### Unique GP visit dates in the year prior to 1OAD
analysis = cprd$analysis(analysis_prefix)
countGP <- cprd$tables$observation %>% 
  select(patid, obsdate) %>%
  distinct() %>%
  inner_join(index_dates) %>%
  filter(obsdate < index_date) %>%
  mutate(
    indexminus1yr = sql("DATE_SUB(index_date, INTERVAL 1 YEAR)")
  ) %>%
  filter(obsdate >= indexminus1yr) %>%
  group_by(patid) %>%
  summarise(
    countGP1yr = n()
  ) %>%
  ungroup() %>%
  analysis$cached("gpcount1yr1oad", indexes=c("patid"))

cohort2011comorbs <- cohort2011comorbs %>%
  left_join(countGP) %>%
  left_join(pret2ddep_code1y) %>%
  mutate(
    countGP1yr = ifelse(is.na(countGP1yr), 0, countGP1yr),
    n_dep_dates1y = ifelse(is.na(n_dep_dates1y), 0, n_dep_dates1y)
  ) %>% 

  analysis$cached("cohort2011_bio_comorbs_int4", indexes = "patid")

############################################################################################
### Section 8: Provide non-numeric levels for categorical variables
############################################################################################
analysis_ds <- cohort2011comorbs 

analysis_ds <- analysis_ds %>%
    mutate(
        sex = case_when(
            gender == 1 ~ "Male",
            gender == 2 ~ "Female",
            TRUE ~ NA_character_
        ),
        ethn_5cat_v1 = case_when(
            ethnicity_5cat == 0 ~ "White",
            ethnicity_5cat == 1 ~ "South Asian",
            ethnicity_5cat == 2 ~ "Black",
            ethnicity_5cat == 3 ~ "Other",
            ethnicity_5cat == 4 ~ "Mixed",
            TRUE ~ "Unknown"
        ),
        ethn_5cat_v2 = case_when(
            ethnicity_5cat == 0 ~ "White",
            ethnicity_5cat == 1 ~ "South Asian",
            ethnicity_5cat == 2 ~ "Black",
            ethnicity_5cat == 3 ~ "Other",
            ethnicity_5cat == 4 ~ "Mixed",
            ethnicity_5cat == NA ~ "White",
            TRUE ~ "White"
        ),
        ethn_qrisk2_v1 = case_when(
            ethnicity_qrisk2 == "1" ~ "White",
            ethnicity_qrisk2 == "2" ~ "Indian",
            ethnicity_qrisk2 == "3" ~ "Pakistani",
            ethnicity_qrisk2 == "4" ~ "Bangladeshi",
            ethnicity_qrisk2 == "5" ~ "Other Asian",
            ethnicity_qrisk2 == "6" ~ "Black Caribbean",
            ethnicity_qrisk2 == "7" ~ "Black African",
            ethnicity_qrisk2 == "8" ~ "Chinese",
            ethnicity_qrisk2 == "9" ~ "Other",
            ethnicity_qrisk2 == "0.0" ~ "Unknown",
            TRUE ~ "Unknown"
        ),
        smk_qrisk2_acg = case_when(
            qrisk2_smoking_cat_uncoded %in% c("Light smoker", "Moderate smoker", "Heavy smoker") ~ "Current smoker",
            TRUE ~ qrisk2_smoking_cat_uncoded
        ),
        smk_qrisk2_impv1 = case_when(
            is.na(qrisk2_smoking_cat_uncoded) ~ "Non-smoker",
            qrisk2_smoking_cat_uncoded %in% c("Light smoker", "Moderate smoker", "Heavy smoker") ~ "Current smoker",
            TRUE ~ qrisk2_smoking_cat_uncoded
        ),
        alcohol_cat_impv1 = case_when(
            is.na(alcohol_cat) ~ "Within limits",
            TRUE ~ alcohol_cat
        )
    ) %>% 
  analysis$cached("analysisds_withmissing", indexes = "patid")


############################################################################################
### Section 9: Collect data and reduce to those initiating medication before 2022
############################################################################################
med_ds1 <- analysis_ds %>% collect()
med_ds1 <- med_ds1 %>% 
    filter(first_start != gp_death_end_date) %>%
    filter(year1oad <= 2022)

############################################################################################
### Section 10: Create exposure variable: recency depression 
############################################################################################
q_dep <- quantile(med_ds1$lastpre1oaddep_time[med_ds1$depr_yn == 1])
q_dep
###     0%     25%     50%     75%    100% 
### 0.0000  1.6564  5.6126 12.8186 70.9514 
med_ds1 <- med_ds1 %>%
    mutate(
        lastpre1oaddep_quant = case_when(
            med_ds1$depr_yn == 0 ~ "No MDD",
            med_ds1$depr_yn == 1 & lastpre1oaddep_time <= q_dep[2] ~ "MDD <= 1.7yrs",
            med_ds1$depr_yn == 1 & lastpre1oaddep_time <= q_dep[4] & lastpre1oaddep_time > q_dep[2] ~ "MDD 1.7-12.8yrs",
            TRUE ~ "MDD > 12.8yrs"
        )
    )

med_ds1 <- med_ds1 %>%
    mutate(
        lastpre1oaddep_quant = factor(lastpre1oaddep_quant, levels = c("No MDD", "MDD <= 1.7yrs", "MDD 1.7-12.8yrs", "MDD > 12.8yrs"))

    )

############################################################################################
### Section 11: Define and add in CKD stage 3 at index
############################################################################################

# Get pointer to longitudinal CKD stage table

analysis = cprd$analysis("all_patid")
ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>% analysis$cached("ckd_stages_from_algorithm")
                  

# Merge with index dates to get CKD stages at index date

# Merge with CKD stages (1 row per patid)
analysis = cprd$analysis(analysis_prefix)
ckd_stage_drug_merge <- index_dates %>%
  left_join(ckd_stages_from_algorithm, by="patid") %>%
  mutate(preckdstage=ifelse(!is.na(stage_5) & datediff(stage_5, index_date)<=7, "stage_5",
                            ifelse(!is.na(stage_4) & datediff(stage_4, index_date)<=7, "stage_4",
                                   ifelse(!is.na(stage_3b) & datediff(stage_3b, index_date)<=7, "stage_3b",
                                          ifelse(!is.na(stage_3a) & datediff(stage_3a, index_date)<=7, "stage_3a",
                                                 ifelse(!is.na(stage_2) & datediff(stage_2, index_date)<=7, "stage_2",
                                                        ifelse(!is.na(stage_1) & datediff(stage_1, index_date)<=7, "stage_1", NA)))))),
         
         preckdstagedate=ifelse(preckdstage=="stage_5", stage_5,
                                ifelse(preckdstage=="stage_4", stage_4,
                                       ifelse(preckdstage=="stage_3b", stage_3b,
                                              ifelse(preckdstage=="stage_3a", stage_3a,
                                                     ifelse(preckdstage=="stage_2", stage_2,
                                                            ifelse(preckdstage=="stage_1", stage_1, NA)))))),
         
         preckdstagedatediff=datediff(preckdstagedate, index_date)) %>%
  
  select(patid, index_date, preckdstage, preckdstagedate, preckdstagedatediff) %>%
  
  analysis$cached("ckd_stages", indexes=c("patid", "index_date"))

analysis = cprd$analysis(analysis_prefix)                           
ckd3plus <- ckd_stage_drug_merge %>%
  mutate(
    ckd3plus = case_when(
      preckdstage %in% c("stage_1", "stage_2") ~ 0,
      preckdstage %in% c("stage_3a", "stage_3b", "stage_4", "stage_5") ~ 1,
      TRUE ~ NA_real_
    ),
    ckd3plus = ifelse(is.na(ckd3plus), 0, ckd3plus)   
  ) %>%
  select(patid, ckd3plus) %>%
  analysis$cached("ckd3plus", indexes=c("patid"))

ckd3 <- collect(ckd3plus)

med_ds1 <- med_ds1 %>% left_join(ckd3)

############################################################################################
### Section 12: Extra covariate creation
############################################################################################

med_ds1 <- med_ds1 %>%
    mutate(
        ethn_5cat_analysis = case_when(
            ethnicity_5cat == 0 ~ "White",
            ethnicity_5cat == 1 ~ "South Asian",
            ethnicity_5cat == 2 ~ "Black",
            ethnicity_5cat == 3 ~ "Other",
            ethnicity_5cat == 4 ~ "Mixed",
            TRUE ~ NA
        ),
        e2019_imd_decilef = as.factor(e2019_imd_decile),
        gpcount1yr_excdep = countGP1yr - n_dep_dates1y,
        n_diabcomps = pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + pre_index_date_neuropathy + ckd3plus +
        pre_index_date_primary_hhf + pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction +
        pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia
    )

med_ds1 <- med_ds1 %>%
    mutate(
      year1oad_c = year1oad - median(year1oad)
    )

############################################################################################
### Section 13: Reduce to columns needed for multiple imputation analysis
############################################################################################
med_ds1 <- med_ds1 %>% 
    select(
        intens_event, time_change_yrs, ins_ever, time_insulin_yrs, lastpre1oaddep_quant, 
        age1oad, sex, disease_dur1oad_yrs, ethn_5cat_analysis, prebmi, prehba1c, 
        presbp, predbp,
        smk_qrisk2_acg, alcohol_cat, e2019_imd_decilef, gpcount1yr_excdep, pre_index_date_retinopathy, pre_index_date_diabeticnephropathy, 
        pre_index_date_ckd5_code, pre_index_date_neuropathy, pre_index_date_hypertension, pre_index_date_ihd, pre_index_date_heartfailure, 
        pre_index_date_myocardialinfarction, pre_index_date_pad, pre_index_date_stroke, pre_index_date_tia, pre_index_date_cld, 
        pre_index_date_copd, pre_index_date_bronchiectasis, pre_index_date_asthma, pre_index_date_pulmonaryfibrosis, pre_index_date_pulmonaryhypertension, 
        pre_index_date_solid_cancer, pre_index_date_haem_cancer, pre_index_date_solidorgantransplant, pre_index_date_rheumatoidarthritis, 
        pre_index_date_otherneuroconditions, pre_index_date_dementia, pre_index_date_fh_premature_cvd, pre_index_date_primary_hhf, ckd3plus,
        year1oad_c
    )

### Relevel factors
med_ds1$ethn_5cat_analysis <- as.factor(med_ds1$ethn_5cat_analysis)
med_ds1$ethn_5cat_analysis <- relevel(med_ds1$ethn_5cat_analysis, ref = "White")

med_ds1$smk_qrisk2_acg <- as.factor(med_ds1$smk_qrisk2_acg)
med_ds1$smk_qrisk2_acg <- relevel(med_ds1$smk_qrisk2_acg, ref = "Non-smoker")

med_ds1$alcohol_cat <- as.factor(med_ds1$alcohol_cat)
med_ds1$alcohol_cat <- relevel(med_ds1$alcohol_cat, ref = "Within limits")

med_ds1$sex <- as.factor(med_ds1$sex)
med_ds1$sex <- relevel(med_ds1$sex, ref = "Male")

############################################################################################
### Section 14: Store dataset in secure approved location
############################################################################################

saveRDS(med_ds1, file="path_to_your_secure_area/analysis_dsyyyymmdd.rds")