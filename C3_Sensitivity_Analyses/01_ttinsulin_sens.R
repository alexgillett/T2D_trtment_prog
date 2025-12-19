############################################################################################
############################################################################################
### Sensitivity analyses: Time to insulin
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up (load R libraries, read in analysis dataset)
########################################
### Section 2: Complete case analysis
########################################
### Section 3: Complication count model (complete case)
########################################
### Section 4: Three-month landmark analysis
########################################
### Section 5: Check proportional odds assumption (time-varying depression effects)
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################
library(survival)
library(survminer)
library(rms)
library(dplyr)
library(flexsurv)

med_ds1 <-readRDS(file="path_to_your_secure_area/analysis_dsyyyymmdd.rds")

### ATTENTION: Set your working directory to where you want/ can securely store MI datasets

############################################################################################
### Section 2: Complete case analysis
############################################################################################
### Reduce to complete case
med_ds1 <- med_ds1 %>%
  mutate(
    cc_rm_all = ifelse(is.na(prehba1c) | is.na(prebmi) | is.na(e2019_imd_decile) | is.na(qrisk2_smoking_cat) | is.na(alcohol_cat) | is.na(ethnicity_5cat), 1, 0)
  )

med_cc <- med_ds1 %>% filter(cc_rm_all == 0)

### Reduce to columns needed for analysis
med_cc <- med_cc %>% 
    select(
        ins_ever, time_insulin_yrs, intens_event, time_change_yrs, lastpre1oaddep_quant, 
        age1oad, sex, disease_dur1oad_yrs, ethn_5cat_analysis, prebmi, prehba1c, 
        smk_qrisk2_acg, alcohol_cat, e2019_imd_decilef, gpcount1yr_excdep, pre_index_date_retinopathy, pre_index_date_diabeticnephropathy, 
        ckd3plus, pre_index_date_neuropathy, pre_index_date_hypertension, pre_index_date_ihd, pre_index_date_heartfailure, 
        pre_index_date_myocardialinfarction, pre_index_date_pad, pre_index_date_stroke, pre_index_date_tia, pre_index_date_cld, 
        pre_index_date_copd, pre_index_date_bronchiectasis, pre_index_date_asthma, pre_index_date_pulmonaryfibrosis, pre_index_date_pulmonaryhypertension, 
        pre_index_date_solid_cancer, pre_index_date_haem_cancer, pre_index_date_solidorgantransplant, pre_index_date_rheumatoidarthritis, 
        pre_index_date_otherneuroconditions, pre_index_date_dementia, pre_index_date_fh_premature_cvd, pre_index_date_primary_hhf,
        n_diabcomps, year1oad_c
    )

### Run Royston-Parmar PO model (as in primary analysis but with complete case data) 
mod1_po_cc <- flexsurvspline(Surv(time_insulin_yrs, ins_ever) ~ lastpre1oaddep_quant + 
  age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep +
  pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + ckd3plus + pre_index_date_neuropathy +
  pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction + 
  pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia +
  pre_index_date_cld + pre_index_date_primary_hhf + year1oad_c, 
  data = med_cc, k = 5, scale="odds")

### Save model
saveRDS(mod1_po_cc, file = "ttinsulin_cc_po.rds")

############################################################################################
### Section 3: Complication count model (complete case)
############################################################################################

mod2_po_cc <- flexsurvspline(Surv(time_insulin_yrs, ins_ever) ~ lastpre1oaddep_quant + 
  age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep + n_diabcomps + 
  pre_index_date_cld + year1oad_c, 
  data = med_cc, k = 5, scale="odds") 

AIC(mod1_po_cc, mod2_po_cc)

### Save model
saveRDS(mod2_po_cc, file = "ttinsulin_cc_po_count.rds")

############################################################################################
### Section 4: Three-month landmark analysis
############################################################################################
### remove those who change (add or switch), or are censored, within the first 3 months
med_cc_3m <- med_cc %>%
  filter(time_insulin_yrs > (3/12))

### Update time
med_cc_3m <- med_cc_3m %>%
  mutate(time_insulin_yrs = (time_insulin_yrs - (3/12)))

### Run model
mod3_po_landmark <- flexsurvspline(Surv(time_insulin_yrs, ins_ever) ~ lastpre1oaddep_quant + 
   age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep + 
  pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + ckd3plus + pre_index_date_neuropathy +
  pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction + 
  pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia +
  pre_index_date_cld + pre_index_date_primary_hhf + year1oad_c,
  data = med_cc_3m, k = 5, scale="odds") 

### Save model
saveRDS(mod3_po_landmark, file = "ttinsulin_cc_landmark.rds")

############################################################################################
### Section 5: Check proportional odds assumption (time-varying depression effects)
############################################################################################

### Time-varying depression effect 1: interaction with first log-time spline
mod4_cc_tvc1 <- flexsurvspline(Surv(time_insulin_yrs, ins_ever) ~ lastpre1oaddep_quant + 
   age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep + #n_diabcomps + 
  pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + ckd3plus + pre_index_date_neuropathy +
  pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction + 
  pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia +
  pre_index_date_cld + pre_index_date_primary_hhf + year1oad_c + 
  gamma1(lastpre1oaddep_quant), ### This adds an interaction between depression subgroups and the first spline of log-time
  data = med_cc, k = 5, scale="odds") 
end_timepotvc = Sys.time()

### Save model
saveRDS(mod4_cc_tvc1, file = "ttinsulin_tvc_check1.rds")

### Compare with no tvc model
AIC(mod4_cc_tvc1, mod1_po_cc) 

mod4_cc_tvc1$BIC 
mod1_po_cc$BIC 

### Time-varying depression effect 2: interaction with first and second log-time spline
mod5_cc_tvc2 <- flexsurvspline(Surv(time_insulin_yrs, ins_ever)~ lastpre1oaddep_quant + 
   age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep + 
  pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + ckd3plus + pre_index_date_neuropathy +
  pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction + 
  pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia +
  pre_index_date_cld + pre_index_date_primary_hhf + year1oad_c + 
  gamma1(lastpre1oaddep_quant) + gamma2(lastpre1oaddep_quant),
  data = med_cc, k = 5, scale="odds") 

### Save model
saveRDS(mod5_cc_tvc2, file = "ttinsulin_tvc_check2.rds")

### BIC
mod5_cc_tvc2$BIC 
mod4_cc_tvc1$BIC
mod1_po_cc$BIC 