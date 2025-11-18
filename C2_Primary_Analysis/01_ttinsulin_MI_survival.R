############################################################################################
############################################################################################
### Multiple imputation script: time to insulin initiation
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up (load R libraries, read in analysis dataset)
########################################
### Section 2: Set up to use Nelson-Aalen estimators
########################################
### Section 3: Perform Multiple Imputation and Fit Royston–Parmar Models (Loop Over 30 Imputations)
########################################
### Section 4: Pool imputed models and extract pooled estimates (PO & PH) for insulin outcome
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################
library(survival)
library(survminer)
library(dplyr)
library(flexsurv)
library(mice)
library(mitools)

mi_data <-readRDS(file="path_to_your_secure_area/analysis_dsyyyymmdd.rds")

### Remove other outcome
mi_data <- mi_data %>% select(-intens_event, -time_change_yrs)

### ATTENTION: Set your working directory to where you want/ can securely store MI datasets

############################################################################################
### Section 2: Set up to use Nelson-Aalen estimators
############################################################################################
mi_data$H0 <- nelsonaalen(mi_data, timevar = time_insulin_yrs, statusvar = ins_ever)
# setup run
imp0 <- mice(mi_data, maxit = 0)
pred <- imp0$predictorMatrix
# remove event time from predictor (high correlation with H0)
pred[, "time_insulin_yrs"] <- 0

############################################################################################
### Section 3: Perform multiple imputation and fit Royston–Parmar models (loop over 30 imputations)
############################################################################################
# This loop performs 30 independent imputations (m=1 each), saves each completed dataset, and fits 
# both proportional-odds and proportional-hazards Royston–Parmar models for the insulin initiation outcome.

imputation_range <- 1:30

# START LOOP
for (i in imputation_range) {
  message("Running imputation ", i)
  
  # Run 1 imputation with a unique seed
  imp_i <- mice(mi_data, m = 1, maxit = 20, seed = 1000 + i)
  
  # Extract the completed dataset
  data_i <- complete(imp_i, 1)
  
  # Save imputed dataset (optional)
  saveRDS(data_i, paste0("mi_na_insulin_", i, ".rds"))
  
  # Fit Royston-Parmar model
  rp_model_i <- flexsurvspline(
    Surv(time_insulin_yrs, ins_ever) ~ lastpre1oaddep_quant + 
  age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep + #n_diabcomps + 
  pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + ckd3plus + pre_index_date_neuropathy +
  pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction + 
  pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia +
  pre_index_date_cld + pre_index_date_primary_hhf + year1oad_c,
    data = data_i,
    k = 5,
    scale = "odds"
  )
  
  # Save model
  saveRDS(rp_model_i, paste0("mi_na_YEARrpm_insulin_", i, ".rds"))
  
  # Fit Royston-Parmar model
  rpph_model_i <- flexsurvspline(
    Surv(time_insulin_yrs, ins_ever) ~ lastpre1oaddep_quant + 
  age1oad + sex + disease_dur1oad_yrs + ethn_5cat_analysis + 
  prebmi + prehba1c + smk_qrisk2_acg + alcohol_cat + e2019_imd_decilef + gpcount1yr_excdep + #n_diabcomps + 
  pre_index_date_retinopathy + pre_index_date_diabeticnephropathy + ckd3plus + pre_index_date_neuropathy +
  pre_index_date_hypertension + pre_index_date_ihd + pre_index_date_heartfailure + pre_index_date_myocardialinfarction + 
  pre_index_date_pad + pre_index_date_stroke + pre_index_date_tia +
  pre_index_date_cld + pre_index_date_primary_hhf + year1oad_c,
    data = data_i,
    k = 5,
    scale = "hazard"
  )

  # Save model
  saveRDS(rpph_model_i, paste0("mi_na_YEARrpm_PH_insulin_", i, ".rds"))

  # Clean up
  rm(imp_i, data_i, rp_model_i, rpph_model_i)
  gc()
}

############################################################################################
### Section 4: Pool imputed models and extract pooled estimates (PO & PH) for insulin outcome
############################################################################################

### p-value function (taken from Stats Geek: https://thestatsgeek.com/2020/11/05/p-values-after-multiple-imputation-using-mitools-in-r/)
MIcombineP <- function(MIcombineRes) {
  tStat <- MIcombineRes$coefficients/sqrt(diag(MIcombineRes$variance))
  2*pt(-abs(tStat),df=MIcombineRes$df)
}

############################################
### Read in PO models and extract estimates
############################################
n_imputations <- 30
rp_models <- list()
coef_list <- list()
vcov_list <- list()
bic_po_list <- list()
for (i in 1:n_imputations) {
  model_i <- readRDS(paste0("mi_na_YEARrpm_insulin_", i, ".rds"))
  rp_models[[i]] <- model_i
  coef_list[[i]] <- coef(model_i)
  vcov_list[[i]] <- vcov(model_i)
  bic_po_list[[i]] <- BIC(model_i)
}

############################################
### Pool PO models using Rubin's Rules
############################################
pooled <- MIcombine(results = coef_list, variances = vcov_list)
summary_pooled <- summary(pooled)
pooled_p <- MIcombineP(pooled)

### Add p-values to to summary
summary_with_p <- cbind(summary_pooled, p_value = pooled_p)
summary_with_p

### Save the final results
saveRDS(pooled, "na_rp_insulin_YEARpooled.rds")
saveRDS(summary_pooled, "na_pooled_insulin_YEARsummary.rds")

############################################
### Extract and exponentiate depression terms (ORs)
############################################

# Rows corresponding to depression recency levels
depression_rows <- summary_pooled[grep("lastpre1oaddep_quantMDD", rownames(summary_pooled)), ]

# Compute pooled ORs and 95% CIs
depression_ORs <- data.frame(
  Variable = rownames(depression_rows),
  OR       = exp(depression_rows$results),
  LowerCI  = exp(depression_rows$results - 1.96 * depression_rows$se),
  UpperCI  = exp(depression_rows$results + 1.96 * depression_rows$se)
)

# Print and save ORs
print(depression_ORs)
saveRDS(depression_ORs, "na_pooled_depOR_insulin_YEARmi.rds")

############################################
### Read in PH models and extract estimates
############################################
rph_models   <- list()
coef_list    <- list()
vcov_list    <- list()
bic_ph_list  <- list()
for (i in 1:n_imputations) {
  # Read Royston–Parmar proportional hazards model for imputation i
  model_i <- readRDS(paste0("mi_na_YEARrpm_PH_insulin_", i, ".rds"))
  
  rph_models[[i]]  <- model_i
  coef_list[[i]]   <- coef(model_i)
  vcov_list[[i]]   <- vcov(model_i)
  bic_ph_list[[i]] <- BIC(model_i)
}

############################################
### Pool PH models using Rubin's Rules
############################################
pooled <- MIcombine(results = coef_list, variances = vcov_list)
summary_pooled <- summary(pooled)

### Save the final results
saveRDS(pooled, "na_rpPH_insulin_YEARpooled.rds")
saveRDS(summary_pooled, "na_PH_pooled_insulin_YEARsummary.rds")

############################################
### Extract and exponentiate depression terms (HRs)
############################################

depression_rows <- summary_pooled[grep("lastpre1oaddep_quantMDD", rownames(summary_pooled)), ]

depression_PHs <- data.frame(
  Variable = rownames(depression_rows),
  HR       = exp(depression_rows$results),
  LowerCI  = exp(depression_rows$results - 1.96 * depression_rows$se),
  UpperCI  = exp(depression_rows$results + 1.96 * depression_rows$se)
)

print(depression_PHs)
saveRDS(depression_PHs, "na_pooled_depPH_insulin_YEARmi.rds")

############################################
### Export BIC values for PO vs PH comparison
############################################
bic_ph <- unlist(bic_ph_list)
bic_po <- unlist(bic_po_list)

out <- cbind(bic_po, bic_ph)
write.csv(out, file = "insulinbic_po_ph.csv")