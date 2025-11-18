############################################################################################
############################################################################################
### Define the index dates (first oral GLM monotherapy prescription) AND
### Define (initial) outcome variables (time to intensification, time to insulin)
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
# Load R libraries, link to your CPRD data via the aurum package, load in existing tables needed in script
########################################
### Section 2: Create tables for longitudinal medication episode information for patients
# Create longitudinal drug- and chemical-class–level treatment histories for each patient. 
# It standardises substance names, assigns chemical classes, and derives continuous treatment episodes (t2dmeds_drgepisode, t2dmeds_chemepisode, t2dmeds_chemepisode30days) with start and end dates.
########################################
### Section 3: Define index_dates
########################################
### Section 4: Define initial outcomes
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

analysis = cprd$analysis(cohort_prefix)
# Assign name and location to previously extracted raw T2D prescription data
raw_oha    <- raw_oha_prodcodes %>% analysis$cached("raw_oha_prodcodes")
raw_insulin<- raw_insulin_prodcodes %>% analysis$cached("raw_insulin_prodcodes")

# Note 'ids2004apr_pracex' is a MySql table containing IDs (patid) of patients with T2D excluding those with a T2D diagnosis before 
# 01/04/2004 and excluding those from a merged GP practice (to avoid double counting)
analysis = cprd$analysis(cohort_prefix)
ids2004apr_pracex <- ids2004apr_pracex %>% analysis$cached("ids2004apr_pracex")

# Note- created in script 01_create_death_dat
analysis = cprd$analysis("all_patid")
death_end_dat <- death_end_dat %>% analysis$cached("death_end_dat")

# Clean prescription data:
# Join prescription tables with death_dat_end and remove events occuring before DOB and after end date. 
# Join with the IDs initially be be included (T2D patients diagnosed after 1st April 2004, within valid practices).
# Join with the product dictionary,
analysis = cprd$analysis("cohort_prefix")
clean_oha_proddict <- raw_oha %>% 
    inner_join(death_end_dat, by="patid") %>%
    filter(issuedate>=min_dob & issuedate<=gp_death_end_date) %>%
    select(-min_dob, -gp_death_end_date, -gp_end_date, -cprd_ddate) %>%
    inner_join(ids2004apr_pracex, by="patid") %>%
    inner_join(cprd$tables$productDict, by="prodcodeid") %>%
    analysis$cached("clean_oha_proddict", indexes=c("patid", "issuedate", "prodcodeid", "dosageid"))

clean_insulin_proddict <- raw_insulin %>% 
    inner_join(death_end_dat, by="patid") %>%
    filter(issuedate>=min_dob & issuedate<=gp_death_end_date) %>%
    select(-min_dob, -gp_death_end_date, -gp_end_date, -cprd_ddate) %>%
    inner_join(ids2004apr_pracex, by="patid") %>%
    inner_join(cprd$tables$productDict, by="prodcodeid") %>%
    analysis$cached("clean_insulin_proddict", indexes=c("patid", "issuedate", "prodcodeid", "dosageid"))

analysis = cprd$analysis("cohort_prefix")
clean_oha_proddict <- clean_oha_proddict %>% analysis$cached("clean_oha_proddict")
clean_insulin_proddict <- clean_insulin_proddict %>% analysis$cached("clean_insulin_proddict")

############################################################################################
### Section 2: Create tables for longitudinal medication episode information for patients
############################################################################################

########################################
### Handle missing drugsubstancename for OHA
### (Lyxumia pack -> Lixisenatide, Subcutaneous)
########################################

analysis <- cprd$analysis("alex_t2")

drugsub_oha_tmp <- clean_oha_proddict %>% 
  select(drugsubstancename, routeofadministration) %>%
  distinct() %>% 
  mutate(
    drugsubstancename2     = ifelse(is.na(drugsubstancename), "Lixisenatide", drugsubstancename),
    routeofadministration2 = ifelse(is.na(drugsubstancename), "Subcutaneous", routeofadministration)
  ) %>% 
  select(drugsubstancename2, routeofadministration2) %>%
  rename(
    drugsubstancename     = drugsubstancename2,
    routeofadministration = routeofadministration2
  ) %>%
  distinct()

########################################
### Map substance to chemical class (OHA)
### -> oha_chemname_ref  (chemname1, chemname2)
########################################

analysis = cprd$analysis("cohort_prefix")
oha_chemname_ref <- drugsub_oha_tmp %>% 
    mutate(chemname1 = case_when(
      drugsubstancename %in% c("Gliclazide", "Tolbutamide", "Glimepiride", "Glibenclamide", "Glipizide", "Gliquidone", "Chlorpropamide") ~ "Sulfonylureas",
      drugsubstancename %in% c("Metformin hydrochloride", "Metformin hydrochloride/ Vildagliptin", "Metformin hydrochloride/ Rosiglitazone maleate", 
      "Metformin hydrochloride/ Pioglitazone hydrochloride", "Metformin hydrochloride/ Saxagliptin hydrochloride", 
      "Linagliptin/ Metformin hydrochloride", "Metformin hydrochloride/ Sitagliptin", "Empagliflozin/ Metformin hydrochloride",
      "Alogliptin benzoate/ Metformin hydrochloride", "Canagliflozin hemihydrate/ Metformin hydrochloride", 
      "Dapagliflozin propanediol monohydrate/ Metformin hydrochloride") ~ "Metformin",
      drugsubstancename %in% c("Pioglitazone hydrochloride", "Rosiglitazone maleate") ~ "Thiazolidinedione",
      drugsubstancename %in% c("Sitagliptin", "Sitagliptin hydrochloride", "Saxagliptin hydrochloride", "Alogliptin benzoate", "Linagliptin", 
      "Vildagliptin", "Empagliflozin/ Linagliptin", "Dapagliflozin propanediol monohydrate/ Saxagliptin hydrochloride") ~ "DPP4",
      drugsubstancename %in% c("Exenatide", "Semaglutide", "Liraglutide", "Lixisenatide", "Dulaglutide", 
      "Insulin degludec/ Liraglutide", "Insulin glargine/ Lixisenatide", "Albiglutide") ~ "GLP1",
      drugsubstancename %in% c("Nateglinide", "Repaglinide") ~ "Meglitinides",
      drugsubstancename %in% c("Dapagliflozin propanediol monohydrate", "Empagliflozin", "Canagliflozin hemihydrate", 
      "Ertugliflozin L-pyroglutamic acid") ~ "SGLT2",
      drugsubstancename %in% c("Acarbose") ~ "Acarbose",
      drugsubstancename %in% c("Tirzepatide") ~ "GIP_GLP1",
      TRUE ~ "Unknown" # Fallback for unmatched cases
    ),
    chemname2 = case_when(
      drugsubstancename %in% c("Metformin hydrochloride/ Vildagliptin", "Metformin hydrochloride/ Saxagliptin hydrochloride", "Linagliptin/ Metformin hydrochloride", "Metformin hydrochloride/ Sitagliptin", "Alogliptin benzoate/ Metformin hydrochloride") ~ "DPP4",
      drugsubstancename %in% c("Metformin hydrochloride/ Rosiglitazone maleate", "Metformin hydrochloride/ Pioglitazone hydrochloride") ~ "Thiazolidinedione",
      drugsubstancename %in% c("Insulin degludec/ Liraglutide", "Insulin glargine/ Lixisenatide") ~ "Insulin",
      drugsubstancename %in% c("Empagliflozin/ Metformin hydrochloride", "Canagliflozin hemihydrate/ Metformin hydrochloride", 
      "Empagliflozin/ Linagliptin", "Dapagliflozin propanediol monohydrate/ Metformin hydrochloride", 
      "Dapagliflozin propanediol monohydrate/ Saxagliptin hydrochloride") ~ "SGLT2",
      TRUE ~ "None" # Fallback for unmatched cases
    )) %>% 
    analysis$cached("oha_chemname_ref")

########################################
### Create clean_oha_chemname (join OHA with chem class mapping)
########################################

clean_oha_chemname <- clean_oha_proddict %>% 
    mutate(
        drugsubstancename2 = ifelse(is.na(drugsubstancename), "Lixisenatide", drugsubstancename),
        routeofadministration=ifelse(is.na(drugsubstancename), "Subcutaneous", routeofadministration)
        ) %>%
    select(-drugsubstancename) %>%
    rename(drugsubstancename=drugsubstancename2) %>%
    inner_join(oha_chemname_ref, by=c("drugsubstancename", "routeofadministration")) %>%
    analysis$cached("clean_oha_chemname", indexes=c("patid", "issuedate", "prodcodeid", "dosageid"))

########################################
### Create clean_insulin_chemname (simple mapping)
########################################

clean_insulin_chemname <- clean_insulin_proddict %>% 
  mutate(
    chemname1 = case_when(
      drugsubstancename %in% c(
        "Insulin degludec/ Liraglutide",
        "Insulin glargine/ Lixisenatide "
      ) ~ "GLP1",
      TRUE ~ "Insulin"
    )
  ) %>%
  mutate(
    chemname2 = ifelse(chemname1 != "Insulin", "Insulin", "None")
  ) %>%
  analysis$cached("clean_insulin_chemname",
                  indexes = c("patid", "issuedate", "prodcodeid", "dosageid"))

########################################
### Combine OHA + insulin for downstream prescription episode creation
########################################

clean_oha_chemname1 <- clean_oha_chemname %>% 
  select(patid, issuedate, drugsubstancename, chemname1, chemname2)
clean_insulin_chemname1 <- clean_insulin_chemname %>% 
  select(patid, issuedate, drugsubstancename, chemname1, chemname2)

t2dmeds_combi <- clean_oha_chemname1 %>% 
  union_all(clean_insulin_chemname1) %>%
  distinct()

######################################################################################################
### 2A. Drug-level prescription episodes: t2dmeds_drgepisode
######################################################################################################

t2dmeds_drgepisode <- t2dmeds_combi %>%
    #filter(row_number() <= 5000) %>%  # Limited to 5000 rows during testing phase
    mutate(
        prev_drug_dt = sql("LAG(issuedate) OVER (PARTITION BY patid, drugsubstancename ORDER BY issuedate)"),
        diff_days_drug = sql("DATEDIFF(issuedate, LAG(issuedate) OVER (PARTITION BY patid, drugsubstancename ORDER BY issuedate))"),
        new_episode_flag = case_when(
            is.na(diff_days_drug) | diff_days_drug > 183 ~ 1,
            TRUE ~ 0
        )
    ) %>%
    mutate(
        prescription_episode = sql("SUM(new_episode_flag) OVER (PARTITION BY patid, drugsubstancename ORDER BY issuedate)")
    ) %>%
    group_by(patid, drugsubstancename, prescription_episode) %>%
    mutate(
        drug_start = sql("MIN(issuedate) OVER (PARTITION BY patid, drugsubstancename, prescription_episode)"),
        drug_end = sql("MAX(issuedate) OVER (PARTITION BY patid, drugsubstancename, prescription_episode)")
    ) %>%
    ungroup() %>%
    select(patid, drugsubstancename, chemname1, chemname2, prescription_episode, drug_start, drug_end) %>%
    distinct() %>%
    arrange(patid, drug_start) %>%
    analysis$cached("t2dmeds_drgepisode", indexes=c("patid", "drug_start", "drug_end"))

######################################################################################################
### 2B. Chemical-level prescription episodes: t2dmeds_chemepisode
######################################################################################################

t2dmeds_chemepisode <- t2dmeds_combi %>%
    mutate(chemname = chemname1) %>%
    select(patid, chemname, issuedate) %>%
    union_all(
        t2dmeds_combi %>%
            filter(chemname2 != "None") %>%
            mutate(chemname = chemname2) %>%
            select(patid, chemname, issuedate)
    ) %>%
    distinct() %>%
    mutate(
        prev_chem_dt = sql("LAG(issuedate) OVER (PARTITION BY patid, chemname ORDER BY issuedate)"),
        diff_days_chem = sql("DATEDIFF(issuedate, LAG(issuedate) OVER (PARTITION BY patid, chemname ORDER BY issuedate))")
    ) %>%
    mutate(
        new_chem_episode_flag = case_when(
            is.na(diff_days_chem) | diff_days_chem > 183 ~ 1,
            TRUE ~ 0
        )
    ) %>%
    mutate(
        chem_prescription_episode = sql("
            SUM(new_chem_episode_flag) OVER (PARTITION BY patid, chemname ORDER BY issuedate)
        ")
    ) %>%
    group_by(patid, chemname, chem_prescription_episode) %>%
    mutate(
        chem_start = sql("MIN(issuedate) OVER (PARTITION BY patid, chemname, chem_prescription_episode)"),
        chem_end = sql("MAX(issuedate) OVER (PARTITION BY patid, chemname, chem_prescription_episode)")
    ) %>%
    ungroup() %>%
    select(patid, chemname, chem_prescription_episode, chem_start, chem_end) %>%
    distinct() %>%
    arrange(patid, chem_start) %>%
    analysis$cached("t2dmeds_chemepisode", indexes = c("patid", "chem_start", "chem_end"))

######################################################################################################
### 2C. Chemical-level episodes with padded end date: t2dmeds_chemepisode30days
######################################################################################################

t2dmeds_chemepisode30days <- t2dmeds_chemepisode %>%
   mutate(
        chem_end2 = sql("DATE_ADD(chem_end, INTERVAL 30 DAY)")  # Adds 30 days to chem_end
    ) %>%
    select(-chem_end) %>%
    rename(chem_end = chem_end2) %>%
    arrange(patid, chem_start) %>%
    analysis$cached("t2dmeds_chemepisode30days", indexes = c("patid", "chem_start", "chem_end"))

######################################################################################################
### Section 3: Define index_dates
######################################################################################################

##########################################################################
### 3A. From longitudinal chemical table identify the first prescription
##########################################################################
# At this stage, this could include injections
analysis = cprd$analysis("analysis_prefix")
first_oad_prescription <- t2dmeds_chemepisode30days %>% 
    mutate(first_dt = sql("MIN(chem_start) OVER (PARTITION BY patid)")) %>%
    mutate(is_min = as.numeric(chem_start == first_dt)) %>%
    filter(is_min == 1) %>%
    mutate(
        num_chems = sql("COUNT(*) OVER (PARTITION BY patid)")
    ) %>%
    filter(num_chems == 1) %>%
    select(-is_min, -first_dt) %>%
    analysis$cached("first1med_dts_interim", indexes=c("patid", "chemname", "chem_start"))

# Look at the distinct chemical names 
first_oad_prescription %>%
    distinct(chemname)
#1 Acarbose         
#2 DPP4             
#3 GLP1             
#4 Insulin          
#5 Meglitinides     
#6 Metformin        
#7 SGLT2            
#8 Sulfonylureas    
#9 Thiazolidinedione

# Want oral OADs only
# Note. For GLP1 only Rybelsus is oral and whilst in the prodcode list used none appear to be this drug so...
#glp1_check <- t2dmeds_drgepisode %>% 
#    filter(chemname1 == "GLP1" | chemname2 == "GLP1") %>%
#    summarise(out = distinct(drugsubstancename))

########################################
### To remove: Insulin, GLP1
########################################
analysis = cprd$analysis("analysis_prefix")
first_oad_prescription <- first_oad_prescription %>% 
    filter(chemname != "Insulin") %>%
    filter(chemname != "GLP1") %>%
    analysis$cached("first1OAD_date", indexes=c("patid", "chemname", "chem_start"))

########################################
### Store IDs with first prescription being a single oral antidiabetic medication
########################################
ids1OAD <- first_oad_prescription %>% select(patid) %>% distinct()

##########################################################################
### 3B. Store index dates- just for clarity/ align with Exeter code
##########################################################################
analysis = cprd$analysis("analysis_prefix")
index_dates <- first_oad_prescription %>%
  select(patid, index_date = chem_start) %>%
  distinct()%>%
    analysis$cached("index_dates1oad", indexes=c("patid", "index_date"))
  

######################################################################################################
### Section 4: Define outcomes
######################################################################################################

##########################################################################
### 4A. The first insulin date (where applicable)
##########################################################################
analysis = cprd$analysis("alex_t2")
insulin_firstdate <- t2dmeds_drgepisode %>% 
    mutate(
        ins_ever = ifelse(
            chemname1 == "Insulin", 
            1,
            ifelse(
                chemname2 == "Insulin",
                1,
                0
                )
            )
        ) %>%
    filter(ins_ever == 1) %>%
    mutate(first_ins_dt = sql("MIN(drug_start) OVER (PARTITION BY patid)")) %>%
    filter(drug_start == first_ins_dt) %>%
    rename(ins1st_drugtype = drugsubstancename, ins1stchem1=chemname1, ins1stchem2=chemname2, ins1st_date = first_ins_dt) %>%
    select(-prescription_episode, -drug_start, -drug_end) %>%
    distinct() %>%
    analysis$cached("insulin_firstdate", indexes=c("patid", "ins1st_date"))
### Note, some people have multiple insulin entries on same day.
#insulin_firstdate %>% summarise(out = distinct(ins1stchem2))

### IDs and date of first insulin prescriptions (for merging later)
insulin_tojoin <- insulin_firstdate %>% 
    select(patid, ins_ever, ins1st_date) %>%
    distinct()

##########################################################################
### 4B. Define intensification 
##########################################################################

########################################
### Add 6 month buffer to episode length and define the first OAD treatment
########################################

# We restrict to patients starting on one OAD, add a 6-month buffer to the end of the first chemical episode, and 
# identify the first treatment (chemname) and its start/end dates (which define the index date and the window for 
# detecting treatment changes).
meds_buffered <- t2dmeds_chemepisode %>% 
    inner_join(ids1OAD) %>%
    mutate(
    chem_end_183 = sql("DATE_ADD(chem_end, INTERVAL 183 DAY)")
  )  

# What is the first drug someone gets...
first_drug <- meds_buffered %>%
  mutate(
    rownum = sql("ROW_NUMBER() OVER (PARTITION BY patid ORDER BY chem_start)")
  ) %>%
  filter(rownum == 1) %>%
  select(patid, first_chem = chemname, first_start = chem_start, first_end = chem_end_183)

########################################
### Identify subsequent medications and their overlap with the first OAD
########################################

# Get all subsequent meds...
# Capture all chemical episodes that occur after the first OAD and quantify their timing relative to the buffered first episode.
For each subsequent medication, we calculate how many days it overlaps with the first OAD course and flag whether it overlaps, occurs within 6 months, or starts after the buffer.
subsequent_meds <- meds_buffered %>%
  inner_join(first_drug, by = "patid") %>%
  filter(chem_start > first_start) %>%
  mutate(
    overlap_days = sql("DATEDIFF(first_end, chem_start)")
  ) %>% mutate( 
    is_overlap = overlap_days >= 0,
    is_within_6mo = overlap_days >= 0 & overlap_days <= 183,
    is_after_buffer = chem_start > first_end
  )

########################################
### Classify treatment changes among patients with overlapping meds
########################################
# Among patients with overlapping second-line medications (with the 6 month buffer), classify the type of treatment change.
# For each patient we summarise the number of new chemical classes and maximum overlap, and use this to label 
# the change as add-on, switch, or multi-drug variants, with date_change capturing the earliest intensification date.

### Classify candidate chemicals- keep chemicals with any overlap with first chemical
# Note: switch and add on with multi here doesn't mean that an individual switched to multi or immediately added multi, 
# just that overlapping with prescription 1, multiple drugs were tried either together or sequentially
analysis = cprd$analysis("tt1oad")
candidate_meds <- subsequent_meds %>%
  filter(is_overlap==1) %>%
  group_by(patid) %>%
  summarise(
    n_new_drugs = n_distinct(chemname),
    max_overlap_days = max(overlap_days),  # Maximum overlap — for classification
    date_change = min(chem_start),
    .groups = "drop" # ungroups
  ) %>%
  mutate(
    change_type = case_when(
      n_new_drugs == 1 & max_overlap_days >= 183 ~ "add on",
      n_new_drugs == 1 & max_overlap_days < 183 ~ "switch",
      n_new_drugs > 1 & max_overlap_days >= 183 ~ "add on multi",
      n_new_drugs > 1 & max_overlap_days < 183 ~ "switch multi"
    )
  ) %>% analysis$cached("change_switch_interim", indexes=c("patid"))

########################################
### Define “drop” and “no change” patterns and combine all change types
########################################

# Classify patients who do not meet the overlap-based intensification criteria as either dropping treatment or 
# remaining on the initial OAD only.
# We label patients with subsequent non-overlapping meds as “drop” and those with no further meds as “no change”, 
# then combine these with the add-on/switch categories into a single change dataset per patient.

### initial dropped medication IDs
candidate_meds_patid <- candidate_meds %>% select(patid) %>% distinct()

analysis = cprd$analysis("tt1oad")
drop_dat <- subsequent_meds %>%
    anti_join(candidate_meds_patid, by = "patid") %>%
    mutate(change_type = "drop",
    date_change = first_end) %>%
    select(patid, change_type, date_change) %>% distinct() %>%
    analysis$cached("chem_drop_interim", indexes=c("patid"))


### initial no change dataset (some will be drops)
patients_sub_drugs <- subsequent_meds %>% select(patid) %>% distinct()

analysis = cprd$analysis("tt1oad")
no_change_interim <- first_drug %>%
  anti_join(patients_sub_drugs, by = "patid") %>%
  mutate(change_type = "no change",
  date_change = first_end) %>%
  select(patid, change_type, date_change) %>% distinct() %>%
  analysis$cached("nochange_interim", indexes=c("patid"))

medication_changes_interim <- candidate_meds %>% select(patid, change_type, date_change) %>%
  union_all(drop_dat) %>%
  union_all(no_change_interim)

########################################
### Incorporate death timing and derive the intensification event indicator
########################################

# Align treatment change definitions with death and derive the intensification event indicator.
# We check whether death occurs before the planned change date, adjust date_change to death where appropriate, 
# and define intens_event (1 = any add-on/switch, 0 = no change/drop) with date_change acting as either the event date or censoring time.

medication_changes_interim2 <- medication_changes_interim %>% 
    inner_join(death_end_dat %>% select(patid, gp_death_end_date)) %>%
    mutate(check = sql("DATEDIFF(gp_death_end_date, date_change)")) %>%
    mutate(death_first = ifelse(check < 0, 1, 0)) ### less than 0 death before change

medication_changes <- medication_changes_interim2 %>%
    mutate(
        change_type2 = ifelse(
            change_type != "no change",
            change_type,
            ifelse(
                death_first == 1,
                "no change",
                "drop"
            )
        )
    ) %>%
    mutate(date_change = ifelse(death_first == 1, gp_death_end_date, date_change)) %>%
    select(patid, change_type2, date_change, gp_death_end_date) %>%
    rename(change_type=change_type2) %>%
    mutate(intens_event = case_when(
        change_type %in% c("no change", "drop") ~ 0,
        !change_type %in% c("no change", "drop") ~ 1
    )) %>%
    analysis$cached("chem_intens", indexes=c("patid", "date_change", "gp_death_end_date"))

########################################
### Create patient-level time-to-intensification and time-to-insulin outcomes
########################################
# Construct a preliminary survival dataset with time-to-event outcomes.
# We join insulin initiation and intensification data to the first OAD episode, 
# define event/censor dates for both treatment intensification and insulin initiation, and 
# compute time from index (first_start) to each outcome in days and years, storing the result as ttprogression_v1.

analysis = cprd$analysis("tt1oad")
first_change <- first_drug %>% 
    mutate(
        first_chem_wks = sql("DATEDIFF(first_end, first_start) / 7")
    ) %>%
    left_join(insulin_tojoin) %>%
    left_join(medication_changes) %>%
    mutate(
        ins_ever = ifelse(is.na(ins_ever), 0, ins_ever)
    ) %>%
    mutate(date_ins_cens = ifelse(
        ins_ever == 1,
        ins1st_date,
        gp_death_end_date
    )) %>%
    mutate(
       tte_change_days = sql("DATEDIFF(date_change, first_start)"),
       tte_insulin_days = sql("DATEDIFF(date_ins_cens, first_start)")
    ) %>%
    mutate(
       time_change_yrs = tte_change_days/365.25,
       time_insulin_yrs = tte_insulin_days/365.25
    ) %>%
    analysis$cached("ttprogression_v1", indexes = c("patid", "first_start", "first_end", "ins1st_date"))
