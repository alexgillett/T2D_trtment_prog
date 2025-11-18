############################################################################################
############################################################################################
### Depression at index (y/n) and time between last depression code and index
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
# Load R libraries, link to your CPRD data via the aurum package, load in existing tables needed in script
########################################
### Section 2: Extract and clean raw depression data
########################################
### Section 3: Identify the closest depression code to index
########################################
### Section 4: Count the number of depression-code dates in the year prior to index
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################

### Before usual set up to work with aurum, we need to add a table to the analysis database.
# We did this using RMariaDB due to data being on a university server.
# Launch R and connect via RMariaDB
library(data.table)
library(DBI)
library(RMariaDB)
library(dbplyr)
library(dplyr)

# Create a connection... e.g.
con2 <- dbConnect(MariaDB(), 
                 dbname = "Update to your db name", 
                 host = "127.0.0.1", 
                 port = 3307,
                 user = "Update to your user name",
                 password = "Update to your password"
               )

# Included in the github is the file: 'https://github.com/alexgillett/T2D_trtment_prog/Depression_codes/depression_broad_medcode_terms.txt'
# We need to add this file
# Read the tab-delimited file into R
# E.g.
file_path <- "Update to your path/dep_broad_medcode_terms.txt"
data <- fread(file_path)

# Create the MySql structure for this table
dbExecute(con2, "
CREATE TABLE medcode_term (
    medcodeid BIGINT,
    term TEXT
);
")

# Insert data into the MySQL table
dbWriteTable(con2, name = "medcode_term", value = data, append = TRUE, row.names = FALSE)

# Verify data
med_terms <- dbReadTable(con2, "medcode_term")
print(head(med_terms))

# Close the connection
dbDisconnect(con2)

# Then work with aurum as usual...
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

analysis = cprd$analysis("all_patid")
death_end_dat <- death_end_dat %>% analysis$cached("death_end_dat")

analysis = cprd$analysis(analysis_prefix)
index_dates <- index_dates1oad %>% analysis$cached("index_dates1oad")

analysis = cprd$analysis("medcode")
medcode_term <- medcode_term %>% analysis$cached("term")

############################################################################################
### Section 2: Extract and clean broad depression
############################################################################################

# depression broad codes extraction
analysis = cprd$analysis(cohort_prefix)
raw_dep_broad_medcodes <- cprd$tables$observation %>%
  inner_join(codes$depression_broad, by="medcodeid") %>%
  analysis$cached("raw_dep_broad_medcodes", indexes=c("patid", "obsdate", "depression_broad_cat"))

# Take raw depression codes, remove FH codes and those occuring outside of birth-death dates...
dep_codes1 <- raw_dep_broad_medcodes %>%
  select(patid, obsdate, obstypeid, medcodeid) %>%
  filter(obstypeid!=4) %>% ### gets rid of FH codes...
  inner_join(death_end_dat, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_death_end_date) ### filters to be between birth and death...

# Create a list of depression codes to keep when identifying last code before index
# Want to exclude remission codes, resolved codes, review invite related codes, history of depression codes, 
# removed from depression register codes
codes_keep <- codes$depression_broad %>% 
    left_join(medcode_term)
codes_keep <- codes_keep %>% 
    filter(!depression_broad_cat %in% c("remission", "resolved", "review_invite")) %>%  
    filter(!term %in% c("H/O: depression", "History of depression", 
                        "History of postnatal depression", "Removed from depression register")) 

### For individuals in the raw depression codes table, restrict to those with: 
# at least one depression code in the restricted list and
# are individuals with an index date (so initiate T2D treatment with a single OAD)

### And we only want codes in the codes_keep list...
### This does mean that individuals with only remission, resolved, etc codes will be considered no history of depression as
### we cannot put a date on the last episode.
analysis = cprd$analysis(analysis_prefix)
dep_codes1 <- dep_codes1 %>% 
    inner_join(codes_keep) %>%
    inner_join(index_dates) %>%
    analysis$cached("tmp1oad_depression", indexes=c("patid", "index_date", "obsdate", "medcodeid", "gp_death_end_date", "cprd_ddate"))


############################################################################################
### Section 3: Identify the closest depression code to index
############################################################################################

pre1oad_date <- dep_codes1 %>% 
    filter(obsdate <= index_date) %>% ### pre-existing depression code dates...
    mutate(
        max_obsdate = sql("MAX(obsdate) OVER (PARTITION BY patid)")
    ) %>%
    filter(obsdate == max_obsdate) %>%
    distinct() %>%
    mutate(
        row_count = sql("COUNT(*) OVER (PARTITION BY patid)")  # Count rows per patid
    ) %>%
    analysis$cached("pre1oad_lastdate", indexes=c("patid", "obsdate", "index_date"))

### Some individuals will have multiple rows in this dataset (multiple depression codes on the same date)
### This will be shown in the row_count column
### For our dataset the maximum number per individual is 8
### If your version has a different max you will need to update the following code

### The following code takes the long table of depression codes and converts it into a wide, per-patient summary containing 
# the first 8 depression terms and categories, then keeps only one summary row per patient 
# I.e. converts to wide format
pre1oad_date <- pre1oad_date %>%
  mutate(
    row_num = sql("ROW_NUMBER() OVER (PARTITION BY patid ORDER BY obsdate)")
  ) %>%
  mutate(
    term1 = sql("MAX(CASE WHEN row_num = 1 THEN term END) OVER (PARTITION BY patid)"),
    term2 = sql("MAX(CASE WHEN row_num = 2 THEN term END) OVER (PARTITION BY patid)"),
    term3 = sql("MAX(CASE WHEN row_num = 3 THEN term END) OVER (PARTITION BY patid)"),
    term4 = sql("MAX(CASE WHEN row_num = 4 THEN term END) OVER (PARTITION BY patid)"),
    term5 = sql("MAX(CASE WHEN row_num = 5 THEN term END) OVER (PARTITION BY patid)"),
    term6 = sql("MAX(CASE WHEN row_num = 6 THEN term END) OVER (PARTITION BY patid)"),
    term7 = sql("MAX(CASE WHEN row_num = 7 THEN term END) OVER (PARTITION BY patid)"),
    term8 = sql("MAX(CASE WHEN row_num = 8 THEN term END) OVER (PARTITION BY patid)")
  ) %>%
  mutate(
    codecat1 = sql("MAX(CASE WHEN row_num = 1 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat2 = sql("MAX(CASE WHEN row_num = 2 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat3 = sql("MAX(CASE WHEN row_num = 3 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat4 = sql("MAX(CASE WHEN row_num = 4 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat5 = sql("MAX(CASE WHEN row_num = 5 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat6 = sql("MAX(CASE WHEN row_num = 6 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat7 = sql("MAX(CASE WHEN row_num = 7 THEN depression_broad_cat END) OVER (PARTITION BY patid)"),
    codecat8 = sql("MAX(CASE WHEN row_num = 8 THEN depression_broad_cat END) OVER (PARTITION BY patid)")
  ) %>%
    filter(row_num == 1) %>%
    analysis$cached("pre1oad_lastdate_v2", indexes=c("patid", "obsdate"))


pre1oaddep_date <- pre1oad_date %>%
  select(patid, index_date, max_obsdate, term1, term2, term3, term4, term5, term6, term7, term8,
    codecat1:codecat8) %>%
  rename(pre1oaddep_lastdate = max_obsdate) %>%
  mutate(lastpre1oaddep_time = sql("DATEDIFF(index_date, pre1oaddep_lastdate)")/365.25) %>%
  mutate(depr_yn = 1) %>%
  select(-index_date) %>% analysis$cached("pre1oad_lastdate_v3", indexes=c("patid", "pre1oaddep_lastdate"))

############################################################################################
### Section 4: Count the number of depression-code dates in the year prior to index
############################################################################################

  pret2ddep_code1y <- dep_codes1 %>% 
    mutate(
        dm_diag_minus1y = sql("DATE_SUB(dm_diag_date, INTERVAL 1 YEAR)")  # Create a date 1 year prior to T2D
    ) %>%
    filter(
        obsdate >= dm_diag_minus1y & obsdate <= dm_diag_date  # Filter depression codes within the 1-year window
    ) %>%
    select(
        patid, gender, dob, dm_diag_date, dm_diag_age, obsdate  
    ) %>%
    distinct() %>%  # Focus on unique dates
    mutate(
        n_dep_dates1y = sql("COUNT(*) OVER (PARTITION BY patid)")
    ) %>%
    select(patid, n_dep_dates1y) %>% 
        distinct() %>%
    analysis$cached(
        "pret2ddep_code1y",  # Cache the result with a new name
        indexes = c("patid")
    )
