-- Create a table 'all_patid_death_end_dat' within the project database (denoted 'prj' below) 
-- that is equivalent to r_valid_date_lookup (created during the CPRD-MySql set up) but with:
-- 1) cprd_ddate added in and 
-- 2) a column called gp_death_end_date which is the minimum of cprd_ddate and gp_end_date (end of follow-up from GP practice)

-- 'gp_death_end_date' is then used instead of the 'gp_ons_end_date' column from r_valid_date_lookup (validDateLookup in R code) in this study

-- Step 1: Copy the original look-up table table
CREATE TABLE prj.all_patid_death_end_dat AS
SELECT * FROM cprd_aurum_data.r_valid_date_lookup;

-- Step 2: Add the cprd_ddate column
ALTER TABLE prj.all_patid_death_end_dat ADD COLUMN cprd_ddate DATE;

UPDATE prj.all_patid_death_end_dat t1
JOIN cprd_aurum_data.patient t2 ON t1.patid = t2.patid
SET t1.cprd_ddate = t2.cprd_ddate;

-- Step 3: Add the gp_death_end_date column
ALTER TABLE prj.all_patid_death_end_dat ADD COLUMN gp_death_end_date DATE;

UPDATE prj.all_patid_death_end_dat
SET gp_death_end_date = CASE
    WHEN cprd_ddate IS NOT NULL AND gp_end_date IS NOT NULL THEN LEAST(cprd_ddate, gp_end_date)
    WHEN cprd_ddate IS NOT NULL THEN cprd_ddate
    WHEN gp_end_date IS NOT NULL THEN gp_end_date
    ELSE NULL
END;

-- Step 4: Create indexes for the new date table

CREATE INDEX x_patid_all_patid_death_end_dat
ON prj.all_patid_death_end_dat (patid);

CREATE INDEX x_gp_end_date_all_patid_death_end_dat
ON prj.all_patid_death_end_dat (gp_end_date);

CREATE INDEX x_gp_death_end_date_all_patid_death_end_dat
ON prj.all_patid_death_end_dat (gp_death_end_date);

CREATE INDEX x_cprd_ddate_all_patid_death_end_dat
ON prj.all_patid_death_end_dat (cprd_ddate);

CREATE INDEX x_min_dob_all_patid_death_end_dat
ON prj.all_patid_death_end_dat (min_dob);
