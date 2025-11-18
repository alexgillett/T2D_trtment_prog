############################################################################################
############################################################################################
### Standardised survival curve code: Insulin
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
########################################
### Section 2: Generate standardised survival across MI datasets and pool results
########################################
### Section 3: Generate risk difference in standardised survival across MI datasets and pool results
########################################
### Section 4: Generate restricted mean survival time (RMST) and difference in RMST across MI datasets and pool results
########################################

############################################################################################
### Section 1: Setup
############################################################################################

library(survival)
library(survminer)
library(dplyr)
library(flexsurv)
library(mice)
library(mitools)

library(purrr)
library(tibble)

### ATTENTION: Make sure you set you working directory to where you want your output saved

########################################
### Get IDs for random subsample of 10,000 participants
########################################
### Only run once:
### Select row IDs
#dat_1 <- readRDS(paste0("mi_na_intens", "1", ".rds"))
#N <- dim(dat_1)[1]
#
#set.seed(202510)                      
#keep_idx <- sort(sample.int(N, 10000, replace = FALSE))
#saveRDS(keep_idx, "mi_keep_idx_10k.rds")

### Read in the 10,000 IDs
keep_idx <- readRDS("mi_keep_idx_10k.rds")
### Follow-up time points at which to calculate survival
times <- c(0, 0.05, seq(0.1, 0.5, by = 0.1), 0.7, 0.9, seq(1, 12, by = 0.25)) ### 54
### Follow-up time points at which to calculate risk difference (survival)
times_r <- c(1,3,5,10,12)
### Follow-up time points at which to calculate RMST
times_rmst <- c(1,5,10)

############################################################################################
### Section 2: Generate standardised survival 
############################################################################################

range <- 1:30

for (i in range) {
message("Running surv rmst ", i)

  model_i <- readRDS(paste0("mi_na_YEARrpm_insulin_", i, ".rds"))
  vars_needed <- all.vars(formula(model_i))
  dat_i <- readRDS(paste0("mi_na_insulin_", i, ".rds"))
  dat_i <- dat_i[keep_idx, ]
  dat_i <- dat_i[, intersect(names(dat_i), vars_needed), drop = FALSE]

  standsurv_i <- standsurv(model_i, newdata = dat_i, t=times, type = "survival", 
                    se = TRUE, ci=TRUE, boot=TRUE, B=300,
                    at = list(list(lastpre1oaddep_quant = "No MDD"),
                                                        list(lastpre1oaddep_quant = "MDD > 12.8yrs"),
                                                       list(lastpre1oaddep_quant = "MDD 1.7-12.8yrs"),
                                                       list(lastpre1oaddep_quant = "MDD <= 1.7yrs")))
  
  ### save survival data i
  saveRDS(standsurv_i, paste0("standsurv_mi_insulin_", i, ".rds"))
  rm(dat_i, model_i, standsurv_i)
}

########################################
### Read in the survival data across MI datasets
########################################
n_imps <- 30
survcurve_out <- list()
for(i in 1:n_imps){
    out_i <- readRDS(paste0("standsurv_mi_insulin_", i, ".rds"))
    survcurve_out[[i]] <- out_i
}

########################################
### Function: Rubin's rules for a single scalar estimate across M imputations
########################################
# Inputs:
#   q : numeric vector of length M with point estimates (one per imputation)
#   u : numeric vector of length M with within-imputation variances
# Output:
#   named vector with lower/estimate/upper CI, SE, and MI degrees of freedom
########################################
pool_scalar <- function(q, u) {
  M    <- length(q)                 # number of imputations
  qbar <- mean(q)                   # pooled point estimate
  Ubar <- mean(u)                   # mean within-imputation variance
  B    <- var(q)                    # between-imputation variance
  Tvar <- Ubar + (1 + 1/M) * B      # total variance
  r    <- ((1 + 1/M) * B) / Ubar    # relative increase in variance due to MI
  df   <- (M - 1) * (1 + 1/r)^2     # Barnard–Rubin degrees of freedom
  se   <- sqrt(Tvar)                # pooled standard error
  tcrt <- qt(0.975, df = df)        # two-sided 95% t critical value
  c(lo = qbar - tcrt * se,
    est = qbar,
    hi = qbar + tcrt * se,
    se = se,
    df = df)
}

########################################
### Function: Pool MI survival for each depression group and time
########################################
# Assumptions:
#   - mi_tbls is a list of tibbles with identical 'time' grids
#   - Each tibble has columns: time, at1, at1_se, ..., at4, at4_se
#   - at1..at4 are survival estimates in [0,1]; at*_se are their SEs
# Arguments:
#   mi_tbls : list of tibbles (one per imputation)
#   groups  : names corresponding to at1..at4 (customize as needed)
#   prefix  : the shared column prefix ("at")
#   eps     : small clamp to avoid 0/1 when logit-transforming
# Returns:
#   tibble with columns: time, group, est, lo, hi, se, df
#     where est/lo/hi are pooled survival probabilities
########################################
pool_mi_survival_wide <- function(mi_tbls,
                                  groups = c("None", "Distant", "Intermediate", "Recent"),
                                  prefix = "at",
                                  eps = 1e-8) {
  stopifnot(length(mi_tbls) >= 2)

  # Use the time grid from the first imputation; check others match
  time_vec <- mi_tbls[[1]]$time
  ntime    <- length(time_vec)
  M        <- length(mi_tbls)

  # Function to pool one group's estimates (e.g., at1 + at1_se) across imputations
  pool_one_group <- function(g_idx) {
    est_col <- paste0(prefix, g_idx)      # e.g., "at1"
    se_col  <- paste0(prefix, g_idx, "_se")  # e.g., "at1_se"

    # Pre-allocate matrices of size M x ntime for estimates and variances
    Q <- matrix(NA_real_, nrow = M, ncol = ntime)  # survival estimates
    U <- matrix(NA_real_, nrow = M, ncol = ntime)  # within-imputation variances

    # Collect survival + variance at each time for each imputation
    for (m in seq_len(M)) {
      dfm <- mi_tbls[[m]]
      # Ensure the time grid matches across imputations
      stopifnot(all(dfm$time == time_vec))

      # Point estimates and SEs for this group/this imputation
      q <- dfm[[est_col]]
      s <- dfm[[se_col]]

      # Clamp survival to (0,1) to stabilize logit transform
      Q[m, ] <- pmin(pmax(q, eps), 1 - eps)
      U[m, ] <- s^2
    }

    # ---- Pool on the logit scale for better normality/bounds handling ----
    logit <- function(p) log(p / (1 - p))

    # Delta method: Var(logit S) ≈ Var(S) / [S^2 * (1-S)^2]
    Q_tr <- logit(Q)
    U_tr <- U / (Q^2 * (1 - Q)^2)

    # Apply Rubin's rules per time point independently
    pooled_mat <- t(vapply(seq_len(ntime), function(j) {
      pool_scalar(q = Q_tr[, j], u = U_tr[, j])
    }, FUN.VALUE = numeric(5))) # lo, est, hi, se, df

    pooled <- as_tibble(pooled_mat) |>
      mutate(
        time = time_vec,
        # Back-transform CIs and estimate from logit to probability
        est = 1 / (1 + exp(-est)),
        lo  = 1 / (1 + exp(-lo)),
        hi  = 1 / (1 + exp(-hi)),
        group = groups[g_idx]
      ) |>
      select(time, group, est, lo, hi, se, df)

    pooled
  }

  # Bind pooled results for each depression group
  purrr::map_dfr(seq_along(groups), pool_one_group)
}

########################################
### Pool survival results
########################################
pooled_surv <- pool_mi_survival_wide(mi_tbls = survcurve_out)

dim(pooled_surv)
pooled_surv$lo[pooled_surv$time == 0] <- 1
pooled_surv$hi[pooled_surv$time == 0] <- 1

saveRDS(pooled_surv, file="mi_insulin_pooled_surv.rds")

rm(out_i)

############################################################################################
### Section 3: Generate risk difference in standardised survival
############################################################################################

for (i in range) {

message("Running survdiff ", i)

  model_i <- readRDS(paste0("mi_na_YEARrpm_insulin_", i, ".rds"))
  vars_needed <- all.vars(formula(model_i))
  dat_i <- readRDS(paste0("mi_na_insulin_", i, ".rds"))
  dat_i <- dat_i[keep_idx, ]
  dat_i <- dat_i[, intersect(names(dat_i), vars_needed), drop = FALSE]

standsurv_i <- standsurv(model_i, newdata = dat_i, t=times_r, type = "survival", 
                    se = TRUE, ci=TRUE, boot=TRUE, B=300, contrast = "difference",
                    at = list(list(lastpre1oaddep_quant = "No MDD"),
                                                        list(lastpre1oaddep_quant = "MDD > 12.8yrs"),
                                                       list(lastpre1oaddep_quant = "MDD 1.7-12.8yrs"),
                                                       list(lastpre1oaddep_quant = "MDD <= 1.7yrs")))
  

  ### save risk diff data i
  saveRDS(standsurv_i, paste0("standsurvdiff_mi_insulin_", i, ".rds"))
  rm(dat_i, model_i, standsurv_i)

}

########################################
### Read in all risk difference data across MI datasets
########################################
survdiff_out <- list()

for(i in 1:n_imps){
    out_i <- readRDS(paste0("standsurvdiff_mi_insulin_", i, ".rds"))
    survdiff_out[[i]] <- out_i
}

########################################
### Function: Pool contrasts function
########################################
# Pools on the identity scale (appropriate for differences in probabilities).
# mi_tbls: list of tibbles with columns:
#   time, contrast2_1, contrast2_1_se, contrast3_1, contrast3_1_se, contrast4_1, contrast4_1_se
# comp_map: named character vector: names = column suffixes, values = nice labels
#           e.g., c("2_1"="Distant - None","3_1"="Intermediate - None","4_1"="Recent - None")
########################################
pool_mi_contrasts <- function(mi_tbls,
                              comp_map = c("2_1"="Distant - None",
                                           "3_1"="Intermediate - None",
                                           "4_1"="Recent - None"),
                              prefix = "contrast") {
  stopifnot(length(mi_tbls) >= 2)
  time_vec <- mi_tbls[[1]]$time
  ntime    <- length(time_vec)
  M        <- length(mi_tbls)

  pool_one_contrast <- function(suffix) {
    est_col <- paste0(prefix, suffix)          # e.g., "contrast2_1"
    se_col  <- paste0(prefix, suffix, "_se")   # e.g., "contrast2_1_se"

    Q <- matrix(NA_real_, M, ntime)  # estimates (differences)
    U <- matrix(NA_real_, M, ntime)  # variances

    for (m in seq_len(M)) {
      dfm <- mi_tbls[[m]]
      stopifnot(all(dfm$time == time_vec))
      Q[m, ] <- dfm[[est_col]]
      U[m, ] <- (dfm[[se_col]])^2
    }

    pooled_mat <- t(vapply(seq_len(ntime), function(j) {
      pool_scalar(Q[, j], U[, j])
    }, FUN.VALUE = numeric(5)))  # lo, est, hi, se, df on identity scale

    as_tibble(pooled_mat) |>
      mutate(time = time_vec,
             comparison = unname(comp_map[[suffix]])) |>
      select(time, comparison, est, lo, hi, se, df)
  }

  map_dfr(names(comp_map), pool_one_contrast)
}

########################################
### Pool risk differences across MI datasets
########################################
comp_map <- c("2_1"="Distant - None",
              "3_1"="Intermediate - None",
              "4_1"="Recent - None")
pooled_contrasts <- pool_mi_contrasts(survdiff_out, comp_map = comp_map)

head(pooled_contrasts)
cbind(out_i$contrast2_1, out_i$contrast2_1_lci, out_i$contrast2_1_uci)

pooled_surv_restricted <- pool_mi_survival_wide(survdiff_out)

saveRDS(pooled_surv_restricted, "mi_insulin_pooled_survreduced.rds")
saveRDS(pooled_contrasts, "mi_insulin_pooled_riskdiff.rds")

############################################################################################
### Section 4: Generate restricted mean survival time (RMST)
############################################################################################

times_r <- c(1,3,5,10)
for (i in range) {

message("Running rmst ", i)

  model_i <- readRDS(paste0("mi_na_YEARrpm_insulin_", i, ".rds"))
  vars_needed <- all.vars(formula(model_i))
  dat_i <- readRDS(paste0("mi_na_insulin_", i, ".rds"))
  dat_i <- dat_i[keep_idx, ]
  dat_i <- dat_i[, intersect(names(dat_i), vars_needed), drop = FALSE]

  ### rmst...
  rmst_i <- standsurv(model_i, newdata = dat_i, t=times_r, type = "rmst", 
                     se = TRUE, ci=TRUE, boot=TRUE, B=100, contrast = "difference",
                    at = list(list(lastpre1oaddep_quant = "No MDD"),
                                                        list(lastpre1oaddep_quant = "MDD > 12.8yrs"),
                                                       list(lastpre1oaddep_quant = "MDD 1.7-12.8yrs"),
                                                       list(lastpre1oaddep_quant = "MDD <= 1.7yrs")))
  ### save survival data i
  saveRDS(rmst_i, paste0("rmst_mi_insulin_", i, ".rds"))

  rm(dat_i, model_i, rmst_i)
}

########################################
### Read in RMST data from across MI datasets
########################################
rmst_out <- list()

for(i in 1:n_imps){
    out_i <- readRDS(paste0("rmst_mi_insulin_", i, ".rds"))
    rmst_out[[i]] <- out_i
}

########################################
### Function: Pool RMST results
########################################
pool_mi_rmst <- function(mi_tbls,
                         groups = c("None","Distant","Intermediate","Recent"),
                         prefix = "at") {
  stopifnot(length(mi_tbls) >= 2)
  time_vec <- mi_tbls[[1]]$time
  ntime    <- length(time_vec)
  M        <- length(mi_tbls)

  pool_one_group <- function(g_idx) {
    est_col <- paste0(prefix, g_idx)        # e.g., "at1"
    se_col  <- paste0(prefix, g_idx, "_se") # e.g., "at1_se"

    Q <- matrix(NA_real_, M, ntime)   # RMST estimates (years)
    U <- matrix(NA_real_, M, ntime)   # variances

    for (m in seq_len(M)) {
      dfm <- mi_tbls[[m]]
      stopifnot(all(dfm$time == time_vec))
      Q[m, ] <- dfm[[est_col]]
      U[m, ] <- (dfm[[se_col]])^2
    }

    pooled_mat <- t(vapply(seq_len(ntime), function(j) {
      pool_scalar(Q[, j], U[, j])
    }, FUN.VALUE = numeric(5)))  # lo, est, hi, se, df

    as_tibble(pooled_mat) |>
      mutate(time = time_vec,
             group = groups[g_idx]) |>
      select(time, group, est, lo, hi, se, df)
  }

  map_dfr(seq_along(groups), pool_one_group)
}

########################################
### Function: Pool difference in RMST
########################################
pool_mi_rmst_contrasts <- function(mi_tbls,
                                   comp_map = c("2_1"="Distant - None",
                                                "3_1"="Intermediate - None",
                                                "4_1"="Recent - None"),
                                   prefix = "contrast") {
  stopifnot(length(mi_tbls) >= 2)
  time_vec <- mi_tbls[[1]]$time
  ntime    <- length(time_vec)
  M        <- length(mi_tbls)

  pool_one_contrast <- function(suffix) {
    est_col <- paste0(prefix, suffix)
    se_col  <- paste0(prefix, suffix, "_se")

    Q <- matrix(NA_real_, M, ntime)   # RMST differences (years)
    U <- matrix(NA_real_, M, ntime)   # variances

    for (m in seq_len(M)) {
      dfm <- mi_tbls[[m]]
      stopifnot(all(dfm$time == time_vec))
      Q[m, ] <- dfm[[est_col]]
      U[m, ] <- (dfm[[se_col]])^2
    }

    pooled_mat <- t(vapply(seq_len(ntime), function(j) {
      pool_scalar(Q[, j], U[, j])
    }, FUN.VALUE = numeric(5)))

    as_tibble(pooled_mat) |>
      mutate(time = time_vec,
             comparison = unname(comp_map[[suffix]])) |>
      select(time, comparison, est, lo, hi, se, df)
  }

  map_dfr(names(comp_map), pool_one_contrast)
}

########################################
### Pool RMST results
########################################
rmst_pooled <- pool_mi_rmst(rmst_out)
rmst_diff_pooled <- pool_mi_rmst_contrasts(rmst_out)

saveRDS(rmst_pooled, "mi_insulin_pooled_rmst.rds")
saveRDS(rmst_diff_pooled, "mi_insulin_pooled_rmstdiff.rds")

### RMST in months
pooled_rmst_months <- rmst_diff_pooled %>%
  mutate(
    est = est * 12,
    lo  = lo * 12,
    hi  = hi * 12
  )
