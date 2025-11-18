##########################################################################
##########################################################################
### Standardised survival plots for:
### (a) Treatment intensification and (b) Insulin initiation
##########################################################################
##########################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
########################################
### Section 2: Generate plotting data
########################################
### Section 3: Plot standardised survival curves (probability event-free)
########################################
### Section 4: Convert survival to cumulative incidence and plot again
########################################
### Section 5:
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################
library(survival)
library(survminer)
library(rms)
library(splines)
library(flexsurv)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(scales)

### Set path to plots
out_path <- "/path_to_plots/"

### Load pooled standardised survival results
pooled_intens <- readRDS("/path_to_model_output/mi_intens_pooled_surv.rds")
pooled_insulin <- readRDS("/path_to_model_output/mi_insulin_pooled_surv.rds")

############################################################################################
### Section 2: Generate plotting data
############################################################################################
### Combine outcomes into one plotting dataset and tidy factors
plot_dat <- bind_rows(
  pooled_intens  %>% mutate(outcome = "Treatment intensification"),
  pooled_insulin %>% mutate(outcome = "Insulin initiation")
) %>%
  mutate(
    group = factor(group, levels = c("None","Distant","Intermediate","Recent"))
  )
plot_dat <- plot_dat %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c("Treatment intensification", "Insulin initiation") 
    )
  )
plot_dat <- plot_dat %>%
  mutate(
    outcome = fct_recode(
      outcome,
      "(a) Treatment intensification" = "Treatment intensification",
      "(b) Insulin initiation"         = "Insulin initiation"
    )
  )

############################################################################################
### Section 3: Plot standardised survival curves (probability event-free)
############################################################################################

pdf(paste0(out_path, "pooled_standardised_survcurves.pdf"))
ggplot(plot_dat, aes(x = time, y = est, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ outcome, ncol = 1, scales = "fixed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Years since monotherapy initiation",
    y = "Survival (probability remaining event-free)",
    colour = "Depression recency",
    fill   = "Depression recency"
    #title  = "Risk-standardised survival by depression recency"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
dev.off()

############################################################################################
### Section 4: Convert survival to cumulative incidence and plot again
############################################################################################
plot_dat_ci <- plot_dat %>%
  transmute(
    time, group, outcome,
    est = 1 - est, lo = 1 - hi, hi = 1 - lo  # invert and keep CIs in order
  )

pdf(paste0(out_path, "pooled_standardised_cumincidence.pdf"))
ggplot(plot_dat_ci, aes(time, est, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, linewidth = 0) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ outcome, ncol = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Years since monotherapy initiation",
    y = "Cumulative incidence",
    colour = "Depression recency", fill = "Depression recency"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
dev.off()