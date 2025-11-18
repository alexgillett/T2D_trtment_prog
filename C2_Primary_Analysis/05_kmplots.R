############################################################################################
############################################################################################
### KM plots for both outcomes (strata = depression recency)
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
########################################
### Section 2: KM plot for time to treatment intensification
########################################
### Section 3: KM plot for time to insulin initiation
########################################
### Section 4: Combined plot
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################
### Libraries:
library(flexsurv)
library(survival)
library(ggplot2)
library(dplyr)
library(forcats)
library(survminer)

### Data:
med_ds1 <- readRDS(file="path_to_your_secure_area/analysis_dsyyyymmdd.rds")

### Recode the depression recency categories into cleaner labels
med_ds1 <- med_ds1 %>%
  mutate(lastpre1oaddep_quant = fct_recode(lastpre1oaddep_quant,
    "None"         = "No MDD",
    "Recent"       = "MDD <= 1.7yrs",
    "Intermediate" = "MDD 1.7-12.8yrs",
    "Distant"      = "MDD > 12.8yrs"
  ))

### Check levels
levels(med_ds1$lastpre1oaddep_quant)

### Sets the order of the factor levels
med_ds1 <- med_ds1 %>%
  mutate(lastpre1oaddep_quant = factor(lastpre1oaddep_quant,
    levels = c("None", "Distant", "Intermediate", "Recent")))

med_ds1$Dep <- med_ds1$lastpre1oaddep_quant


############################################################################################
### Section 2: KM plot for time to treatment intensification
############################################################################################
### KM stratified by depression recency
km_fit_intensfull <- survfit(Surv(time_change_yrs, intens_event) ~ Dep,
                  data = med_ds1)

### Rename strata
names(km_fit_intensfull$strata) <- c("None", "Distant", "Intermediate", "Recent")

### Plot KM
p3 <- ggsurvplot(
  km_fit_intensfull,
  data = med_ds1,
  risk.table = TRUE,
  xlab = "Years since monotherapy start",
  ylab = "Survival (no intensification)",
  legend.title = "Depression",
  tables.y.text = FALSE
)

###print(p3)

### Save as PDF
pdf("/path_to_plots/intens_kmplot_dep.pdf",
    width = 7, height = 6)   # sizes are in inches for PDF
print(p3)
dev.off()

############################################################################################
### Section 3: KM plot for time to insulin initiation
############################################################################################
### KM stratified by depression recency
km_fit_insulin <- survfit(Surv(time_insulin_yrs, ins_ever) ~ Dep,
                  data = med_ds1)

### Rename strata
names(km_fit_insulin$strata) <- c("None", "Distant", "Intermediate", "Recent")

### Plot KM
p4 <- ggsurvplot(
  km_fit_insulin,
  data = med_ds1,
  risk.table = TRUE,
  xlab = "Years since monotherapy start",
  ylab = "Survival (no insulin)",
  legend.title = "Depression",
  tables.y.text = FALSE
)

###print(p4)

### Save as PDF
pdf("/path_to_plots/insulin_kmplot_dep.pdf",
    width = 7, height = 6)   # sizes are in inches for PDF
print(p4)
dev.off()

############################################################################################
### Section 4: Combined plot
############################################################################################
### Both intens and insulin together...
library(patchwork)
fit_intens  <- km_fit_intensfull
fit_insulin <- km_fit_insulin

# Plot options
leg_title <- "Depression recency"
leg_labs  <- levels(med_ds1$lastpre1oaddep_quant)

# Shared theme/ style
km_theme <- theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
    plot.title = element_text(face = "bold", size = 10)
  )

# Create plots without risk tables
p1 <- ggsurvplot(
  fit_intens,
  data = med_ds1,
  conf.int = FALSE,
  censor = FALSE,
  palette = "Set1",
  legend.title = "Depression recency",
  legend.labs = levels(med_ds1$lastpre1oaddep_quant),
  xlab = "",
  ylab = "Proportion without treatment intensification",
  title = "(a) Treatment intensification",
  ggtheme = km_theme,
  risk.table = FALSE
)$plot

p2 <- ggsurvplot(
  fit_insulin,
  data = med_ds1,
  conf.int = FALSE,
  censor = FALSE,
  palette = "Set1",
  legend.title = "Depression recency",
  legend.labs = levels(med_ds1$lastpre1oaddep_quant),
  xlab = "Years since monotherapy initiation",
  ylab = "Proportion without insulin initiation",
  title = "(b) Insulin initiation",
  ggtheme = km_theme,
  risk.table = FALSE
)$plot

# Combine vertically using patchwork 
combined_plot <- p1 / p2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save 
pdf(file.path(out_path, "KM_intens_and_insulin_combined_clean.pdf"))
print(combined_plot)
dev.off()