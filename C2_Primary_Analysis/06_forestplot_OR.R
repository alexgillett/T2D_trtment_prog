############################################################################################
############################################################################################
### Forest plots for ORs
############################################################################################
############################################################################################
### Sections in this script
############################################################################################
### Section 1: Set up
########################################
### Section 2: Preparing data for regression plotting (adjusted OR forest plot)
########################################
### Section 3: Generate plot
############################################################################################

############################################################################################
### Section 1: Setup
############################################################################################
library(survival)
library(survminer)
library(splines)
library(flexsurv)
library(dplyr)
library(tidyverse)
library(ggplot2)

### Data:
insulin_dat <-  readRDS("/path_to_model_output/na_pooled_depOR_insulin_YEARmi.rds")
intens_dat <- readRDS("/path_to_model_output/na_pooled_depOR_intens_YEARmi.rds")

############################################################################################
### Section 2: Preparing data for regression plotting (adjusted OR forest plot)
############################################################################################
### This code loads the pooled MI OR results for both outcomes, cleans and standardises the variable names, 
### recodes the depression exposure categories, sets the plotting order, and prepares a tidy combined dataset 
### for use in the final OR forest plot.

### Add an outcome label and standardises column names
adj_ORs_intens <- intens_dat %>%
  mutate(Outcome = "Intensification") %>%
  rename(upper = UpperCI, lower = LowerCI)

adj_ORs_ins <- insulin_dat %>%
  mutate(Outcome = "Insulin initiation") %>%
  rename(upper = UpperCI, lower = LowerCI)

### Combine both outcomes into one dataset
OR_data <- bind_rows(adj_ORs_intens, adj_ORs_ins)

### Recode the depression category variable into readable labels
OR_data <- OR_data %>%
  mutate(term = recode(Variable,
    "lastpre1oaddep_quantMDD <= 1.7yrs" = "Recent",
    "lastpre1oaddep_quantMDD 1.7-12.8yrs" = "Intermediate",
    "lastpre1oaddep_quantMDD > 12.8yrs" = "Distant"
  ))
### Set the factor order for depression categories
OR_data <- OR_data %>%
mutate(term = factor(term, levels = c(
    "Recent",
    "Intermediate",
    "Distant"
  )))

### Set the order of the outcomes
OR_data <- OR_data %>%
  mutate(Outcome = factor(Outcome, levels = c("Insulin initiation", "Intensification")))

### Sort the data
OR_data <- OR_data %>%
  arrange(term, Outcome)  

############################################################################################
### Section 3: Generate plot
############################################################################################

out_pathORplots <- "/path_to_plots/"

pdf(paste0(out_pathORplots, "MIorplotdeponly_intens_insulin.pdf"))
ggplot(OR_data, aes(x = OR, y = term, color = Outcome, group = Outcome)) +
  geom_vline(xintercept = 1, color = "gray75", linetype = "dashed") +
  geom_linerange(aes(xmin = lower, xmax = upper), size = 1.5, position = position_dodge(width = 0.5)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  theme_minimal() +
  scale_color_manual(values = c("Intensification" = "#0072B2", "Insulin initiation" = "#D55E00"),
  breaks = c("Intensification", "Insulin initiation")) +
  xlim(c(0.8, max(OR_data$upper) * 1.1)) +  # Auto-scale x-axis based on CI range
  labs(
       x = "Adjusted Odds Ratio (95% CI)", y = "Depression recency", color = "Outcome") +
  theme(
    legend.position = c(0.9, 0.9),  # Moves legend to the top-right
    legend.justification = c(1, 1),  # Anchors legend to the top-right
    legend.key.size = unit(0.4, "cm"),  # Reduce legend size
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )   + geom_text(
  aes(label = sprintf("%.2f", OR), color = Outcome),
  position = position_dodge(width = 0.5),
  vjust = -1,  # moves label above the CI
  size = 4,
  show.legend = FALSE
)
dev.off()