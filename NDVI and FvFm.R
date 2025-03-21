
# Bayesian models of Fv/Fm and NDVI across timepoints
# This script uses rstanarm with timepoint-specific data subsets, fits interaction and additive models, compares them with LOO, and reports summaries

# Load packages
library(here)
library(lme4)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(bayestestR)
library(tidyverse)

# Load data
field_dat <- read.delim(here("Processed data for R.txt"))

# Format variables
field_dat <- field_dat %>%
  mutate(
    Treatment = as.factor(Treatment),
    Glyphosate = as.factor(Glyphosate),
    TiO2 = as.factor(TiO2),
    Experiment.ID = as.factor(Experiment.ID),
    stress_applied = as.numeric(Time > 1),
    NDVI = ifelse(NDVI < 0, NA, NDVI)
  )

# Split data by timepoint
split_by_time <- split(field_dat, field_dat$Time)

# Function to fit and evaluate models
fit_and_report <- function(timepoint, data) {
  message(paste0("----- Timepoint ", timepoint, " -----"))
  
  # --- Fv/Fm ---
  m_fv_int <- stan_lmer(FvFm ~ Glyphosate * TiO2 + (1|Patch),
                        data = data, seed = 38,
                        prior_intercept = normal(0, 0.32, autoscale = FALSE),
                        prior = normal(0, 0.32, autoscale = FALSE),
                        prior_aux = exponential(1),
                        prior_covariance = decov(1,1,1,1))
  m_fv_add <- stan_lmer(FvFm ~ Glyphosate + TiO2 + (1|Patch),
                        data = data, seed = 38,
                        prior_intercept = normal(0, 0.32, autoscale = FALSE),
                        prior = normal(0, 0.32, autoscale = FALSE),
                        prior_aux = exponential(1),
                        prior_covariance = decov(1,1,1,1))
  
  loo_fv_int <- loo(m_fv_int, k_threshold = 0.7)
  loo_fv_add <- loo(m_fv_add, k_threshold = 0.7)
  comp_fv <- loo_compare(loo_fv_int, loo_fv_add)
  
  best_fv <- if (comp_fv[1, 1] < 0) m_fv_int else m_fv_add
  print(best_fv, digits = 3)
  hdi(best_fv, ci = 0.89)
  pp_check(best_fv)
  describe_posterior(best_fv)
  
  # --- NDVI ---
  m_ndvi_int <- stan_lmer(NDVI ~ Glyphosate * TiO2 + (1|Patch),
                          data = data, seed = 38,
                          prior_intercept = normal(0, 0.32, autoscale = FALSE),
                          prior = normal(0, 0.4, autoscale = FALSE),
                          prior_aux = exponential(1),
                          prior_covariance = decov(1,1,1,1))
  m_ndvi_add <- stan_lmer(NDVI ~ Glyphosate + TiO2 + (1|Patch),
                          data = data, seed = 38,
                          prior_intercept = normal(0, 0.32, autoscale = FALSE),
                          prior = normal(0, 0.4, autoscale = FALSE),
                          prior_aux = exponential(1),
                          prior_covariance = decov(1,1,1,1))
  
  loo_ndvi_int <- loo(m_ndvi_int, k_threshold = 0.7)
  loo_ndvi_add <- loo(m_ndvi_add, k_threshold = 0.7)
  comp_ndvi <- loo_compare(loo_ndvi_int, loo_ndvi_add)
  
  best_ndvi <- if (comp_ndvi[1, 1] < 0) m_ndvi_int else m_ndvi_add
  print(best_ndvi, digits = 3)
  hdi(best_ndvi, ci = 0.89)
  pp_check(best_ndvi)
  describe_posterior(best_ndvi)
}

# Loop through timepoints T1–T7
for (i in 1:7) {
  fit_and_report(i, split_by_time[[i]])
}
