# Total carbohydrates analysis – pre, post, and adjusted values
# This script models sediment total carbohydrate concentrations from contact core sediment samples
# collected before and after exposure to glyphosate and TiO₂ nanoparticles.
# Bayesian hierarchical models are fitted using rstanarm with weakly regularising priors.

# Load packagess
library(plyr)
library(here)
library(tidyverse)
library(lme4)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(bayestestR)

# Load processed data
contact_core <- read.delim(here("contact_core_data.txt"))

# Convert treatment variables to factors
contact_core$Stressor.status <- as.factor(contact_core$Stressor.status)
contact_core$Glyphosate <- as.factor(contact_core$Glyphosate)
contact_core$TiO2 <- as.factor(contact_core$TiO2)

# Split into pre- and post-treatment subsets
contact_pre <- subset(contact_core, Stressor.status == "Pre")
contact_post <- subset(contact_core, Stressor.status == "Post")

# Compute adjusted difference in total carbohydrates (post - pre)
tot_diff <- contact_post
tot_diff$tot_diff <- contact_post$Total.carbs - contact_pre$Total.carbs

##############################################
## Pre-treatment models
##############################################

m1_pre_tot <- stan_lmer(
  Total.carbs ~ Glyphosate + TiO2 + (1 | Patch),
  data = contact_pre,
  seed = 48,
  prior_intercept = normal(4682, 2958, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(6200, 6200, 5900), autoscale = FALSE),
  prior_aux = exponential(0.00085, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_pre_tot_loo <- loo(m1_pre_tot, k_threshold = 0.7)

m1_pre_tot_w <- stan_lmer(
  Total.carbs ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = contact_pre,
  seed = 48,
  prior_intercept = normal(4682, 2958, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(6200, 6200, 5900, 5900), autoscale = FALSE),
  prior_aux = exponential(0.00085, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_pre_tot_w_loo <- loo(m1_pre_tot_w, k_threshold = 0.7)

loo_compare(m1_pre_tot_loo, m1_pre_tot_w_loo)

##############################################
## Post-treatment models
##############################################

m1_post_tot <- stan_lmer(
  Total.carbs ~ Glyphosate + TiO2 + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(5006, 4586, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(8100, 8100, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00055, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_post_tot_loo <- loo(m1_post_tot, k_threshold = 0.7)

m1_post_tot_w <- stan_lmer(
  Total.carbs ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(5006, 4586, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(8100, 8100, 8000, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00055, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_post_tot_w_loo <- loo(m1_post_tot_w, k_threshold = 0.7)

m2_post_tot <- stan_lmer(
  Total.carbs ~ Glyphosate * TiO2 + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(5006, 4586, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0), c(8100, 8100, 8000, 8000, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00055, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m2_post_tot_loo <- loo(m2_post_tot)

m2_post_tot_w <- stan_lmer(
  Total.carbs ~ Glyphosate * TiO2 + Percent.water + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(5006, 4586, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0, 0), c(8100, 8100, 8000, 8000, 8000, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00055, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m2_post_tot_w_loo <- loo(m2_post_tot_w, k_threshold = 0.7)

loo_compare(m1_post_tot_loo, m1_post_tot_w_loo, m2_post_tot_loo, m2_post_tot_w_loo)

##############################################
## Adjusted (post - pre) models
##############################################

m1_diff_tot <- stan_lmer(
  tot_diff ~ Glyphosate + TiO2 + (1 | Patch),
  data = tot_diff,
  seed = 48,
  prior_intercept = normal(320, 5600, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(11000, 11000, 11000), autoscale = FALSE),
  prior_aux = exponential(0.00044, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_diff_tot_loo <- loo(m1_diff_tot, k_threshold = 0.7)

m1_diff_tot_w <- stan_lmer(
  tot_diff ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = tot_diff,
  seed = 48,
  prior_intercept = normal(320, 5600, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(11000, 11000, 11000, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00044, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_diff_tot_w_loo <- loo(m1_diff_tot_w, k_threshold = 0.7)

m2_diff_tot <- stan_lmer(
  tot_diff ~ Glyphosate * TiO2 + (1 | Patch),
  data = tot_diff,
  seed = 48,
  prior_intercept = normal(320, 5600, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0), c(11000, 11000, 11000, 11000, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00044, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m2_diff_tot_loo <- loo(m2_diff_tot, k_threshold = 0.7)

m2_diff_tot_w <- stan_lmer(
  tot_diff ~ Glyphosate * TiO2 + Percent.water + (1 | Patch),
  data = tot_diff,
  seed = 48,
  prior_intercept = normal(5006, 4586, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0, 0), c(11000, 11000, 11000, 11000, 8000, 8000), autoscale = FALSE),
  prior_aux = exponential(0.00044, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m2_diff_tot_w_loo <- loo(m2_diff_tot_w, k_threshold = 0.7)

loo_compare(m1_diff_tot_loo, m1_diff_tot_w_loo, m2_diff_tot_loo, m2_diff_tot_w_loo)
