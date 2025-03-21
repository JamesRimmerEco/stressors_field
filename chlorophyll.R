# Chlorophyll a analysis – pre, post, and adjusted values
# This script models sediment chlorophyll a concentrations from contact core sediment samples
# collected before and after exposure to glyphosate and TiO₂ nanoparticles (TiO2 NP).
# Bayesian hierarchical models are fitted using rstanarm.

# Load packages
library(plyr)
library(here)
library(tidyverse)
library(lme4)
library(rstanarm)
library(shinystan)
library(bayesplot)
library(bayestestR)

# Load data
contact_core <- read.delim(here("contact_core_data.txt"))

# Factorise key variables
contact_core$Stressor.status <- as.factor(contact_core$Stressor.status)
contact_core$Glyphosate <- as.factor(contact_core$Glyphosate)
contact_core$TiO2 <- as.factor(contact_core$TiO2)

# Split data into pre/post exposure
contact_pre <- subset(contact_core, Stressor.status == "Pre")
contact_post <- subset(contact_core, Stressor.status == "Post")

# Calculate adjusted change in chlorophyll (post - pre)
chlor_diff <- contact_post
chlor_diff$chlor_diff <- contact_post$Chlorophyll - contact_pre$Chlorophyll

##############################################
## Pre-exposure models
##############################################

m1_pre <- stan_lmer(
  Chlorophyll ~ Glyphosate + TiO2 + (1 | Patch),
  data = contact_pre,
  seed = 48,
  prior_intercept = normal(14, 5.4, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(11.38, 11.38, 11.38), autoscale = FALSE),
  prior_aux = exponential(0.46, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_pre_loo <- loo(m1_pre, k_threshold = 0.7)

m_w_pre <- stan_lmer(
  Chlorophyll ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = contact_pre,
  seed = 48,
  adapt_delta = 0.99,
  prior_intercept = normal(14, 5.4, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(11.38, 11.38, 11.38, 11.38), autoscale = FALSE),
  prior_aux = exponential(0.46, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
mw_pre_loo <- loo(m_w_pre, k_threshold = 0.7)

loo_compare(m1_pre_loo, mw_pre_loo)

# Posterior draws for pre-stress models
chlor_pre_gly1 <- as.matrix(m_w_pre, pars = "Glyphosate1")
chlor_pre_gly2 <- as.matrix(m_w_pre, pars = "Glyphosate2")
chlor_pre_tio2 <- as.matrix(m_w_pre, pars = "TiO21")

mean(chlor_pre_gly1 < 0)
mean(chlor_pre_gly2 > 0)
mean(chlor_pre_tio2 < 0)

##############################################
## Post-exposure models
##############################################

m1_post <- stan_lmer(
  Chlorophyll ~ Glyphosate + TiO2 + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(20, 9.1, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(19, 19, 19), autoscale = FALSE),
  prior_aux = exponential(0.27, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_1_post <- loo(m1_post)

m1_post_w <- stan_lmer(
  Chlorophyll ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(20, 9.1, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(19, 19, 19, 19), autoscale = FALSE),
  prior_aux = exponential(0.27, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_1w_post <- loo(m1_post_w, k_threshold = 0.7)

# Posterior summaries
chlor_samp_gly1 <- as.matrix(m1_post_w, pars = "Glyphosate1")
chlor_samp_gly2 <- as.matrix(m1_post_w, pars = "Glyphosate2")
chlor_samp_tio2 <- as.matrix(m1_post_w, pars = "TiO21")

mean(chlor_samp_gly1 < 0)
mean(chlor_samp_gly2 < 0)
mean(chlor_samp_tio2 > 0)

# Interaction model
m2_post <- stan_lmer(
  Chlorophyll ~ Glyphosate * TiO2 + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(20, 9.1, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0), c(19, 19, 19, 19, 19), autoscale = FALSE),
  prior_aux = exponential(0.27, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_2_post <- loo(m2_post, k_threshold = 0.7)

m2_post_w <- stan_lmer(
  Chlorophyll ~ Glyphosate * TiO2 + Percent.water + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(20, 9.1, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0, 0), c(19, 19, 19, 19, 19, 19), autoscale = FALSE),
  prior_aux = exponential(0.27, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_2_post_w <- loo(m2_post_w, k_threshold = 0.7)

loo_compare(loo_1_post, loo_1w_post, loo_2_post, loo_2_post_w)

##############################################
## Adjusted chlorophyll difference models
##############################################

m1_adjust <- stan_lmer(
  chlor_diff ~ Glyphosate + TiO2 + (1 | Patch),
  data = chlor_diff,
  seed = 48,
  prior_intercept = normal(5.6, 11, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(24, 24, 24), autoscale = FALSE),
  prior_aux = exponential(0.22, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_1_adj <- loo(m1_adjust)

m1_adjust_w <- stan_lmer(
  chlor_diff ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = chlor_diff,
  seed = 48,
  prior_intercept = normal(5.6, 11, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(24, 24, 24, 24), autoscale = FALSE),
  prior_aux = exponential(0.22, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_1w_adj <- loo(m1_adjust_w, k_threshold = 0.7)

# Interaction model for adjusted chlorophyll
m2_adjust <- stan_lmer(
  chlor_diff ~ Glyphosate * TiO2 + (1 | Patch),
  data = chlor_diff,
  seed = 48,
  prior_intercept = normal(5.6, 11, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0), c(24, 24, 24, 24, 24), autoscale = FALSE),
  prior_aux = exponential(0.22, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_2_adj <- loo(m2_adjust)

m2_adjust_w <- stan_lmer(
  chlor_diff ~ Glyphosate * TiO2 + Percent.water + (1 | Patch),
  data = chlor_diff,
  seed = 48,
  prior_intercept = normal(5.6, 11, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0, 0), c(24, 24, 24, 24, 24, 24), autoscale = FALSE),
  prior_aux = exponential(0.22, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_2w_adj <- loo(m2_adjust_w, k_threshold = 0.7)

loo_compare(loo_1_adj, loo_1w_adj, loo_2_adj, loo_2w_adj)

# Posterior probabilities for adjusted interaction model
chlor_adj_gly1 <- as.matrix(m2_adjust, pars = "Glyphosate1")
chlor_adj_gly2 <- as.matrix(m2_adjust, pars = "Glyphosate2")
chlor_adj_tio2 <- as.matrix(m2_adjust, pars = "TiO21")
chlor_adj_gly1ti <- as.matrix(m2_adjust, pars = "Glyphosate1:TiO21")
chlor_adj_gly2ti <- as.matrix(m2_adjust, pars = "Glyphosate2:TiO21")

mean(chlor_adj_gly1 < 0)
mean(chlor_adj_gly2 < 0)
mean(chlor_adj_tio2 > 0)
mean(chlor_adj_gly1ti < 0)
mean(chlor_adj_gly2ti > 0)

# Group mean changes (manual calculation)
chlor_diff$Treatment <- as.factor(chlor_diff$Treatment)

mean(subset(chlor_diff, Treatment == "a")$chlor_diff)
mean(subset(chlor_diff, Treatment == "b")$chlor_diff)
mean(subset(chlor_diff, Treatment == "c")$chlor_diff)
mean(subset(chlor_diff, Treatment == "d")$chlor_diff)
mean(subset(chlor_diff, Treatment == "e")$chlor_diff)
mean(subset(chlor_diff, Treatment == "f")$chlor_diff)
