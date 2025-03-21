# EPS: Colloidal carbohydrate analysis – pre, post, and adjusted values
# This script models colloidal carbohydrate content (EPS proxy) from contact core sediment samples
# collected before and after exposure to glyphosate and TiO₂ nanoparticles.
# Bayesian hierarchical models are fitted using rstanarm with weakly regularising priors.

# Load packages
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

# Split dataset into pre- and post-treatment subsets
contact_pre <- subset(contact_core, Stressor.status == "Pre")
contact_post <- subset(contact_core, Stressor.status == "Post")

# Compute adjusted difference in EPS (post - pre)
coll_diff <- contact_post
coll_diff$coll_diff <- contact_post$Colloidal.carbs - contact_pre$Colloidal.carbs

##################################
## Pre-treatment models
##################################

m1_pre_coll <- stan_lmer(
  Colloidal.carbs ~ Glyphosate + TiO2 + (1 | Patch),
  data = contact_pre,
  seed = 48,
  prior_intercept = normal(203, 138, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(290, 290, 290), autoscale = FALSE),
  prior_aux = exponential(0.018, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_pre_loo <- loo(m1_pre_coll, k_threshold = 0.7)

m1_pre_coll_w <- stan_lmer(
  Colloidal.carbs ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = contact_pre,
  seed = 48,
  prior_intercept = normal(203, 138, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(290, 290, 290, 290), autoscale = FALSE),
  prior_aux = exponential(0.018, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
m1_pre_loo_w <- loo(m1_pre_coll_w, k_threshold = 0.7)

loo_compare(m1_pre_loo, m1_pre_loo_w)

##################################
## Post-treatment models
##################################

m1_post_coll <- stan_lmer(
  Colloidal.carbs ~ Glyphosate + TiO2 + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(215, 252, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(530, 530, 500), autoscale = FALSE),
  prior_aux = exponential(0.01, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_col_post1 <- loo(m1_post_coll, k_threshold = 0.7)

m1_post_coll_w <- stan_lmer(
  Colloidal.carbs ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(215, 252, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(530, 530, 500, 500), autoscale = FALSE),
  prior_aux = exponential(0.01, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_col_post1_w <- loo(m1_post_coll_w, k_threshold = 0.7)

m2_post_coll <- stan_lmer(
  Colloidal.carbs ~ Glyphosate * TiO2 + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(215, 252, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0), c(530, 530, 500, 500, 500), autoscale = FALSE),
  prior_aux = exponential(0.01, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_col_post2 <- loo(m2_post_coll, k_threshold = 0.7)

m2_post_coll_w <- stan_lmer(
  Colloidal.carbs ~ Glyphosate * TiO2 + Percent.water + (1 | Patch),
  data = contact_post,
  seed = 48,
  prior_intercept = normal(215, 252, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0, 0), c(530, 530, 500, 500, 500, 500), autoscale = FALSE),
  prior_aux = exponential(0.01, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_col_post2_w <- loo(m2_post_coll_w, k_threshold = 0.7)

loo_compare(loo_col_post1, loo_col_post1_w, loo_col_post2, loo_col_post2_w)

##################################
## Adjusted (post - pre) models
##################################

m1_adju_coll <- stan_lmer(
  coll_diff ~ Glyphosate + TiO2 + (1 | Patch),
  data = coll_diff,
  seed = 48,
  prior_intercept = normal(12, 260, autoscale = FALSE),
  prior = normal(c(0, 0, 0), c(540, 540, 500), autoscale = FALSE),
  prior_aux = exponential(0.001, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_adj_coll1 <- loo(m1_adju_coll, k_threshold = 0.7)

m1_adju_coll_w <- stan_lmer(
  coll_diff ~ Glyphosate + TiO2 + Percent.water + (1 | Patch),
  data = coll_diff,
  seed = 48,
  prior_intercept = normal(12, 260, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0), c(540, 540, 500, 500), autoscale = FALSE),
  prior_aux = exponential(0.001, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_adj_coll1_w <- loo(m1_adju_coll_w, k_threshold = 0.7)

m2_adju_coll <- stan_lmer(
  coll_diff ~ Glyphosate * TiO2 + (1 | Patch),
  data = coll_diff,
  seed = 48,
  prior_intercept = normal(12, 260, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0), c(540, 540, 500, 500, 500), autoscale = FALSE),
  prior_aux = exponential(0.001, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_adj_coll2 <- loo(m2_adju_coll, k_threshold = 0.7)

m2_adju_coll_w <- stan_lmer(
  coll_diff ~ Glyphosate * TiO2 + Percent.water + (1 | Patch),
  data = coll_diff,
  seed = 48,
  prior_intercept = normal(12, 260, autoscale = FALSE),
  prior = normal(c(0, 0, 0, 0, 0, 0), c(540, 540, 500, 500, 500, 500), autoscale = FALSE),
  prior_aux = exponential(0.001, autoscale = FALSE),
  prior_covariance = decov(1, 1, 1, 1)
)
loo_adj_coll2_w <- loo(m2_adju_coll_w, k_threshold = 0.7)

loo_compare(loo_adj_coll1, loo_adj_coll1_w, loo_adj_coll2, loo_adj_coll2_w)
