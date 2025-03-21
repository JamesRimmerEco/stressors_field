# Erosion modelling – treatment effects on erosion type and shear stress thresholds
# This script analyses erosion responses from contact core samples subjected to glyphosate and TiO₂ NP stressors.
# Models include categorical erosion type frequencies (Type I, I.5, II) and ordinal/logistic models of shear stress thresholds.

# Load packages
library(plyr)
library(here)
library(tidyverse)
library(lme4)
library(MASS)        # For polr()
library(DHARMa)      # For residual checks
library(brant)       # For proportional odds test
library(broom)

# Load formatted erosion dataset
CSM.dat <- read.delim(here("CSM formatted.txt"))

# Clean and format variables
CSM.dat$Treatment <- as.factor(CSM.dat$Treatment)
CSM.dat$Stressor.status <- as.factor(CSM.dat$Stressor.status)
CSM.dat$Erosion.type <- as.factor(CSM.dat$Erosion.type)
CSM.dat$Glyphosate <- as.factor(CSM.dat$Glyphosate)
CSM.dat$TiO2 <- as.factor(CSM.dat$TiO2)

# Subset by stressor exposure status
CSM.pre <- subset(CSM.dat, Stressor.status == "Pre")
CSM.post <- subset(CSM.dat, Stressor.status == "Post")

############################################################
# Section 0: Descriptive erosion type frequencies
############################################################

# Counts of erosion types before and after exposure
table(CSM.pre$Erosion.type)
table(CSM.post$Erosion.type)

# Cross-tabulations by treatment factors (pre-exposure)
table(CSM.pre$Erosion.type, CSM.pre$Glyphosate)
table(CSM.pre$Erosion.type, CSM.pre$TiO2)

# Three-way breakdowns of erosion types (pre/post)
with(CSM.pre, table(Erosion.type, Glyphosate, TiO2))
with(CSM.post, table(Erosion.type, Glyphosate, TiO2))

############################################################
# Section 1: Erosion type frequencies (Type I, I.5, II)
############################################################

# Raw erosion type counts by treatment group (manually coded)
pois_erosion <- as.data.frame(matrix(nrow = 12, ncol = 0))
pois_erosion$Type_I <- c(6,6,6,7,7,5, 5,6,5,2,5,4)
pois_erosion$Type_II <- c(0,1,1,0,0,2,1,1,1,3,2,3)
pois_erosion$Type_I.5 <- c(0,0,0,0,0,0, 1,0,1,2,0,0)
pois_erosion$TiO2 <- c(0,0,0,1,1,1, 0,0,0,1,1,1)
pois_erosion$Glyphosate <- factor(c(0,1,2,0,1,2, 0,1,2,0,1,2))
pois_erosion$stressor_status <- c(0,0,0,0,0,0, 1,1,1,1,1,1)
pois_erosion$Treatment <- factor(c("a", "b", "c", "a", "b", "c", "d", "e", "f", "d", "e", "f"))

# Poisson GLMs for each erosion type
TI <- glm(Type_I ~ stressor_status + TiO2 + Glyphosate, family = poisson(link = "log"), data = pois_erosion)
summary(TI)

TII <- glm(Type_II ~ Glyphosate + TiO2 + stressor_status, family = poisson(link = "log"), data = pois_erosion)
summary(TII)

TI.5a <- glm(Type_I.5 ~ Glyphosate + TiO2, family = poisson(link = "log"), data = pois_erosion)
summary(TI.5a)

TI.5b <- glm(Type_I.5 ~ stressor_status, family = poisson(link = "log"), data = pois_erosion)
summary(TI.5b)

# Note: Type I.5 is sparse and may require zero-inflated modelling

############################################################
# Section 2: Shear stress threshold modelling (Type I only)
############################################################

# Subset pre/post samples with erosion type I only
csm.pre.t1 <- subset(CSM.pre, Erosion.type == "1")
csm.post.t1 <- subset(CSM.post, Erosion.type == "1")

csm.pre.t1$TiO2 <- as.factor(csm.pre.t1$TiO2)
csm.pre.t1$Shear.stress_2 <- as.factor(csm.pre.t1$Shear.stress)

csm.post.t1$TiO2 <- as.factor(csm.post.t1$TiO2)
csm.post.t1$Shear.stress_2 <- as.factor(csm.post.t1$Shear.stress)

# Ordinal logistic regression (pre-stress)
pre_ss_a <- polr(Shear.stress_2 ~ Glyphosate + TiO2, data = csm.pre.t1, Hess = TRUE)
pre_ss_b <- polr(Shear.stress_2 ~ Glyphosate * TiO2, data = csm.pre.t1, Hess = TRUE)
pre_ss_c <- polr(Shear.stress_2 ~ Glyphosate + TiO2 + perc.water, data = csm.pre.t1, Hess = TRUE)

# Model selection
anova(pre_ss_a, pre_ss_b, pre_ss_c)
AIC(pre_ss_a, pre_ss_b, pre_ss_c)

# Summary and diagnostics
summary(pre_ss_a)
brant(pre_ss_a)  # Proportional odds assumption test
exp(cbind(OR = coef(pre_ss_a), confint(pre_ss_a)))

# Post-stress logistic regression (binary response)
post_ss_a <- glm(Shear.stress_2 ~ Glyphosate + TiO2, data = csm.post.t1, family = binomial)
post_ss_b <- glm(Shear.stress_2 ~ Glyphosate * TiO2, data = csm.post.t1, family = binomial)
post_ss_c <- glm(Shear.stress_2 ~ Glyphosate + TiO2 + perc.water, data = csm.post.t1, family = binomial)

# Model comparisons
anova(post_ss_a, post_ss_b, test = "Chisq")
anova(post_ss_a, post_ss_c, test = "Chisq")
AIC(post_ss_a, post_ss_b, post_ss_c)

# Summary and effect sizes
summary(post_ss_c)
exp(cbind(OR = coef(post_ss_c), confint(post_ss_c)))
