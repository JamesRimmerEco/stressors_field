# Macrofauna analysis: community composition, stressor effects, and links to MPB Chl a
# This script merges multivariate analysis (NMDS, PERMANOVA) and univariate models linking MPB biomass (Chl a)
# to macrofaunal abundance and total counts.

# Load packages
library(here)
library(vegan)
library(ggplot2)
library(performance)
library(jtools)

################################################################################
#           Part 1: Community composition analysis (NMDS, PERMANOVA, SIMPER)   #
################################################################################

# Load macrofauna community data
macrofauna <- read.delim(here("Macrofauna.txt"))

# Treat missing values (NAs) as zeroes — these are true absences
macrofauna[is.na(macrofauna)] <- 0

# Convert variables to factors
macrofauna$Treatment <- as.factor(macrofauna$Treatment)
macrofauna$Glyphosate <- as.factor(macrofauna$Glyphosate)
macrofauna$TiO2 <- as.factor(macrofauna$TiO2)

# Subset species abundance data
species_abund <- macrofauna[, 7:17]

# Non-metric Multidimensional Scaling (NMDS)
set.seed(1)
nmds_model <- metaMDS(species_abund, k = 2, trymax = 100, trace = FALSE, autotransform = FALSE, distance = "bray")

# PERMANOVA (adonis2)
perm_model_main <- adonis2(species_abund ~ TiO2 + Glyphosate, data = macrofauna, permutations = 999, method = "bray")
perm_model_inter <- adonis2(species_abund ~ TiO2 * Glyphosate, data = macrofauna, permutations = 999, method = "bray")

# Homogeneity of dispersion assumption (betadisper)
distances <- vegdist(species_abund)
dispersion <- betadisper(distances, macrofauna$Treatment)
dispersion_test <- permutest(dispersion, pairwise = TRUE)

# SIMPER analysis
simper_gly <- simper(species_abund, group = macrofauna$Glyphosate)
simper_tio2 <- simper(species_abund, group = macrofauna$TiO2)

################################################################################
#           Part 2: Abundance–biomass link (MPB Chlorophyll a)                #
################################################################################

# Load macrofauna biomass + Chl a dataset
macrofauna_biomass <- read.delim(here("Macrofauna 2.txt"))

# Ensure numeric columns are treated correctly
numeric_cols <- c("Chlorophyll", "Polychaetes", "Oligochaetes", "Crustaceans", "Bivalves", "Individual.count")
macrofauna_biomass[numeric_cols] <- lapply(macrofauna_biomass[numeric_cols], as.numeric)

# GLMs and LM for abundance vs Chl a
mod_poly <- glm(Polychaetes ~ Chlorophyll, family = poisson(link = "log"), data = macrofauna_biomass)
mod_oli <- lm(Oligochaetes ~ Chlorophyll, data = macrofauna_biomass)
mod_crus <- lm(Crustaceans ~ Chlorophyll, data = macrofauna_biomass)
mod_biv <- glm(Bivalves ~ Chlorophyll, family = poisson(link = "log"), data = macrofauna_biomass)

# Two options for total count model (GLM and LM)
mod_ind_glm <- glm(Individual.count ~ Chlorophyll, family = poisson(link = "log"), data = macrofauna_biomass)
mod_ind_lm <- lm(Individual.count ~ Chlorophyll, data = macrofauna_biomass)

################################################################################
#           Part 3: Stressor effects on total macrofaunal abundance            #
################################################################################

# Subset first 42 rows where stressor levels are assigned
macrofauna_stress <- macrofauna_biomass[1:42, c("Glyphosate", "TiO2", "Individual.count")]
macrofauna_stress$Glyphosate <- as.factor(macrofauna_stress$Glyphosate)
macrofauna_stress$TiO2 <- as.factor(macrofauna_stress$TiO2)

# Test interaction and additive models
mod_stress_int <- lm(Individual.count ~ Glyphosate * TiO2, data = macrofauna_stress)
mod_stress_add <- lm(Individual.count ~ Glyphosate + TiO2, data = macrofauna_stress)

# Model comparison
AIC(mod_stress_int, mod_stress_add)
anova(mod_stress_int, mod_stress_add)
