# README: Response of natural estuarine microphytobenthic biofilms to multiple anthropogenic stressors.

This repository accompanies the manuscript *Response of natural estuarine microphytobenthic biofilms to multiple anthropogenic stressors*, which investigates the effects of glyphosate and titanium dioxide nanoparticles (TiO₂ NP) on MPB communities in a temperate estuary. The study includes non-destructive monitoring, destructive sampling, macrofauna surveys, and modelling using Bayesian and frequentist techniques in R.

---

## Repository Structure

### Data

- `Processed data for R.txt` — Main experimental data file. 294 observations of Fv/Fm and NDVI over 7 timepoints.
- `Processed data for R.xlsx` — Excel version of the above.
- `Macrofauna.txt` — Macrofauna species data with 42 observations across 17 variables.
- `Macrofauna 2.txt` — Macrofauna species and derived metrics including richness and abundance. 84 observations.
- `contact_core_data.txt` — Sediment biotic variables (chlorophyll, carbohydrates, water content).
- `CSM formatted.txt` — Cohesive Strength Meter data, erosion thresholds, shear stress and erosion type.
- `Logger_1.xlsx`, `Logger_7.xlsx`, `Moleman_1.xlsx`, `Moleman_2.xlsx` — Light and temperature logger outputs. Each has 3,700+ rows with timestamp, temperature, and lux.

### Supplementary Figure S1

Photograph of the experimental field site setup with labels (in gitignore, available with the paper). 

### Supplementary Figure S2
Bar plot of sediment shear stress thresholds (Type I erosion) by treatment following exposure, 
generated using shear_stress_models.R. This figure corresponds to Supplementary Figure 2 in the 
manuscript and shows means and standard deviations for each glyphosate and TiO₂ NP treatment combination
(in gitignore).

### Supplementary Figure S3

nMDS ordination plot of macrofauna community structures (`Figure_S3_NMDS_macrofauna.png`) showing 
relative position of experimental units based on Bray-Curtis dissimilarity. Generated using `NMDS plots.R`. 
Corresponds to analysis described in the manuscript macrofauna section (in gitignore). 


### Scripts

- `NDVI and FvFm.R` — Fits Bayesian mixed models to Fv/Fm and NDVI at each timepoint.
- `environmental_models.R` — Assesses environmental drivers of NDVI and Fv/Fm using GLMMs.
- `chlorophyll.R` — Models chlorophyll a content and interaction with water content.
- `coll_carbohydrates.R` — Models colloidal carbohydrates, EPS responses, and pre/post-treatment contrasts.
- `total_carbohydrates.R` — Total carbohydrate analyses (colloidal + intracellular).
- `shear_stress_models.R` — Erosion type and sediment stability analyses from CSM output.
- `macrofauna_analysis.R` — Analyses macrofaunal responses and community composition (nMDS, PERMANOVA).

### Project & Metadata

- `stressors_field.Rproj` — R project file for consistent working environment.
- `README.md` — Summary of repository purpose and contents (this file).

---

## Reproducing Manuscript Results

This repository includes all R scripts and data required to reproduce figures 
and statistical results reported in the manuscript. Bayesian model results can be checked 
against the paper by re-running:

```r
source("NDVI and FvFm.R")            # Timepoint-specific mixed models
source("chlorophyll.R")             # Chlorophyll a
source("coll_carbohydrates.R")      # Colloidal EPS responses
source("total_carbohydrates.R")     # Total carbohydrates
source("shear_stress_models.R")     # CSM erosion threshold models
source("macrofauna_analysis.R")     # Macrofauna
source("environmental_models.R")    # Light/temperature effects on MPB
```

To ensure functionality:

- Open the `stressors_field.Rproj` project file
- Use `here::here()` for relative file paths
- Ensure packages are installed: `rstanarm`, `bayestestR`, `brms`, `vegan`, `lme4`, `readxl`, `tidyverse`, `bayesplot`, etc.

---
