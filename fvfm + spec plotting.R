###############################################################################
# Reading data, formatting, and creating bar charts for FvFm 
# and NDVI across time, grouped by glyphosate and TiO2. 
###############################################################################

# 1) Load libraries
library(readxl)   # For reading XLSX
library(dplyr)    # Data wrangling
library(plyr)     # 'ddply' for summarySE
library(ggplot2)  # Plotting
library(patchwork)

# 2) Read data
my_data <- read_excel("Processed data for R.xlsx")
names(my_data)[names(my_data) == "Experiment ID"] <- "Experiment.ID"
my_data$Experiment.ID <- as.factor(my_data$Experiment.ID)

# 3) Rename it for convenience
field_dat <- my_data

# 4) Convert relevant columns to factor, remove negative NDVI, etc.
field_dat$Treatment      <- as.factor(field_dat$Treatment)
field_dat$Glyphosate     <- as.factor(field_dat$Glyphosate)
field_dat$TiO2           <- as.factor(field_dat$TiO2)
field_dat$Experiment.ID  <- as.factor(field_dat$Experiment.ID)

# Remove negative NDVI entries (replace with NA)
field_dat$NDVI[field_dat$NDVI < 0] <- NA

# Specify whether stress was applied (time > 1 = "post-stress"):
field_dat$stress_applied <- as.numeric(field_dat$Time > 1)

# 5) Define the summarySE function (from Cookbook for R)
summarySE <- function(data = NULL, measurevar, groupvars = NULL, 
                      na.rm = FALSE, conf.interval = 0.95, .drop = TRUE) {
  # This function calculates count (N), mean, sd, se, and CI for each group
  library(plyr)
  
  # Helper for counting non-NA
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) else length(x)
  }
  
  # Summarise the data by groups
  datac <- ddply(
    data, groupvars, .drop = .drop,
    .fun = function(xx, col) {
      c(N    = length2(xx[[col]], na.rm = na.rm),
        mean = mean(xx[[col]], na.rm = na.rm),
        sd   = sd(xx[[col]], na.rm = na.rm))
    },
    measurevar
  )
  
  # Rename "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  
  # Calculate standard error
  datac$se <- datac$sd / sqrt(datac$N)
  
  # Calculate confidence intervals
  ciMult <- qt(conf.interval / 2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###############################################################################
# 6) Prepare data frames for plotting FvFm and NDVI with Time, Glyphosate, TiO2
###############################################################################

# (a) FvFm subset
field_dat_fvfm <- field_dat %>%
  select(Time, Glyphosate, TiO2, FvFm) %>%
  filter(!is.na(FvFm))  # remove rows lacking FvFm

# (b) NDVI subset
field_dat_ndvi <- field_dat %>%
  select(Time, Glyphosate, TiO2, NDVI) %>%
  filter(!is.na(NDVI))  # remove rows lacking NDVI

###############################################################################
# 7) Create bar charts for FvFm (Figure 3)
###############################################################################

# 7.1 Summaries of FvFm by (Time, Glyphosate)
fvfm_sum_gly <- summarySE(
  data       = field_dat_fvfm,
  measurevar = "FvFm",
  groupvars  = c("Time", "Glyphosate")
)

# 7.2 Summaries of FvFm by (Time, TiO2)
fvfm_sum_tio2 <- summarySE(
  data       = field_dat_fvfm,
  measurevar = "FvFm",
  groupvars  = c("Time", "TiO2")
)

fvfm_sum_gly$Gly_values <- as.numeric(as.character(fvfm_sum_gly$Glyphosate))
fvfm_sum_gly$Gly_values[fvfm_sum_gly$Gly_values == 1] <- 12000
fvfm_sum_gly$Gly_values[fvfm_sum_gly$Gly_values == 2] <- 20000
# 0 stays 0
fvfm_sum_gly$Gly_values <- as.factor(fvfm_sum_gly$Gly_values)

fvfm_sum_tio2$TiO2_values <- as.numeric(as.character(fvfm_sum_tio2$TiO2))
fvfm_sum_tio2$TiO2_values[fvfm_sum_tio2$TiO2_values == 1] <- 1000

# 0 stays 0
fvfm_sum_tio2$TiO2_values <- as.factor(fvfm_sum_tio2$TiO2_values)

# 7.3 Plot Figure 3(a): FvFm vs Time, grouped by Glyphosate
fig3a_fvfm_gly <- ggplot(fvfm_sum_gly, aes(
  x    = factor(Time),
  y    = FvFm,
  fill = Gly_values
)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  # ± SD error bars
  geom_errorbar(aes(ymin = FvFm - sd, ymax = FvFm + sd),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    # Grey for 0, Light Blue for 12000, Dark Blue for 20000
    values = c("0" = "#999999", "12000" = "#56B4E9", "20000" = "#0072B2"),
    name   = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x     = "Time",
    y     = expression(paste("F"[v], "/", "F"[m]))
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# 7.4 Plot Figure 3(b): FvFm vs Time, grouped by TiO2
fig3b_fvfm_tio2 <- ggplot(fvfm_sum_tio2, aes(
  x    = factor(Time),
  y    = FvFm,
  fill = TiO2_values
)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(ymin = FvFm - sd, ymax = FvFm + sd),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    # Grey for 0, Brown for 1000
    values = c("0" = "#999999", "1000" = "#994F00"),
    name   = expression(paste("TiO"[2], " (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x     = "Time",
    y     = expression(paste("F"[v], "/", "F"[m]))
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# Print or view these:
fig3a_fvfm_gly
fig3b_fvfm_tio2

###############################################################################
# 8) Create bar charts for NDVI (Figure 4)
###############################################################################

# 8.1 Summaries by (Time, Glyphosate)
ndvi_sum_gly <- summarySE(
  data       = field_dat_ndvi,
  measurevar = "NDVI",
  groupvars  = c("Time", "Glyphosate")
)

# 8.2 Summaries by (Time, TiO2)
ndvi_sum_tio2 <- summarySE(
  data       = field_dat_ndvi,
  measurevar = "NDVI",
  groupvars  = c("Time", "TiO2")
)

# Convert factor levels for glyphosate, TiO2 if needed:
ndvi_sum_gly$Gly_values <- as.numeric(as.character(ndvi_sum_gly$Glyphosate))
ndvi_sum_gly$Gly_values[ndvi_sum_gly$Gly_values == 1] <- 12000
ndvi_sum_gly$Gly_values[ndvi_sum_gly$Gly_values == 2] <- 20000
ndvi_sum_gly$Gly_values <- as.factor(ndvi_sum_gly$Gly_values)

ndvi_sum_tio2$TiO2_values <- as.numeric(as.character(ndvi_sum_tio2$TiO2))
ndvi_sum_tio2$TiO2_values[ndvi_sum_tio2$TiO2_values == 1] <- 1000
ndvi_sum_tio2$TiO2_values <- as.factor(ndvi_sum_tio2$TiO2_values)

# 8.3 Figure 4(a): NDVI vs Time, grouped by Glyphosate
fig4a_ndvi_gly <- ggplot(ndvi_sum_gly, aes(
  x    = factor(Time),
  y    = NDVI,
  fill = Gly_values
)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(ymin = NDVI - sd, ymax = NDVI + sd),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("0" = "#999999", "12000" = "#56B4E9", "20000" = "#0072B2"),
    name   = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x     = "Time",
    y     = "NDVI"
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# 8.4 Figure 4(b): NDVI vs Time, grouped by TiO2
fig4b_ndvi_tio2 <- ggplot(ndvi_sum_tio2, aes(
  x    = factor(Time),
  y    = NDVI,
  fill = TiO2_values
)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(ymin = NDVI - sd, ymax = NDVI + sd),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  scale_fill_manual(
    values = c("0" = "#999999", "1000" = "#994F00"),
    name   = expression(paste("TiO"[2], " (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x     = "Time",
    y     = "NDVI"
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

fig4a_ndvi_gly
fig4b_ndvi_tio2

Figure2 <- (fig3a_fvfm_gly / fig3b_fvfm_tio2) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 22, face = "bold")
  )

# Save the stacked figure
ggsave("Figure2.png", Figure2, width = 8, height = 10, dpi = 300)

# Similarly for Figure 4
Figure3 <- (fig4a_ndvi_gly / fig4b_ndvi_tio2) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 22, face = "bold")
  )

ggsave("Figure3.png", Figure3, width = 8, height = 10, dpi = 300)
