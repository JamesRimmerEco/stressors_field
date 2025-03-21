# Environmental drivers of Fv/Fm and NDVI – logger-derived temperature and light data
# This script processes surface and subsurface logger data to generate daytime and nighttime means for light and temperature.
# These environmental variables are linked to repeated measures of Fv/Fm and NDVI, and analysed using mixed models.

# Load packages
library(here)
library(readxl)
library(lubridate)
library(lme4)
library(jtools)
library(performance)
library(nlme)

# Read in logger data and biofilm responses
Logger_1 <- read_excel(here("Logger_1.xlsx"), col_types = c("numeric", "date", "date", "skip", "numeric", "numeric"))
Logger_7 <- read_excel(here("Logger_7.xlsx"), col_types = c("numeric", "date", "date", "skip", "numeric", "numeric"))
Moleman_1 <- read_excel(here("Moleman_1.xlsx"), col_types = c("numeric", "date", "date", "skip", "numeric", "numeric"))
Moleman_2 <- read_excel(here("Moleman_2.xlsx"), col_types = c("numeric", "date", "date", "skip", "numeric", "numeric"))
field_dat <- read_excel(here("Processed data for R.xlsx"))

# Subset logger data to match experiment duration
Logger_1 <- subset(Logger_1, Date_time < "2021-08-24 07:35:00")
Logger_7 <- subset(Logger_7, Date_time < "2021-08-24 07:35:00")
Moleman_1 <- subset(Moleman_1, Date_time < "2021-08-24 07:35:00")
Moleman_2 <- subset(Moleman_2, Date_time < "2021-08-24 07:35:00")

# Average surface logger readings (temperature and light)
average_surface <- Logger_1[1:3]
average_surface$Temp <- rowMeans(cbind(Logger_1[4], Logger_7[4]))
average_surface$Lux <- rowMeans(cbind(Logger_1[5], Logger_7[5]))

# Average subsurface temperature
average_sub <- Moleman_2[1:3]
average_sub$Temp <- rowMeans(cbind(Moleman_1[4], Moleman_2[4]))

# Helper function: calculate daily means across specified intervals
daily_means <- function(intervals, data, cols = 4:5) {
  do.call(rbind, lapply(intervals, function(int) {
    colMeans(data[data$Date_time %within% int, cols], na.rm = TRUE)
  }))
}

# Define day and night intervals (12–23 August)
day_intervals <- lapply(12:23, function(d) interval(
  ymd_hm(paste0("2021-08-", d, " 05:35")),
  ymd_hm(paste0("2021-08-", d, " 20:55"))
))
night_intervals <- lapply(12:23, function(d) interval(
  ymd_hm(paste0("2021-08-", d - 1, " 20:55")),
  ymd_hm(paste0("2021-08-", d, " 05:35"))
))

# Surface day and night means
day_means <- as.data.frame(daily_means(day_intervals, average_surface))
night_means <- as.data.frame(daily_means(night_intervals, average_surface))

# Subsurface day and night means
day_means_s <- as.data.frame(daily_means(day_intervals, average_sub, cols = 4))
night_means_s <- as.data.frame(daily_means(night_intervals, average_sub, cols = 4))

# Attach environmental data to repeated measures dataset
field_dat$day_temp <- rep(day_means[c(1,3,5,7,10,12,14), 1], each = 42)
field_dat$lux <- rep(day_means[c(1,3,5,7,10,12,14), 2], each = 42)
field_dat$night_temp <- rep(night_means[c(1,3,5,7,10,12,14), 1], each = 42)
field_dat$sub_day_temp <- rep(day_means_s[c(1,3,5,7,10,12,14), 1], each = 42)
field_dat$sub_night_temp <- rep(night_means_s[c(1,3,5,7,10,12,14), 1], each = 42)

# Identify sampling period (AM/PM)
field_dat$Period <- factor(rep(c("AM", "AM", "AM", "PM", "PM", "PM", "AM"), each = 42))

# Transformations for modelling
field_dat$fvfm_pos <- field_dat$FvFm + 0.001
field_dat$log_day_temp <- log(field_dat$day_temp)
field_dat$lux_scale <- scale(field_dat$lux)
field_dat$scaled_day_temp_sub <- scale(field_dat$sub_day_temp)

# Models: Fv/Fm with different predictors
m_fv_1 <- glmer(fvfm_pos ~ day_temp + Period + (1|Time), family = inverse.gaussian(link = "log"), data = field_dat)
m_fv_2 <- glmer(fvfm_pos ~ lux_scale + Period + (1|Time), family = inverse.gaussian(link = "log"), data = field_dat)
m_fv_3 <- glmer(fvfm_pos ~ night_temp + Period + (1|Time), family = inverse.gaussian(link = "log"), data = field_dat)
m_fv_4 <- glmer(fvfm_pos ~ scaled_day_temp_sub + Period + (1|Time), family = inverse.gaussian(link = "log"), data = field_dat)
m_fv_5 <- glmer(fvfm_pos ~ sub_night_temp + Period + (1|Time), family = inverse.gaussian(link = "log"), data = field_dat)

# Model comparisons for Fv/Fm
AIC(m_fv_1, m_fv_2, m_fv_3, m_fv_4, m_fv_5)
anova(m_fv_1, m_fv_2, m_fv_3, m_fv_4, m_fv_5)

# Models: NDVI with different predictors
field_dat$NDVI_sqrt <- sqrt(field_dat$NDVI)

m_ndvi_1 <- lmer(NDVI_sqrt ~ day_temp + Period + (1|Time), data = field_dat)
m_ndvi_2 <- lmer(NDVI_sqrt ~ lux_scale + Period + (1|Time), data = field_dat)
m_ndvi_3 <- lmer(NDVI_sqrt ~ night_temp + Period + (1|Time), data = field_dat)
m_ndvi_4 <- lmer(NDVI_sqrt ~ sub_day_temp + Period + (1|Time), data = field_dat)
m_ndvi_5 <- lmer(NDVI_sqrt ~ sub_night_temp + Period + (1|Time), data = field_dat)

# Model comparisons for NDVI
AIC(m_ndvi_1, m_ndvi_2, m_ndvi_3, m_ndvi_4, m_ndvi_5)
anova(m_ndvi_1, m_ndvi_2, m_ndvi_3, m_ndvi_4, m_ndvi_5)
