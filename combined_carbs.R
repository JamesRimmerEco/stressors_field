##############################################################################
# Combined Figure: Colloidal and Total Carbohydrates by Treatment
##############################################################################

# 1) Load packages
library(dplyr)
library(plyr)
library(ggplot2)
library(here)
library(patchwork)

# 2) Read data
contact_core <- read.delim(here("contact_core_data.txt"))

# 3) Define summarySE
summarySE <- function(data = NULL, measurevar, groupvars = NULL,
                      na.rm = FALSE, conf.interval = 0.95, .drop = TRUE) {
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x)) else length(x)
  }
  datac <- ddply(
    data, groupvars, .drop = .drop,
    .fun = function(xx, col) {
      c(N    = length2(xx[[col]], na.rm = na.rm),
        mean = mean(xx[[col]], na.rm = na.rm),
        sd   = sd(xx[[col]], na.rm = na.rm))
    },
    measurevar
  )
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)
  ciMult <- qt(conf.interval / 2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# 4) Convert columns to factors
contact_core$Stressor.status <- as.factor(contact_core$Stressor.status)
contact_core$Glyphosate      <- as.factor(contact_core$Glyphosate)
contact_core$TiO2            <- as.factor(contact_core$TiO2)

# 5) Subset to post-stressor data only
contact_post <- subset(contact_core, Stressor.status == "Post")

# 6a) Colloidal carbohydrate summary
sum_coll_post <- summarySE(contact_post, "Colloidal.carbs", c("Glyphosate", "TiO2"))
sum_coll_post$Gly_values <- factor(recode(as.numeric(sum_coll_post$Glyphosate), `1`=0, `2`=12000, `3`=20000))
sum_coll_post$TiO2_values <- factor(recode(as.numeric(sum_coll_post$TiO2), `1`=0, `2`=1000))

# 6b) Total carbohydrate summary
sum_tot_post <- summarySE(contact_post, "Total.carbs", c("Glyphosate", "TiO2"))
sum_tot_post$Gly_values <- factor(recode(as.numeric(sum_tot_post$Glyphosate), `1`=0, `2`=12000, `3`=20000))
sum_tot_post$TiO2_values <- factor(recode(as.numeric(sum_tot_post$TiO2), `1`=0, `2`=1000))

# 7) Plot colloidal carbs
fig_coll <- ggplot(sum_coll_post, aes(x = TiO2_values, y = Colloidal.carbs, fill = Gly_values)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(ymin = Colloidal.carbs - sd, ymax = Colloidal.carbs + sd),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  coord_cartesian(ylim = c(0, 360)) +
  scale_fill_manual(
    values = c("0"="#999999", "12000"="#56B4E9", "20000"="#0072B2"),
    name   = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x = expression(paste("TiO"[2], " (", mu, "g L"^-1, ")")),
    y = expression(paste("Colloidal carbohydrates (", mu, "g g"^-1, ")"))
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# 8) Plot total carbs
fig_tot <- ggplot(sum_tot_post, aes(x = TiO2_values, y = Total.carbs, fill = Gly_values)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(ymin = Total.carbs - sd, ymax = Total.carbs + sd),
                width = 0.3,
                position = position_dodge(width = 0.9)) +
  coord_cartesian(ylim = c(0, 8000)) +
  scale_fill_manual(
    values = c("0"="#999999", "12000"="#56B4E9", "20000"="#0072B2"),
    name   = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x = expression(paste("TiO"[2], " (", mu, "g L"^-1, ")")),
    y = expression(paste("Total carbohydrates (", mu, "g g"^-1, ")"))
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# 9) Combine and export
Figure5 <- (fig_coll / fig_tot) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 22, face = "bold")
  )

# Save
ggsave("Figure5.png", Figure5, width = 8, height = 10, dpi = 300)
