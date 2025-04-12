##############################################################################
# Figure 8 - Mean sediment shear stress of critical erosion (Type I) 
# following stressor application, by TiO2 and Glyphosate.
##############################################################################

# 1) Load packages
library(dplyr)
library(plyr)
library(ggplot2)
library(here)

# 2) Read data
CSM.dat <- read.delim(here("CSM formatted.txt"))

# 3) Define summarySE
summarySE <- function(data = NULL, measurevar, groupvars = NULL,
                      na.rm = FALSE, conf.interval = 0.95, .drop = TRUE) {
  library(plyr)
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
  ciMult   <- qt(conf.interval / 2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

# 4) Convert key columns to factors
CSM.dat$Treatment       <- as.factor(CSM.dat$Treatment)
CSM.dat$Stressor.status <- as.factor(CSM.dat$Stressor.status)
CSM.dat$Erosion.type    <- as.factor(CSM.dat$Erosion.type)
CSM.dat$Glyphosate      <- as.factor(CSM.dat$Glyphosate)
CSM.dat$TiO2            <- as.factor(CSM.dat$TiO2)

# 5) Subset to post-stress data with Type I erosion only
CSM.post.I <- subset(CSM.dat, Stressor.status == "Post" & Erosion.type == "1")

# 6) Summarise shear stress by glyphosate and TiO2
sum_ss_post <- summarySE(
  data       = CSM.post.I,
  measurevar = "Shear.stress",
  groupvars  = c("Glyphosate", "TiO2")
)

# 7) Convert factor levels to numeric-coded strings
sum_ss_post$Gly_values <- as.numeric(sum_ss_post$Glyphosate)
sum_ss_post$Gly_values[sum_ss_post$Gly_values == 1] <- 0
sum_ss_post$Gly_values[sum_ss_post$Gly_values == 2] <- 12000
sum_ss_post$Gly_values[sum_ss_post$Gly_values == 3] <- 20000
sum_ss_post$Gly_values <- as.factor(sum_ss_post$Gly_values)

sum_ss_post$TiO2_values <- as.numeric(sum_ss_post$TiO2)
sum_ss_post$TiO2_values[sum_ss_post$TiO2_values == 1] <- 0
sum_ss_post$TiO2_values[sum_ss_post$TiO2_values == 2] <- 1000
sum_ss_post$TiO2_values <- as.factor(sum_ss_post$TiO2_values)

# 8) Create the bar chart
fig8_shear_post <- ggplot(sum_ss_post, aes(
  x    = TiO2_values,
  y    = Shear.stress,
  fill = Gly_values
)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(
    ymin = Shear.stress - sd,
    ymax = Shear.stress + sd
  ),
  width = 0.3,
  position = position_dodge(width = 0.9)
  ) +
  coord_cartesian(ylim = c(0, 0.12)) +  # Adjust as needed
  scale_fill_manual(
    values = c("0" = "#999999", "12000" = "#56B4E9", "20000" = "#0072B2"),
    name   = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  xlab(expression(paste("TiO"[2], " (", mu, "g L"^-1, ")"))) +
  ylab(expression(paste("Shear stress (N m"^-2, ")"))) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 20),
    axis.title  = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title= element_text(size = 20)
  )

# 9) Print or save
print(fig8_shear_post)
ggsave("Figure8.png", fig8_shear_post, width = 8, height = 6, dpi = 300)
