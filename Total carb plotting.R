##############################################################################
# Figure 7 - Mean total carbohydrate content of sediment 
# after stressor application, by TiO2 and Glyphosate.
##############################################################################

# 1) Load packages
library(dplyr)
library(plyr)
library(ggplot2)
library(here)

# 2) Read data
contact_core <- read.delim(here("contact_core_data.txt"))

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
  
  # Use plyr::rename to avoid conflicts with dplyr::rename:
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)
  ciMult   <- qt(conf.interval / 2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# 4) Convert columns to factors
contact_core$Stressor.status <- as.factor(contact_core$Stressor.status)
contact_core$Glyphosate      <- as.factor(contact_core$Glyphosate)
contact_core$TiO2            <- as.factor(contact_core$TiO2)

# 5) Subset to post-stressor data only
contact_post <- subset(contact_core, Stressor.status == "Post")

# 6) Summarise total carbohydrate by glyphosate and TiO2
sum_tot_post <- summarySE(
  data       = contact_post, 
  measurevar = "Total.carbs", 
  groupvars  = c("Glyphosate", "TiO2")
)

# 7) Convert factor levels to numeric-coded strings
sum_tot_post$Gly_values <- as.numeric(sum_tot_post$Glyphosate)
sum_tot_post$Gly_values[sum_tot_post$Gly_values == 1] <- 0
sum_tot_post$Gly_values[sum_tot_post$Gly_values == 2] <- 12000
sum_tot_post$Gly_values[sum_tot_post$Gly_values == 3] <- 20000
sum_tot_post$Gly_values <- as.factor(sum_tot_post$Gly_values)

sum_tot_post$TiO2_values <- as.numeric(sum_tot_post$TiO2)
sum_tot_post$TiO2_values[sum_tot_post$TiO2_values == 1] <- 0
sum_tot_post$TiO2_values[sum_tot_post$TiO2_values == 2] <- 1000
sum_tot_post$TiO2_values <- as.factor(sum_tot_post$TiO2_values)

# 8) Create the bar chart
fig7_tot_post <- ggplot(sum_tot_post, aes(
  x    = TiO2_values,
  y    = Total.carbs,
  fill = Gly_values
)) +
  geom_bar(
    position = position_dodge(width = 0.9),
    stat     = "summary",
    fun      = "mean",
    colour   = "black"
  ) +
  geom_errorbar(aes(
    ymin = Total.carbs - sd,
    ymax = Total.carbs + sd
  ),
  width = 0.3,
  position = position_dodge(width = 0.9)
  ) +
  coord_cartesian(ylim = c(0, 8000)) +  # Adjust as needed
  scale_fill_manual(
    values = c("0" = "#999999", "12000" = "#56B4E9", "20000" = "#0072B2"),
    name   = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  xlab(expression(paste("TiO"[2], " (", mu, "g L"^-1, ")"))) +
  ylab(expression(
    paste("Total carbohydrates (", mu, "g g"^-1, ")")
  )) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 20),
    axis.title  = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title= element_text(size = 20)
  )

# 9) Print or export
print(fig7_tot_post)
ggsave("Figure7.png", fig7_tot_post, width = 8, height = 6, dpi = 300)
