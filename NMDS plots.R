##############################################################################
# Figure 9 - NMDS Ordination of Macrofauna Community (Bray-Curtis)
##############################################################################

# 1) Load necessary packages
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)  # to combine two plots side by side

# 2) Read data
# (Adjust file path if needed, or use here::here("Macrofauna.txt") if you prefer)
Macrofauna <- read.delim("Macrofauna.txt")

# 3) Data preparation
Macrofauna[is.na(Macrofauna)] <- 0
Macrofauna$Treatment   <- as.factor(Macrofauna$Treatment)
Macrofauna$Glyphosate  <- as.factor(Macrofauna$Glyphosate)
Macrofauna$TiO2        <- as.factor(Macrofauna$TiO2)

# Extract only the abundance columns (adjust columns 7:17 if your columns differ)
Macrofauna_2 <- Macrofauna[, 7:17]

# 4) Run NMDS with k=2 using Bray-Curtis
# You mentioned stress ~ 0.12 for k=2
NMDS2 <- metaMDS(comm = Macrofauna_2, distance = "bray",
                 k = 2, trymax = 100, trace = FALSE, autotransform = FALSE)

# 5) Extract NMDS site scores
nmds_scores <- scores(NMDS2, display = "sites") %>% as.data.frame()
nmds_scores$Glyphosate <- Macrofauna$Glyphosate
nmds_scores$TiO2       <- Macrofauna$TiO2
nmds_scores$Treatment  <- Macrofauna$Treatment

# 6) Convert factor levels to numeric-coded strings 
#    so that e.g. glyphosate: 1->0, 2->12000, 3->20000
nmds_scores <- nmds_scores %>%
  mutate(
    Gly_values = as.numeric(Glyphosate),
    TiO2_values = as.numeric(TiO2)
  )

# For glyphosate:
nmds_scores$Gly_values[nmds_scores$Gly_values == 1] <- 0
nmds_scores$Gly_values[nmds_scores$Gly_values == 2] <- 12000
nmds_scores$Gly_values[nmds_scores$Gly_values == 3] <- 20000
nmds_scores$Gly_values <- factor(nmds_scores$Gly_values)

# For TiO2:
nmds_scores$TiO2_values[nmds_scores$TiO2_values == 1] <- 0
nmds_scores$TiO2_values[nmds_scores$TiO2_values == 2] <- 1000
nmds_scores$TiO2_values <- factor(nmds_scores$TiO2_values)

# 7) Function to find hull points for each group
#    We can group by glyphosate or TiO2 to draw polygons around replicates
find_hull <- function(df) {
  df[chull(df$NMDS1, df$NMDS2), ]
}

# 8) Prepare hull data for glyphosate grouping
hull_gly <- nmds_scores %>%
  group_by(Gly_values) %>%
  do(find_hull(.))

# 9) Plot (a): NMDS by glyphosate
fig9a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = Gly_values)) +
  geom_point(size = 3) +
  geom_polygon(
    data = hull_gly,
    aes(x = NMDS1, y = NMDS2, fill = Gly_values),  # fill for the polygon
    alpha = 0.2,
    colour = NA,   # or "black" if you want outlines
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("0" = "#999999", "12000" = "#56B4E9", "20000" = "#0072B2"),
    name = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  scale_fill_manual(
    values = c("0" = "#999999", "12000" = "#56B4E9", "20000" = "#0072B2"),
    name = expression(paste("Glyphosate (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "(a) NMDS by Glyphosate"
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# 10) Prepare hull data for TiO2 grouping
hull_tio2 <- nmds_scores %>%
  group_by(TiO2_values) %>%
  do(find_hull(.))

# 11) Plot (b): NMDS by TiO2
fig9b <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, colour = TiO2_values)) +
  geom_point(size = 3) +
  geom_polygon(
    data = hull_tio2,
    aes(x = NMDS1, y = NMDS2, fill = TiO2_values),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("0" = "#999999", "1000" = "#994F00"),
    name = expression(paste("TiO"[2], " (", mu, "g L"^-1, ")"))
  ) +
  scale_fill_manual(
    values = c("0" = "#999999", "1000" = "#994F00"),
    name = expression(paste("TiO"[2], " (", mu, "g L"^-1, ")"))
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "(b) NMDS by TiO2"
  ) +
  theme_classic() +
  theme(
    axis.text   = element_text(size = 12),
    axis.title  = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title= element_text(size = 14)
  )

# 12) Combine side by side
Figure9 <- fig9a + fig9b

# 13) Print or save
print(Figure9)
ggsave("Figure9.png", Figure9, width = 12, height = 5, dpi = 300)
