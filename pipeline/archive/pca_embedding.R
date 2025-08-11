library(tidyverse)
library(ggprism)
library(scales)
library(here)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/pca/cv_eb_editors_on_pca.R")

indir   <- "data"
outdir  <- "adsurkov"
infile  <- "R255X_E488QD_onehot_dummy_demeaned.csv"
outfile <- "R255X_E488QD_dbscan_cluster_of_pca.csv"


# ============================================================================ #
# Runs PCA on one-hot embedded (w/ dummy variable) sequence data to generate
# a 2D approximation of sequence space. Plots editing-capable sequences in PCA
# space and colors them according to maximum a posteriori estimate of editing
# capability.

model_data <- read_csv(here(indir, infile))

pca_mat <- model_data %>%
  select(
    -N10, -n, -k,
    -map, -upper_cred, -lower_cred,
    -mle, -upper_ci, -lower_ci,
    -alpha, -beta 
  ) %>%
  as.matrix()

pca_fit <- prcomp(pca_mat, center=FALSE, scale.=FALSE)
print(summary(pca_fit))

scores <- as_tibble(pca_fit$x) %>%
  mutate(N10 = model_data$N10) %>%
  left_join(model_data, by="N10") %>%
  filter(map > 0.10)

# fancy plot
plot <- ggplot(scores, aes(PC1, PC2, colour = map)) +
  geom_point(size = 1.0, alpha = 0.40) +
  scale_colour_gradientn(
    colours = c("#b0f3ff", "#003e4c", "#ff0000"),
    name    = "Editing Rate",
    limits  = c(0, 0.4),
    oob     = scales::squish
  ) +
  guides(
    colour = guide_legend(position="right")
  ) +
  theme_prism(
    palette = "winter_bright", 
    base_size = 14,
    base_family = "sans",
    base_fontface = "bold",
    base_line_size = 1,
    base_rect_size = 1,
    axis_text_angle = 0,
    border = FALSE
  ) +
  theme(
    plot.subtitle = element_text(size=12, margin=margin(t=-10)),
    legend.title = element_text(size=12),
    plot.title = element_text(size=16),
    legend.position = c(0.05, 0.95),
    legend.margin = margin(2, 2, 2, 2),
    legend.box.margin = margin(0, 0, 0, 0)
  ) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA of Editing-Capable Guides",
    subtitle = "(maximum a posteriori editing > 2%)"
  ) +
  scale_y_continuous(
    limits = c(-0.77, 1.2),
    guide = "prism_offset_minor"
  ) +
  scale_x_continuous(
    limits = c(-0.77, 1.2),
    guide = "prism_offset_minor"
  ) +
  coord_fixed()

print(plot)

