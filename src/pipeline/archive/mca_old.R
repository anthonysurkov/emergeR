library(tidyverse)
library(ggprism)
library(scales)
library(here)
library(ade4)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/mca/mca+old.R")

indir   <- "data"
outdir  <- "adsurkov"
infile  <- "R255X_E488QD_clean.csv"
outfile <- "R255X_E488QD_mca.csv"


# ============================================================================ #
# Runs MCA on one-hot embedded (w/ dummy variable) sequence data to generate
# a 2D approximation of sequence space. Plots editing-capable sequences in PCA
# space and colors them according to maximum a posteriori estimate of editing
# capability.

L <- 10 # sequence length, fixed

# build factor table p1..pL from N10
Xfac <- model_data$N10 %>%
  strsplit("") %>% do.call(rbind, .) %>% as.data.frame()
colnames(Xfac) <- paste0("p", 1:L)
Xfac[] <- lapply(Xfac, function(x) factor(x, levels = c("A","C","G","T")))

# MCA with ade4
acm <- dudi.acm(Xfac, nf = 50, scannf = FALSE)

scores <- as_tibble(acm$li) %>%
  rename(MC1 = 1, MC2 = 2) %>%
  mutate(N10 = model_data$N10) %>%
  left_join(model_data, by = "N10") %>%
  filter(map > 0.02)

print(acm$eig)

# print 2d plot
plot <- ggplot(scores, aes(MC1, MC2, colour = map)) +
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
    x = "MC1",
    y = "MC2",
    title = "MCA of Editing-Capable Guides",
    subtitle = "(maximum a posteriori editing > 5%)"
  ) +
  scale_y_continuous(
    limits = c(-0.87, 0.87),
    guide = "prism_offset_minor"
  ) +
  scale_x_continuous(
    limits = c(-0.87, 0.87),
    guide = "prism_offset_minor"
  ) +
  coord_fixed()

print(plot)


