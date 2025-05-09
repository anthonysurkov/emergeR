library(tidyverse)
library(here)
library(umap)
library(stats)

setwd("D:/ngs storage/Natalie/R255X E488QD/scripting/")
here::i_am("eda/edits_visualization.R")


# ============================================================================ #


input_file <- "R255X_E488QD_modeling.csv"

# Generate UMAP plots from the lower edit threshold to the upper edit threshold,
# generating plots once every step size.
# E.g. lower_edit_threshold <- 0.01, upper... <- 0.05, step <- 0.01
#      generates UMAP plots with all >0.01, >0.02, >0.03, >0.04, >0.05 editing
#      sequences.
#
lower_edit_threshold <- 0.05
upper_edit_threshold <- 0.30
step_size            <- 0.05
reads_threshold      <- 10

outdir <- "img/exploratory"


# ============================================================================ #


main <- function() {
  model_data   <- read_csv(here("data", input_file))
  model_data   <- model_data %>% 
    filter(edit >= lower_edit_threshold) %>%
    filter(n >= reads_threshold)
  
  umap_edits <- generate_umap_representation(model_data=model_data)
  
  # calculate axes
  x_limits <- range(umap_edits$X1)
  y_limits <- range(umap_edits$X2)
  x_pad    <- 0.05 * diff(x_limits)
  y_pad    <- 0.05 * diff(y_limits)
  x_limits <- c(x_limits[1] - x_pad, x_limits[2] + x_pad)
  y_limits <- c(y_limits[1] - y_pad, y_limits[2] + y_pad)
  
  for (i in seq(lower_edit_threshold, upper_edit_threshold, by=step_size)) {
    plot_path <- here(outdir,
      sprintf("local_seq_landscape_%5.2f.png", i)
    )
    ggsave(
      plot_path,
      plot=plot_umap_edits(
        umap_edits=umap_edits, 
        x_limits=x_limits,
        y_limits=y_limits,
        edit_cutoff=i),
      width=6, height=6, dpi=150
    )
  }   
  
}


generate_umap_representation <- function(model_data) {
  feature_data <- model_data %>% select(-n_other, -N10, -GGA, -GAA, -n, -edit) 
  
  # PCA to reduce feature dimensions
  pca_result <- prcomp(feature_data, center=TRUE, scale.=TRUE, rank.=10)
  pca_matrix <- pca_result$x

  # Generate UMAP reduction from PCA results
  message("UMAP reduction started...")
  feature_umap <- umap(pca_matrix)
  edit_data    <- as_tibble(feature_umap$layout, .name_repair="unique") %>%
    rename(X1 = 1, X2 = 2) %>%
    mutate(edit = model_data$edit)
  
  return(edit_data)
}

plot_umap_edits <- function(umap_edits, edit_cutoff, x_limits, y_limits) {
  edit_plot <- ggplot(umap_edits %>% filter(edit > edit_cutoff),
                      aes(x=X1, y=X2, color=edit)) +
    geom_point(size=1, alpha=0.5) +
    scale_color_viridis_c(limits=c(0,1)) +
    xlim(x_limits) +
    ylim(y_limits) +
    labs(
      title="Local Sequence Landscape by Editing Percentage",
      subtitle=paste0("Colored by editing; for editing > ", edit_cutoff),
      x="UMAP Dimension 1",
      y="UMAP Dimension 2"
    )
  
  return(edit_plot)
}


main()
