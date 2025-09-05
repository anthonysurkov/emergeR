#' Dimensionally reduces MCA embedding of EMERGe data using t-distributed
#' stochastic neighbor embedding (t-SNE).
#'
#' Motivated by insufficient separability of clusters on EMERGe MCA data alone.
#' Calculates hyperparameters that work on EMERGe data if left NULL.
#' Calculates a t-SNE fit and appends t-SNE coordinates to each n10 in X.
#' 
#' @param X_mca (tibble) EMERGe data with at least 'n10' and MCA columns that 
#'          start with "CS".
#' @param alpha (float) Exaggeration factor, handled internally by default.
#' @param h (float) Learning rate, handled internally by default.
#' @param seed (float) RNG seed, handled internally by default.
#' @param perplexity (float) effective number of nearest neighbors,
#'                   handled internally by default.
#' @param verbose (logical) whether to print Rtsne progress messages.
#' @param max_iter (integer) maximum number of iterations for
#'                 t-SNE optimization.
#' @return X (tibble) with 2D t-SNE coordinates attached as 'tsne1'
#'         and 'tsne2' columns.
#' @export
#' 
#' @importFrom Rtsne Rtsne
#' @importFrom dplyr select mutate
#' @importFrom tidyselect starts_with
#' @importFrom magrittr %>%
append_tsne <- function(
    X_mca, alpha = NULL, h = NULL, 
    seed = NULL, perplexity = 30,
    verbose = FALSE, max_iter = 30
) {
  mca <- X_mca %>% select(starts_with("CS"))
  
  n <- nrow(mca)
  if (alpha == NULL) alpha <- n / 10
  if (h == NULL) h <- n / alpha
  if (seed == NULL) set.seed(1)
  
  tsne_fit <- Rtsne(
    as.matrix(mca),
    exaggeration_factor = alpha,
    eta                 = h,
    dims                = 2,
    perplexity          = perplexity,
    verbose             = verbose,
    max_iter            = max_iter,
    PCA                 = FALSE
  )
  
  X_tsne <- X_mca %>%
    mutate(
      tsne1 = tsne_fit$Y[, 1],
      tsne2 = tsne_fit$Y[, 2]
    )
  
  return(X_tsne)
}

#' Plots t-SNE coordinates colored by posterior mean editing.
#' 
#' @param X_tsne (tibble) EMERGe data with at minimum 'tsne1', 'tsne2', and
#' #        'postmean' columns.
#' @return (plot) ggprism-formatted plot.
#' @export
#' 
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradientn guides
#' @importFrom ggplot2 guide_legend theme labs scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous coord_fixed element_text margin
#' @importFrom ggprism theme_prism
#' @importFrom scales squish
plot_tsne_edits <- function(X_tsne) {
  stopifnot(all(c("n10", "tsne1", "tsne2")))
  
  p <- ggplot(X_tsne, aes(tsne1, tsne2, colour = postmean)) +
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
      x = "TSNE1",
      y = "TSNE2",
      title = "t-SNE of Embedded Editing-Capable Guides",
      subtitle = "99.9th percentile Maximum A Posteriori Editors"
    ) +
    scale_y_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    scale_x_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    coord_fixed()
  return(p)
}

#' Plots t-SNE coordinates colored by cluster assignment.
#' 
#' @param X_tsne (tibble) EMERGe data with at minimum 'tsne1', 'tsne2', and
#' #        'cluster' columns.
#' @return (plot) ggprism-formatted plot.
#' @export
#' 
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradientn guides
#' @importFrom ggplot2 guide_legend theme labs scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous coord_fixed element_text margin
#' @importFrom ggprism theme_prism
#' @importFrom scales squish
plot_tsne_clusters <- function(X_tsne) {
  stopifnot(all(c("n10", "tsne1", "tsne2")))
  
  p <- ggplot(X, aes(tsne1, tsne2, colour = cluster)) +
    geom_point(size = 1.0, alpha = 0.40) +
    scale_colour_discrete(
      name = "Cluster ID",
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
      x = "TSNE1",
      y = "TSNE2",
      title = "t-SNE of Embedded Editing-Capable Guides",
      subtitle = "99.9th percentile Maximum A Posteriori Editors"
    ) +
    scale_y_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    scale_x_continuous(
      limits = c(-21, 21),
      guide = "prism_offset_minor"
    ) +
    coord_fixed() 
  return(p)
}