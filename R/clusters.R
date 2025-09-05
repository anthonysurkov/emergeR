#' Clusters the t-SNE embedding of EMERGe data.
#' 
#' Runs Hierarchical Density-Based Clustering with Applications to Noise
#' (HDBSCAN) to cluster the t-SNE embedding of EMERGe data.
#' 
#' @param X (tibble) EMERGe data with at minimum 'n10', 'tsne1', and 'tsne2'
#'          columns.
#' @param minPts (integer) Minimum number of points necessary to define a
#'        cluster.
#' @return Provided X (tibble) with cluster assignment appended as a
#'         'cluster' column.
#' @export
#' 
#' @importFrom dbscan hdbscan
#' @importFrom dplyr select mutate
#' @importFrom magrittr %>%
append_clusters <- function(X, minPts = 10) {
  stopifnot(all(c("n10", "tsne1", "tsne2")))
  
  X_tsne <- X %>%
    select(tsne1, tsne2)
  
  hd <- hdbscan(X_tsne, minPts = minPts)
  
  X_clustered <- X %>%
    mutate(
      cluster_raw = hd$cluster,
      cluster = factor(
        ifelse(hd$cluster == 0, "Noise", as.character(hd$cluster)),
        levels = c("Noise", sort(unique(hd$cluster[hd$cluster != 0])))
      )
    ) %>%
    select(-cluster_raw)
  
  return(X_clustered)
}

#' Plots a sequence logo per each identified EMERGe cluster.
#' 
#' @param X_clustered (tibble) EMERGe data with at minimum 'n10' and 'cluster'
#'                    columns.
#' @return (list of plots) A list of seq logos, sized according to the number
#'         of clusters.
#' @export
#' 
#' @importFrom dplyr group_by group_split pull
#' @importFrom purrr map_int imap
#' @importFrom ggseqlogo ggseqlogo
#' @importFrom ggplot2 ggtitle guides theme element_text
#' @importFrom ggprism theme_prism
plot_motif_logos <- function(X_clustered) {
  stopifnot(all(c("n10", "cluster")))
  
  seq_grouped <- X_clustered %>%
    group_by(cluster) %>%
    group_split(.keep = FALSE) %>%
    { .[ order(map_int(., nrow), decreasing = TRUE) ]}
  
  logo_list <- imap(seq_grouped, ~ {
    ggseqlogo(
      pull(.x, n10),
      method = "bits",
      seq_type = "dna",
      col_scheme = "nucleotide"
    ) +
      ggtitle(
        paste0("99.9th Percentile Editors (MAP), Cluster ", .y, 
               " (n = ", nrow(.x), ")"
        )
      ) +
      guides(colour = "none") +
      theme_prism(
        palette         = "winter_bright",
        base_size       = 14,
        base_family     = "sans",
        base_fontface   = "bold",
        base_line_size  = 1,
        base_rect_size  = 1,
        axis_text_angle = 0,
        border          = FALSE
      ) +
      theme(
        plot.title   = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_text(face = "bold"),
        axis.text.y  = element_text(face = "bold")
      )
  }) 
  return(logo_list)
}