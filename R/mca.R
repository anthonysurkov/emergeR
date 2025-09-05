#' Transforms cleaned EMERGe data into a richer representation using
#' multiple correspondence analysis (MCA).
#' 
#' Similar to PCA, MCA identifies the main axes of variation across
#' categorical variables and projects the data along them.
#' 
#' @param X (tibble) EMERGe dataset with 'n10' and 'postmean' columns.
#'          The 'postmean' column may be obtained from emergeR::stats().
#' @param nf (integer) default 50, maximum number of factors kept from MCA
#'           (dudi.acm)
#' @param top_quantile (float) default 0.999, defines the percentile by which
#'                     sequences are considered to be top performers.
#'                     For a dataset with ~900,000 n10 sequences, the default
#'                     value selects ~900 top sequences for further analysis.
#' @return A list with two components:
#'   \describe{
#'     \item{edit_top}{A tibble containing the subset of EMERGe data
#'       corresponding to sequences in the top quantile, joined with
#'       their MCA coordinates.}
#'     \item{eigens}{A numeric vector of eigenvalues from the MCA,
#'       representing the variance explained by each axis.}
#'   }
#' @export
#' 
#' @importFrom ade4 dudi.acm
#' @importFrom stats quantile
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
append_mca <- function(X, nf = 50, top_quantile = 0.999) {
  Xfac <- factor_sequences(X = X)
  mca  <- dudi.acm(df = Xfac, nf = 50, scannf = FALSE)
  
  top_cutoff <- quantile(X$postmean, top_quantile)
  top_seqs   <- X$n10[X$postmean >= top_cutoff]
  X_top      <- X %>%
    filter(n10 %in% top_seqs)
  
  edit_top  <- get_MCs(X = X_top, mca = mca)
  
  list(edit_top = edit_top, eigens = mca$eig)
}

#' Factorizes n10 sequences into position-wise categorical variables.
#' 
#' Splits each sequence in the 'n10' column into its nucleotide characters,
#' aligns them by position, and returns a data frame with one column per
#' position. Each column is a factor with levels A, G, C, T.
#' 
#' @param X (data frame) with an 'n10' column of DNA sequences.
#' @return data frame with factored nucleotides from n10 sequences.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr %>%
factor_sequences <- function(X) {
  Xfac <- X$n10 %>%
    strsplit("") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(Xfac) <- paste0("p", 1:L)
  Xfac[] <- lapply(Xfac, function(x) factor(x, levels = c("A", "G", "C", "T")))
}

#' Extract MCA scores for EMERGe guides.
#' 
#' Joins MCA coordinates with the underlying EMERGe dataset, keyed by
#' N10 sequence.
#' 
#' @param X (data frame) with an 'n10' column of DNA sequences.
#' @param mca (MCA result) an object returned by dudi.acm or similar,
#'            containing factor scores in the '$li' comoponent.
#' @return (tibble) X data frame with mca data left joined by n10.
#' 
#' @keywords internal
#' @noRd
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate left_join filter select
get_MCs <- function(X, mca) {
  scores <- as_tibble(mca$li) %>%
    mutate(N10 = X$N10) %>%
    left_join(X, by = "N10") %>%
    select(
      n10, n, k,
      postmean, lower_cred, upper_cred,
      mle, lower_ci, upper_ci,
      alpha, beta,
      everything()
    )
}