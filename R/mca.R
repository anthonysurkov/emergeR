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
#' @param get_eigs (logical) default FALSE, return eigs in a list with the
#'                 MCA tibble for debug purposes.
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
append_mca <- function(X, nf = 50, top_quantile = 0.999, get_eigs = FALSE) {
  Xfac <- factor_sequences(X = X)
  mca  <- dudi.acm(df = Xfac, nf = nf, scannf = FALSE)
  
  top_cutoff <- quantile(X$postmean, top_quantile)
  
  X_mca <- bind_cols(X, as_tibble(mca$li)) %>%
    filter(postmean >= top_cutoff) %>%
    select(
      n10, n, k,
      postmean, lower_cred, upper_cred,
      mle, lower_ci, upper_ci,
      alpha, beta,
      everything()
    )
  
  if (get_eigs) return(list(X_mca = X_mca, eigens = mca$eig))
  return(X_mca)
}

#' Factorizes n10 sequences into position-wise categorical variables.
#' 
#' Splits each sequence in the 'n10' column into its nucleotide characters,
#' aligns them by position, and returns a data frame with one column per
#' position. Each column is a factor with levels A, G, C, T.
#' 
#' @param X (data frame) with an 'n10' column of DNA sequences.
#' @param L (integer), default 10, length of n10 sequence
#' @return data frame with factored nucleotides from n10 sequences.
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom magrittr %>%
factor_sequences <- function(
    X, L = 10, drop_unused = TRUE, drop_constant = TRUE
) {
  Xfac <- X$n10 %>%
    strsplit("") %>%
    do.call(rbind, .) %>%
    as.data.frame()
  colnames(Xfac) <- paste0("p", 1:L)
  
  Xfac[] <- lapply(Xfac, function(x) factor(x, levels = c("A", "G", "C", "T")))
  
  if (drop_unused) {
    Xfac[] <- lapply(Xfac, droplevels)
  }
  
  if (drop_constant) {
    keep <- vapply(Xfac, nlevels, integer(1)) > 1L
    Xfac <- Xfac[, keep, drop = FALSE]
  }
  
  Xfac
}
