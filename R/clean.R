#' Clean an EMERGe dataset after HT Stream preprocessing.
#' 
#' Reads a FASTQ file and returns a standardized tibble.
#' 
#' @param fastq_infile (String) Path to a FASTQ file.
#' @param template Character(1) DNA template sequence.
#' @param flanks Character(4) flank subsequences in expected order,
#'               (f1, f2, f3, f4).
#' @param chunk_size (integer) number of sequences to read at a time from the
#'                   provided fastq.
#' @return A tibble with N10-level summaries: number of reads "n," 
#'         number of edits "k," and other codons "other."
#' @export
#' 
#' @examples 
#' \dontrun{
#' # R255XE488Q. Target is GAA edited to GGA. We take 6-mer flanks around GAA 
#' # and the N10. The full RNA hairpin follows:
#' # GGGAATGGCGCTTCCCGGCA GGAAGT GAA AAGCTG AGGCCGACTTTCTTTCTTTCGTCG GCCTCA 
#' # NNNNNNNNNN TCCTGC CGGGTTGTAGGGTC
#' # Flanks, demarcated above by spaces: GGAAGT AAGCTG GCCTCA TCCTGC
#' template = paste0(
#'   "GGGAATGGCGCTTCCCGGCAGGAAGTGAAAAGCTGAGGCCGACTTTCTTTCTTTCGTCGGCCTCA",
#'   "NNNNNNNNNNTCCTGCCGGGTTGTAGGGTC"
#' )
#' flanks = c("GGAAGT", "AAGCTG", "GCCTCA", "TCCTGC")
#' emergeR::clean(
#'   fastq_infile = "example.fastq", template = template, flanks = flanks
#' )
#' }
#' 
#' @importFrom ShortRead readFastq sread
#' @importFrom Biostrings DNAStringSet extractAt subseq
#' @importFrom BiocGenerics width
#' @importFrom IRanges IRangesList start end
#' @importFrom tibble tibble
#' @importFrom dplyr mutate group_by summarise arrange
#' @importFrom magrittr %>%
clean <- function(
    fastq_infile,
    template,
    flanks,
    chunk_size = 500000
) {
  first_flank  <- list(left = 21L, right = 26L)
  second_flank <- list(left = 30L, right = 35L)
  third_flank  <- list(left = 60L, right = 65L)
  fourth_flank <- list(left = 76L, right = 81L)
  n10_region   <- list(left = 66L, right = 75L)
  edit_region  <- list(left = 27L, right = 29L)
  expected_concat <- paste0(flanks, collapse = "")
  
  size_cdna  <- nchar(template)
  size_left  <- size_cdna - 2
  size_right <- size_cdna + 1
  
  reads <- readFastq(fastq_infile)
  w     <- width(reads)
  reads <- reads[w >= size_left & w <= size_right]
  seqs  <- sread(reads)
  L     <- length(seqs)
  
  f1 <- subseq(seqs, first_flank$left,  first_flank$right)
  f2 <- subseq(seqs, second_flank$left, second_flank$right)
  f3 <- subseq(seqs, third_flank$left,  third_flank$right)
  f4 <- subseq(seqs, fourth_flank$left, fourth_flank$right)
  
  perfect <- (xscat(f1, f2, f3, f4) == expected_concat)
  codon <- as.character(
    subseq(seqs[perfect], edit_region$left, edit_region$right)
  )
  n10    <- as.character(
    subseq(seqs[perfect], n10_region$left,  n10_region$right)
  )
    
  clean <- tibble(n10 = toupper(n10), codon = toupper(codon)) %>%
    group_by(n10) %>%
    summarise(
      n = as.numeric(n()),
      k = as.numeric(sum(codon == "GGA")),
      .groups = "drop"
    ) %>%
    arrange(n10) 

  return(clean)
}


# TODO
#
#' Derives the indices for provided flanks within a DNA template.
#' 
#' @param tpl_chr (string) DNA template used in EMERGe screen, with 'NNNNNNNNNN'
#'                filling the N10 region.
#' @param flanks (list of 4 strings) flanks used to identify successful reads.
#'               See clean() for more information.
#' @return indices for flanks, n10 region, and edit region within the 
#'         template DNA.
#' 
#' @keywords internal
#' @noRd
#'
#' @importFrom Biostrings DNAString DNAStringSet matchPattern xscat
#' @importFrom IRanges IRanges start end
derive_regions <- function(tpl_chr, flanks) {
  stopifnot(length(flanks) == 4L)
  tpl <- DNAString(tpl_chr)
  
  hits <- lapply(flanks, function(f) {
    matchPattern(DNAString(f), tpl, max.mismatch = 0L)
  })
  if (!all(lengths(hits) == 1L))
    stop("Flanks not uniquely found in template.")
  
  starts <- vapply(hits, start, 0L)
  ends   <- vapply(hits, end,   0L)
  
  ord    <- order(starts)
  starts <- starts[ord]; ends <- ends[ord]; flanks <- flanks[ord]
  
  # Compute inter-flank regions
  edit_left  <- ends[1] + 1L
  edit_right <- starts[2] - 1L
  
  n10_left   <- ends[3] + 1L
  n10_right  <- starts[4] - 1L
  
  expected <- DNAString(paste0(flanks, collapse = ""))
  list(
    tpl          = tpl,
    first_flank  = list(left = starts[1], right = ends[1]),
    second_flank = list(left = starts[2], right = ends[2]),
    third_flank  = list(left = starts[3], right = ends[3]),
    fourth_flank = list(left = starts[4], right = ends[4]),
    edit_region  = list(left = ends[1] + 1L,  right = starts[2] - 1L),
    n10_region   = list(left = ends[3] + 1L,  right = starts[4] - 1L),
    expected     = expected
  )
}
